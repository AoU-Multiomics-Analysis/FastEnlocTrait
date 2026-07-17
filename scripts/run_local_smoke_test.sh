#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
  cat >&2 <<'USAGE'
Usage: run_local_smoke_test.sh OPEN_TARGETS_ROOT QTL_DIR OUTPUT_DIR [FASTENLOC] [GTF]

Selects non-MTAG GWAS inputs according to SELECTION_MODE (category_max by
default), runs them against eQTL, sQTL, and pQTL, then follows the workflow's
CLPP, consensus, harmonization, summary, and plotting logic. Existing pair
outputs are reused when possible.

SELECTION_MODE may be category_max, trait_max, or all. The *_max modes keep
the study with the most credible sets in each category or distinct trait.
USAGE
  exit 2
fi

open_targets_root=$1
qtl_dir=$2
outdir=$3
fastenloc_bin=${4:-fastenloc}
gtf=${5:-}
selection_mode=${SELECTION_MODE:-category_max}

case "$selection_mode" in
  category_max|trait_max|all) ;;
  *)
    echo "Invalid SELECTION_MODE: $selection_mode (expected category_max, trait_max, or all)" >&2
    exit 2
    ;;
esac

pipeline_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source_manifest="$open_targets_root/manifests/gwas_manifest.tsv"
relocation="$open_targets_root/manifests/gwas_gcs_relocation.tsv"

qtl_labels=(eQTL sQTL pQTL)
qtl_files=(
  "$qtl_dir/eQTL_all.fastenloc.vcf.gz"
  "$qtl_dir/sQTL_all.fastenloc.vcf.gz"
  "$qtl_dir/pQTL_all.fastenloc.vcf.gz"
)

for path in "$source_manifest" "$relocation" "${qtl_files[@]}"; do
  if [[ ! -f "$path" ]]; then
    echo "Missing required input: $path" >&2
    exit 1
  fi
done
if [[ -n "$gtf" && ! -f "$gtf" ]]; then
  echo "Missing GTF: $gtf" >&2
  exit 1
fi
if [[ "$fastenloc_bin" == */* ]]; then
  [[ -x "$fastenloc_bin" ]] || { echo "fastenloc is not executable: $fastenloc_bin" >&2; exit 1; }
elif ! command -v "$fastenloc_bin" >/dev/null 2>&1; then
  echo "fastenloc is not on PATH: $fastenloc_bin" >&2
  exit 1
fi

mkdir -p "$outdir"/{inputs,work,raw,clpp,harmonized,combined,raw_summaries,final_summaries,figures,logs}
selected_manifest="$outdir/inputs/selected_manifest.tsv"

Rscript - "$source_manifest" "$relocation" "$open_targets_root" "$selected_manifest" "$selection_mode" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
manifest <- fread(args[[1]]) |> as_tibble()
paths <- fread(args[[2]]) |>
  as_tibble() |>
  transmute(study_id, local_source_path)

selection_mode <- args[[5]]
normalize_label <- function(x) {
  tolower(trimws(gsub("[[:space:]]+", " ", x)))
}
ranked <- manifest |>
  mutate(n_credible_sets = as.integer(n_credible_sets),
         n_variants = as.integer(n_variants),
         trait = normalize_label(trait),
         trait_category = normalize_label(trait_category))

mtag_rows <- grepl("MTAG", ranked$trait, ignore.case = TRUE)
if (any(mtag_rows)) {
  message(
    "Excluded ", sum(mtag_rows), " MTAG stud", ifelse(sum(mtag_rows) == 1, "y", "ies"),
    ": ", paste(ranked$study_id[mtag_rows], collapse = ", ")
  )
}
ranked <- ranked |>
  filter(!mtag_rows)

selected <- switch(
  selection_mode,
  category_max = ranked |>
    group_by(trait_category) |>
    slice_max(order_by = n_credible_sets, n = 1, with_ties = FALSE) |>
    ungroup(),
  trait_max = ranked |>
    group_by(trait) |>
    slice_max(order_by = n_credible_sets, n = 1, with_ties = FALSE) |>
    ungroup(),
  all = ranked
) |>
  left_join(paths, by = "study_id") |>
  mutate(gwas_path = file.path(args[[3]], local_source_path)) |>
  select(study_id, trait, n_variants, gwas_path, trait_category, n_credible_sets) |>
  arrange(trait_category, trait, study_id)

stopifnot(nrow(selected) > 0, !anyDuplicated(selected$study_id), all(file.exists(selected$gwas_path)))
write_tsv(selected, args[[4]])
message(
  "Selected ", nrow(selected), " studies representing ", n_distinct(selected$trait),
  " traits across ", n_distinct(selected$trait_category), " categories (", selection_mode, ")"
)
RSCRIPT

consensus="$outdir/combined/consensus_loci.tsv.gz"
Rscript "$pipeline_root/scripts/merge_credible_sets.R" \
  --manifest "$selected_manifest" \
  --out "$consensus" \
  --summary_out "$outdir/combined/consensus_loci.summary.tsv" \
  --jaccard 0.90 \
  --group_col trait \
  >"$outdir/logs/merge_credible_sets.log" 2>&1

printf 'study_id\ttrait\ttrait_category\tn_variants\tn_credible_sets\tlayer\tstatus\n' \
  > "$outdir/run_manifest.tsv"

tail -n +2 "$selected_manifest" | while IFS=$'\t' read -r study_id trait n_variants gwas_path trait_category n_credible_sets; do
  study_key=$(printf '%s' "$study_id" | tr -c 'A-Za-z0-9._-' '_')
  gwas_prepped="$outdir/work/${study_key}.gwas.fastenloc.txt"
  if [[ ! -s "$gwas_prepped" ]]; then
    if [[ "$gwas_path" == *.gz ]]; then
      gzip -dc "$gwas_path"
    else
      cat "$gwas_path"
    fi | awk -F'\t' 'NR == 1 && ($6 == "annotation" || $6 == "locus_string") { next } { print }' \
      > "$gwas_prepped"
  fi

  pids=()
  for qtl_index in "${!qtl_labels[@]}"; do
    {
    label=${qtl_labels[$qtl_index]}
    qtl=${qtl_files[$qtl_index]}
    pair_key="${study_key}.${label}"
    raw_prefix="$outdir/raw/$pair_key"
    clpp_out="$outdir/clpp/$pair_key.clpp.tsv.gz"
    harmonized_prefix="$outdir/harmonized/$pair_key.harmonized_coloc"
    log="$outdir/logs/$pair_key.log"

    if [[ ! -s "$raw_prefix.enloc.gene.out" ]]; then
      "$fastenloc_bin" \
        -eqtl "$qtl" \
        -gwas "$gwas_prepped" \
        -total_variants "$n_variants" \
        -thread 1 \
        -prefix "$raw_prefix" \
        >"$log" 2>&1

      for output in "$raw_prefix".enloc.*.out; do
        header=$(head -n 1 "$output")
        {
          printf 'trait\t%s\n' "$header"
          tail -n +2 "$output" | awk -v value="$trait" 'BEGIN{OFS="\t"}{print value,$0}'
        } > "$output.tmp"
        mv "$output.tmp" "$output"
      done
    fi

    if [[ ! -s "$clpp_out" ]]; then
      Rscript "$pipeline_root/scripts/clpp_fastenloc.R" \
        --gwas "$gwas_path" \
        --qtl "$qtl" \
        --out "$clpp_out" \
        --min_clpp 0.01 \
        >>"$log" 2>&1
    fi

    if [[ ! -s "$harmonized_prefix.gene.tsv.gz" ]]; then
      Rscript "$pipeline_root/scripts/harmonize_coloc.R" \
        --sig "$raw_prefix.enloc.sig.out" \
        --gene "$raw_prefix.enloc.gene.out" \
        --clpp "$clpp_out" \
        --gwas "$gwas_path" \
        --consensus_map "$consensus" \
        --out "$harmonized_prefix.signal.tsv.gz" \
        --cs_out "$harmonized_prefix.cs.tsv.gz" \
        --gene_out "$harmonized_prefix.gene.tsv.gz" \
        --fdr_level 0.05 \
        --study "$study_id" \
        --trait "$trait" \
        --trait_category "$trait_category" \
        --n_variants "$n_variants" \
        --n_credible_sets "$n_credible_sets" \
        --layer "$label" \
        >>"$log" 2>&1
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\tcomplete\n' \
      "$study_id" "$trait" "$trait_category" "$n_variants" "$n_credible_sets" "$label" \
      >> "$outdir/run_manifest.tsv"
    echo "Completed $study_id x $label"
    } &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
done

aggregate_raw() {
  local family=$1
  local suffix=$2
  local output=$3
  local compressed=$4
  local current="$output.current"
  local wrote_header=0
  : > "$output"
  while IFS=$'\t' read -r study_id trait trait_category n_variants n_credible_sets layer status; do
    study_key=$(printf '%s' "$study_id" | tr -c 'A-Za-z0-9._-' '_')
    file="$outdir/$family/${study_key}.${layer}${suffix}"
    if [[ "$compressed" -eq 1 ]]; then
      gzip -dc "$file" > "$current"
    else
      cp "$file" "$current"
    fi
    if [[ $wrote_header -eq 0 ]]; then
      {
        printf 'study_id\tgwas_trait\ttrait_category\tn_variants\tn_credible_sets\tqtl_label\t'
        head -n 1 "$current"
      } > "$output"
      wrote_header=1
    fi
    tail -n +2 "$current" | awk \
      -v study="$study_id" -v trait="$trait" -v category="$trait_category" \
      -v variants="$n_variants" -v credible_sets="$n_credible_sets" -v layer="$layer" \
      'BEGIN{OFS="\t"}{print study,trait,category,variants,credible_sets,layer,$0}' \
      >> "$output"
  done < <(tail -n +2 "$outdir/run_manifest.tsv")
  rm -f "$current"
}

aggregate_raw raw '.enloc.gene.out' "$outdir/combined/combined.enloc.gene.out" 0
aggregate_raw raw '.enloc.enrich.out' "$outdir/combined/combined.enloc.enrich.out" 0
aggregate_raw raw '.enloc.mi.out' "$outdir/combined/combined.enloc.mi.out" 0
aggregate_raw raw '.enloc.sig.out' "$outdir/combined/combined.enloc.sig.out" 0
aggregate_raw raw '.enloc.snp.out' "$outdir/combined/combined.enloc.snp.out" 0
aggregate_raw clpp '.clpp.tsv.gz' "$outdir/combined/clpp.combined.tsv" 1

printf '%s\n' "${qtl_labels[@]}" > "$outdir/inputs/qtl_labels.txt"
Rscript "$pipeline_root/scripts/summarize_raw_coloc.R" \
  --gene "$outdir/combined/combined.enloc.gene.out" \
  --enrich "$outdir/combined/combined.enloc.enrich.out" \
  --mi "$outdir/combined/combined.enloc.mi.out" \
  --sig "$outdir/combined/combined.enloc.sig.out" \
  --snp "$outdir/combined/combined.enloc.snp.out" \
  --clpp "$outdir/combined/clpp.combined.tsv" \
  --qtl_labels "$outdir/inputs/qtl_labels.txt" \
  --out "$outdir/raw_summaries/raw_coloc_summary.all_qtl.tsv" \
  --per_qtl_dir "$outdir/raw_summaries/per_qtl" \
  >"$outdir/logs/summarize_raw_coloc.log" 2>&1

aggregate_gz_tsv() {
  local pattern=$1
  local output=$2
  local tmp="$output.plain"
  local current="$output.current"
  local wrote_header=0
  : > "$tmp"
  while IFS= read -r file; do
    gzip -dc "$file" > "$current"
    if [[ $wrote_header -eq 0 ]]; then
      head -n 1 "$current" > "$tmp"
      wrote_header=1
    fi
    tail -n +2 "$current" >> "$tmp"
  done < <(find "$outdir/harmonized" -type f -name "$pattern" | sort)
  [[ $wrote_header -eq 1 ]] || { echo "No files matched $pattern" >&2; exit 1; }
  gzip -c "$tmp" > "$output"
  rm -f "$tmp" "$current"
}

aggregate_gz_tsv '*.signal.tsv.gz' "$outdir/combined/harmonized_coloc.signal.tsv.gz"
aggregate_gz_tsv '*.cs.tsv.gz' "$outdir/combined/harmonized_coloc.cs.tsv.gz"
aggregate_gz_tsv '*.gene.tsv.gz' "$outdir/combined/harmonized_coloc.gene.tsv.gz"

gtf_args=()
if [[ -n "$gtf" ]]; then
  gtf_args=(--gtf "$gtf")
fi
Rscript "$pipeline_root/scripts/summarize_coloc.R" \
  --cs "$outdir/combined/harmonized_coloc.cs.tsv.gz" \
  --gene "$outdir/combined/harmonized_coloc.gene.tsv.gz" \
  "${gtf_args[@]}" \
  --coloc_rate_out "$outdir/final_summaries/coloc_rate_by_trait.tsv" \
  --gene_summary_out "$outdir/final_summaries/gene_summary_by_trait.tsv" \
  --gene_list_out "$outdir/final_summaries/colocalizing_genes_long.tsv" \
  --gene_threshold any \
  >"$outdir/logs/summarize_coloc.log" 2>&1

for min_coloc_genes in 1 5 10; do
  suffix=".min${min_coloc_genes}"
  if [[ "$min_coloc_genes" -eq 1 ]]; then
    suffix=""
  fi
  Rscript "$pipeline_root/scripts/plot_coloc_summary.R" \
    --gene "$outdir/combined/harmonized_coloc.gene.tsv.gz" \
    --coloc_rate "$outdir/final_summaries/coloc_rate_by_trait.tsv" \
    --plot_data_out "$outdir/final_summaries/coloc_summary_plot_data${suffix}.tsv" \
    --png_out "$outdir/figures/coloc_summary${suffix}.png" \
    --pdf_out "$outdir/figures/coloc_summary${suffix}.pdf" \
    --min_coloc_genes "$min_coloc_genes" \
    >>"$outdir/logs/plot_coloc_summary.log" 2>&1
done

Rscript - "$outdir" <<'RSCRIPT'
suppressPackageStartupMessages({library(data.table); library(dplyr); library(readr)})
root <- commandArgs(trailingOnly = TRUE)[[1]]
run <- fread(file.path(root, "run_manifest.tsv"))
cs <- fread(cmd = paste("gzip -dc", shQuote(file.path(root, "combined/harmonized_coloc.cs.tsv.gz"))))
gene <- fread(cmd = paste("gzip -dc", shQuote(file.path(root, "combined/harmonized_coloc.gene.tsv.gz"))))
rates <- fread(file.path(root, "final_summaries/coloc_rate_by_trait.tsv"))
plot_data <- fread(file.path(root, "final_summaries/coloc_summary_plot_data.tsv"))
qc <- tibble(
  metric = c("studies", "trait_categories", "qtl_layers", "completed_pairs",
             "harmonized_cs_rows", "harmonized_gene_rows", "union_rate_rows", "plotted_traits"),
  value = c(n_distinct(run$study_id), n_distinct(run$trait_category), n_distinct(run$layer),
            sum(run$status == "complete"), nrow(cs), nrow(gene), sum(rates$layer == "union"), nrow(plot_data))
)
write_tsv(qc, file.path(root, "QC_SUMMARY.tsv"))
print(qc)
RSCRIPT

echo "Smoke test complete: $outdir"
