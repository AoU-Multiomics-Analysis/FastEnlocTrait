# Helper Scripts

The workflow image copies all files under `scripts/` into `/home/mambauser`, so WDL tasks call them as `~/script_name.R`.

## Split Trait Data

`scripts/SplitTraitData.R` is retained as an optional utility for older multi-trait fastENLOC trait files. The main `RunFastenloc` workflow does not call this script; current manifest rows are expected to point to one trait/study analysis unit.

```bash
Rscript scripts/SplitTraitData.R \
  --input traits.fastenloc.txt \
  --traits-per-chunk 25
```

## Prepare QTL Fine-Mapping Data

`scripts/PrepQTLFinemapping.R` converts fine-mapped QTL data into fastENLOC QTL format.

```bash
Rscript scripts/PrepQTLFinemapping.R \
  --QTLData qtl_finemapping.tsv \
  --QTLType Expression \
  --OutputFile qtl.fastenloc.vcf.gz
```

Supported `--QTLType` values:

- `Expression`
- `Splicing`
- `Protein`

## Prepare MVP Fine-Mapping Data

`scripts/PrepMVPFineMapping.R` converts an MVP fine-mapping Excel file into fastENLOC trait format.

```bash
Rscript scripts/PrepMVPFineMapping.R \
  --TraitData mvp_finemapping.xlsx \
  --OutputFile MVP.all.fastenloc.vcf.gz
```

The current script filters the MVP input to `Trait == "WBC_Mean_INT"` before writing output.

## Compute CLPP

`scripts/clpp_fastenloc.R` computes CLPP directly from fastENLOC-format GWAS/trait and QTL files:

```bash
Rscript scripts/clpp_fastenloc.R \
  --gwas MVP.all.fastenloc.vcf.gz \
  --qtl qtl.fastenloc.vcf.gz \
  --out clpp_pairs.tsv.gz \
  --min_clpp 0.01
```

The WDL runs this script on each localized GWAS analysis unit and QTL input, then aggregates the resulting TSV files per GWAS, per QTL label, and across all manifest rows.

## Merge GWAS Credible Sets

`scripts/merge_credible_sets.R` merges GWAS credible sets across studies of the same trait when their variant-membership Jaccard index meets the threshold:

```bash
Rscript scripts/merge_credible_sets.R \
  --manifest localized_gwas_manifest.tsv \
  --out consensus_loci.tsv.gz \
  --summary_out consensus_loci.summary.tsv \
  --jaccard 0.90 \
  --group_col trait
```

The manifest passed to this script needs `study_id`, `trait`, and `gwas_path`. In the WDL, this manifest is created after GWAS localization, so `gwas_path` points at job-local files. The output maps each original `study_id`/`gwas_cs` pair to a `consensus_locus_id`.

## Harmonize fastENLOC and CLPP

`scripts/harmonize_coloc.R` joins combined fastENLOC signal and gene outputs with the combined CLPP pairs file:

```bash
Rscript scripts/harmonize_coloc.R \
  --sig combined.enloc.sig.out \
  --gene combined.enloc.gene.out \
  --clpp clpp.combined.tsv \
  --gwas MVP.all.fastenloc.vcf.gz \
  --consensus_map consensus_loci.tsv.gz \
  --out harmonized_coloc.signal.tsv.gz \
  --cs_out harmonized_coloc.cs.tsv.gz \
  --gene_out harmonized_coloc.gene.tsv.gz \
  --fdr_level 0.05 \
  --study GCST90027158 \
  --trait "Alzheimer disease" \
  --trait_category neuro \
  --n_variants 1000000 \
  --n_credible_sets 245 \
  --layer eQTL
```

The WDL runs this after per-GWAS/QTL aggregation and emits signal-level, credible-set-level, and gene-level harmonized tables. Metadata flags are stamped into every harmonized output row, and credible-set rollups carry `consensus_locus_id` for de-duplicated trait coverage. The gene-level output keeps every fastENLOC gene-level row, then joins best signal-level RCP and CLPP evidence when a matching signal is available.

## Summarize Raw fastENLOC and CLPP Outputs

`scripts/summarize_raw_coloc.R` summarizes the all-GWAS/all-QTL raw fastENLOC and CLPP combined outputs by QTL label:

```bash
Rscript scripts/summarize_raw_coloc.R \
  --gene combined.enloc.gene.out \
  --enrich combined.enloc.enrich.out \
  --mi combined.enloc.mi.out \
  --sig combined.enloc.sig.out \
  --snp combined.enloc.snp.out \
  --clpp clpp.combined.tsv \
  --qtl_labels qtl_labels.txt \
  --out raw_coloc_summary.all_qtl.tsv \
  --per_qtl_dir raw_coloc_summary_by_qtl
```

`qtl_labels.txt` should contain one validated QTL label per line. The all-QTL output has one row per `qtl_label` and raw output family (`fastenloc_gene`, `fastenloc_enrich`, `fastenloc_mi`, `fastenloc_sig`, `fastenloc_snp`, and `clpp`). The per-QTL directory contains one compact TSV per QTL label.

## Summarize Colocalization

`scripts/summarize_coloc.R` summarizes the all-GWAS/all-QTL harmonized outputs after consensus-locus harmonization:

```bash
Rscript scripts/summarize_coloc.R \
  --cs harmonized_coloc.cs.tsv.gz \
  --gene harmonized_coloc.gene.tsv.gz \
  --gtf gencode.v44.annotation.gtf.gz \
  --coloc_rate_out coloc_rate_by_trait.tsv \
  --gene_summary_out gene_summary_by_trait.tsv \
  --gene_list_out colocalizing_genes_long.tsv \
  --gene_threshold any
```

The WDL runs this once at the end of the workflow. It produces trait x QTL-layer summaries plus a `union` layer that counts a consensus locus or gene once if it colocalizes in any QTL layer. Coverage rates use `consensus_locus_id` as the denominator, so cross-study GWAS credible sets for the same trait are not double counted. Gene summaries from this script are derived from the credible-set-level gene lists; use `harmonized_coloc.gene.tsv.gz` when you want the complete fastENLOC gene-level table, including high-GLCP/GRCP genes without matched signal-level evidence.

## Plot the Final Trait Summary

`scripts/plot_coloc_summary.R` reproduces the paired bar/lollipop summary at the end of the WDL:

```bash
Rscript scripts/plot_coloc_summary.R \
  --gene harmonized_coloc.gene.tsv.gz \
  --coloc_rate coloc_rate_by_trait.tsv \
  --plot_data_out coloc_summary_plot_data.tsv \
  --png_out coloc_summary.png \
  --pdf_out coloc_summary.pdf
```

The left panel counts distinct `(gene, trait, trait_category)` rows with gene-level `any_coloc == TRUE`, so repeated evidence across studies or QTL layers does not inflate the count. The right panel uses `pct_stringent_union` from the `layer == "union"` row of `coloc_rate_by_trait.tsv`; that rate is based on de-duplicated consensus loci. Traits without a union rate or without any colocalizing gene are not plotted. The script writes the joined plotting data to TSV for auditability and assigns fallback colors to categories not present in its built-in palette.
