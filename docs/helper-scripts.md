# Helper Scripts

The workflow image copies all files under `scripts/` into `/home/mambauser`, so WDL tasks call them as `~/script_name.R`.

## Split Trait Data

`scripts/SplitTraitData.R` splits a multi-trait fastENLOC trait file into chunk files under `chunks/` and writes `chunk_manifest.txt`.

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

The WDL runs this script on each split GWAS chunk and QTL input, then aggregates the resulting TSV files per GWAS, per QTL label, and across all manifest rows.

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

The WDL runs this after per-GWAS/QTL aggregation and emits signal-level, credible-set-level, and gene-level harmonized tables. Metadata flags are stamped into every harmonized output row, and credible-set rollups carry `consensus_locus_id` for de-duplicated trait coverage.
