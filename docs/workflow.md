# Workflow Reference

The main WDL entry point is:

```text
workflows/RunFastEnlocTraits.wdl
```

Workflow name:

```text
RunFastenloc
```

The primary WDL imports task modules from `workflows/tasks/`:

| Module | Tasks |
| --- | --- |
| `aggregation.wdl` | `AggregateFiles`, `AggregateFilesWithQTLLabel`, `AggregateFilesWithGWASMetadata`, `AggregateGzTsvFiles` |
| `clpp.wdl` | `CLPPFastEnloc` |
| `consensus.wdl` | `MergeCredibleSets` |
| `fastenloc.wdl` | `FastEnloc` |
| `harmonize.wdl` | `HarmonizeColoc` |
| `input_validation.wdl` | `ValidateGWASManifest`, `ValidateQTLInputs` |
| `localize.wdl` | `LocalizeGWASData` |

`workflows/tasks/split.wdl` is retained as an optional legacy task for pre-splitting multi-trait inputs, but it is not imported by the main workflow.

## Inputs

| Input | Type | Description |
| --- | --- | --- |
| `GWASManifest` | `File` | Tab-delimited manifest describing one or more GWAS analysis units. |
| `QTLData` | `Array[File]` | One or more QTL fastENLOC input files, each formatted for `fastenloc -eqtl`. |
| `QTLLabels` | `Array[String]` | Label for each QTL input, such as `eQTL`, `sQTL`, or `pQTL`. Must have the same length as `QTLData`, be unique, and match `[A-Za-z0-9._-]+`. |
| `min_clpp` | `Float` | Minimum CLPP value to report. Defaults to `0.01`. |
| `clpp_output_prefix` | `String` | Prefix for CLPP output files. Defaults to `clpp`. |
| `consensus_jaccard` | `Float` | Jaccard threshold for merging cross-study GWAS credible sets within a trait. Defaults to `0.90`. |
| `consensus_output_prefix` | `String` | Prefix for consensus locus map and summary outputs. Defaults to `consensus_loci`. |
| `harmonized_fdr_level` | `Float` | Bayesian FDR level used for harmonized pass flags. Defaults to `0.05`. |
| `harmonized_output_prefix` | `String` | Prefix for harmonized signal, credible-set, and gene outputs. Defaults to `harmonized_coloc`. |

## GWAS Manifest

The manifest must contain these columns:

| Column | Type | Description |
| --- | --- | --- |
| `study_id` | `String` | Primary key for the GWAS analysis unit. Must be unique and match `[A-Za-z0-9._-]+`. A leading `#study_id` header is also accepted. |
| `trait` | `String` | Human-readable trait label stamped into harmonized outputs. |
| `n_variants` | `Int` | GWAS variant denominator passed to `fastenloc -total_variants`. |
| `gwas_path` | `String` | URI or in-runtime path to the fastENLOC-format GWAS file. Supports `gs://`, HTTP(S), and local paths visible inside the task. |
| `trait_category` | `String` | One of `immune`, `neuro`, `cardiometabolic`, `cancer`, or `anthropometric`. |
| `n_credible_sets` | `Int` | Total fine-mapped GWAS credible-set denominator for this analysis unit. |

Each manifest row is localized and analyzed independently, so different studies can use different `n_variants`, trait labels, and disease categories. Each row should point to one trait/study analysis unit.

## Example Inputs

```json
{
  "RunFastenloc.GWASManifest": "gwas_manifest.tsv",
  "RunFastenloc.QTLData": [
    "eqtl.fastenloc.vcf.gz",
    "sqtl.fastenloc.vcf.gz",
    "pqtl.fastenloc.vcf.gz"
  ],
  "RunFastenloc.QTLLabels": ["eQTL", "sQTL", "pQTL"],
  "RunFastenloc.min_clpp": 0.01,
  "RunFastenloc.clpp_output_prefix": "clpp",
  "RunFastenloc.consensus_jaccard": 0.90,
  "RunFastenloc.consensus_output_prefix": "consensus_loci",
  "RunFastenloc.harmonized_fdr_level": 0.05,
  "RunFastenloc.harmonized_output_prefix": "harmonized_coloc"
}
```

## Workflow Behavior

The workflow first validates the GWAS manifest and QTL labels. It then runs a two-stage GWAS flow:

1. `LocalizeGWASData` materializes the row's `gwas_path` as a job-local gzip file.
2. `MergeCredibleSets` builds a localized manifest from those files and merges credible sets within each manifest `trait`.
3. The workflow scatters over `QTLData`/`QTLLabels`.
4. For each GWAS/QTL pair, fastENLOC and CLPP run directly on the localized GWAS file.
5. Per-GWAS/QTL outputs are aggregated, then harmonized with the consensus map.
6. Per-GWAS all-QTL outputs are aggregated.
7. Global all-GWAS/all-QTL outputs are aggregated.

`MergeCredibleSets` only merges cross-study credible sets within the same trait when their variant-membership Jaccard index is at least `consensus_jaccard`. Same-study credible sets are kept distinct.

Raw per-GWAS outputs get a leading `qtl_label` column. Raw global outputs prepend `study_id`, `gwas_trait`, `trait_category`, `n_variants`, and `n_credible_sets` before `qtl_label`. Harmonized outputs carry GWAS metadata directly, keep the manifest trait in `trait`, and keep the QTL label in `layer`. Credible-set harmonized outputs also carry `consensus_locus_id`, `n_studies_in_locus`, `n_cs_in_locus`, and `is_merged`.

## Primary Outputs

| Output | Description |
| --- | --- |
| `normalized_gwas_manifest` | Validated manifest normalized to the required column order. |
| `localized_gwas_manifest` | Internal manifest used by `MergeCredibleSets`, with job-local GWAS file paths. |
| `consensus_loci_out` | Per-credible-set consensus map from `merge_credible_sets.R`. |
| `consensus_loci_summary_out` | Per-trait summary of raw credible sets, consensus loci, and merged-away duplicate sets. |
| `combined_gene_out` | All-GWAS/all-QTL combined `*.enloc.gene.out` results. |
| `combined_enrich_out` | All-GWAS/all-QTL combined `*.enloc.enrich.out` results. |
| `combined_mi_out` | All-GWAS/all-QTL combined `*.enloc.mi.out` results. |
| `combined_sig_out` | All-GWAS/all-QTL combined `*.enloc.sig.out` results. |
| `combined_snp_out` | All-GWAS/all-QTL combined `*.enloc.snp.out` results. |
| `combined_clpp_out` | All-GWAS/all-QTL combined CLPP results. Default filename is `clpp.combined.tsv`. |
| `harmonized_signal_out` | All-GWAS/all-QTL signal-level harmonized table. |
| `harmonized_credible_set_out` | All-GWAS/all-QTL credible-set-level harmonized table. |
| `harmonized_gene_out` | All-GWAS/all-QTL gene-level harmonized table. |

## Per-GWAS Outputs

| Output | Description |
| --- | --- |
| `per_gwas_combined_gene_out` | Per-GWAS all-QTL `*.enloc.gene.out` combined files. |
| `per_gwas_combined_enrich_out` | Per-GWAS all-QTL `*.enloc.enrich.out` combined files. |
| `per_gwas_combined_mi_out` | Per-GWAS all-QTL `*.enloc.mi.out` combined files. |
| `per_gwas_combined_sig_out` | Per-GWAS all-QTL `*.enloc.sig.out` combined files. |
| `per_gwas_combined_snp_out` | Per-GWAS all-QTL `*.enloc.snp.out` combined files. |
| `per_gwas_combined_clpp_out` | Per-GWAS all-QTL CLPP combined files. |
| `per_gwas_harmonized_signal_out` | Per-GWAS all-QTL harmonized signal-level files. |
| `per_gwas_harmonized_credible_set_out` | Per-GWAS all-QTL harmonized credible-set-level files. |
| `per_gwas_harmonized_gene_out` | Per-GWAS all-QTL harmonized gene-level files. |

## Per-GWAS/QTL Outputs

These outputs are nested arrays with shape `Array[GWAS][QTL]`.

| Output | Description |
| --- | --- |
| `per_gwas_qtl_combined_gene_out` | Per-GWAS/QTL `*.enloc.gene.out` combined files. |
| `per_gwas_qtl_combined_enrich_out` | Per-GWAS/QTL `*.enloc.enrich.out` combined files. |
| `per_gwas_qtl_combined_mi_out` | Per-GWAS/QTL `*.enloc.mi.out` combined files. |
| `per_gwas_qtl_combined_sig_out` | Per-GWAS/QTL `*.enloc.sig.out` combined files. |
| `per_gwas_qtl_combined_snp_out` | Per-GWAS/QTL `*.enloc.snp.out` combined files. |
| `per_gwas_qtl_combined_clpp_out` | Per-GWAS/QTL CLPP combined files. |
| `per_gwas_qtl_harmonized_signal_out` | Per-GWAS/QTL harmonized signal-level files. |
| `per_gwas_qtl_harmonized_credible_set_out` | Per-GWAS/QTL harmonized credible-set-level files. |
| `per_gwas_qtl_harmonized_gene_out` | Per-GWAS/QTL harmonized gene-level files. |

## Input Format Notes

The core fastENLOC files are expected to use six tab-delimited columns:

```text
chromosome  position  variant_id  ref  alt  annotation
```

For GWAS data, the sixth column stores the credible-set annotation and should correspond to the manifest row's trait/study analysis unit. For QTL data, annotations can contain multiple records joined by `|`.
