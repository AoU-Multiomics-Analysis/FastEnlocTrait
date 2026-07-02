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
| `consensus.wdl` | `MergeCredibleSets` |
| `harmonize.wdl` | `HarmonizeColoc` |
| `input_validation.wdl` | `ValidateGWASManifest`, `ValidateQTLInputs` |
| `localize.wdl` | `LocalizeGWASData` |
| `shard_coloc.wdl` | `CreateGWASShards`, `RunColocShard`, `AggregateHarmonizedByGWAS` |
| `summarize.wdl` | `SummarizeRawColocByQTL`, `SummarizeColoc` |

`workflows/tasks/split.wdl` is retained as an optional legacy task for pre-splitting multi-trait inputs, but it is not imported by the main workflow.

## Inputs

| Input | Type | Description |
| --- | --- | --- |
| `GWASManifest` | `File` | Tab-delimited manifest describing one or more GWAS analysis units. |
| `GTF` | `File` | GENCODE GTF used by the final summary step to annotate protein-coding genes and gene names. |
| `QTLData` | `Array[File]` | One or more QTL fastENLOC input files, each formatted for `fastenloc -eqtl`. |
| `QTLLabels` | `Array[String]` | Label for each QTL input, such as `eQTL`, `sQTL`, or `pQTL`. Must have the same length as `QTLData`, be unique, and match `[A-Za-z0-9._-]+`. |
| `min_clpp` | `Float` | Minimum CLPP value to report. Defaults to `0.01`. |
| `clpp_output_prefix` | `String` | Prefix for CLPP output files. Defaults to `clpp`. |
| `consensus_jaccard` | `Float` | Jaccard threshold for merging cross-study GWAS credible sets within a trait. Defaults to `0.90`. |
| `consensus_output_prefix` | `String` | Prefix for consensus locus map and summary outputs. Defaults to `consensus_loci`. |
| `harmonized_fdr_level` | `Float` | Bayesian FDR level used for harmonized pass flags. Defaults to `0.05`. |
| `harmonized_output_prefix` | `String` | Prefix for harmonized signal, credible-set, and gene outputs. Defaults to `harmonized_coloc`. |
| `summary_gene_threshold` | `String` | Threshold used to define colocalizing genes in the final gene summary. Defaults to `any`; accepted values are `any`, `GLCP_FDR`, `RCP_FDR`, `RCP_0.5`, `CLPP_0.01`, and `CLPP_0.05`. |
| `raw_coloc_summary_prefix` | `String` | Prefix for the all-QTL raw fastENLOC/CLPP summary output. Defaults to `raw_coloc_summary`. |
| `coloc_rate_output_name` | `String` | Filename for the final trait/layer colocalization-rate summary. Defaults to `coloc_rate_by_trait.tsv`. |
| `gene_summary_output_name` | `String` | Filename for the final trait/layer gene summary. Defaults to `gene_summary_by_trait.tsv`. |
| `gene_list_output_name` | `String` | Filename for the final long colocalizing gene table. Defaults to `colocalizing_genes_long.tsv`. |
| `gwas_units_per_shard` | `Int` | Number of GWAS manifest rows to process in each raw fastENLOC/CLPP shard job. Defaults to `10`. |

## GWAS Manifest

The manifest must contain these columns:

| Column | Type | Description |
| --- | --- | --- |
| `study_id` | `String` | Primary key for the GWAS analysis unit. Must be unique and match `[A-Za-z0-9._-]+`. A leading `#study_id` header is also accepted. |
| `trait` | `String` | Human-readable trait label stamped into harmonized outputs. |
| `n_variants` | `Int` | GWAS variant denominator passed to `fastenloc -total_variants`. |
| `gwas_path` | `String` | URI or in-runtime path to the fastENLOC-format GWAS file. Supports `gs://`, HTTP(S), and local paths visible inside the task. |
| `trait_category` | `String` | Non-empty user-defined grouping label for downstream summaries and figures. |
| `n_credible_sets` | `Int` | Total fine-mapped GWAS credible-set denominator for this analysis unit. |

Each manifest row is localized and analyzed independently, so different studies can use different `n_variants`, trait labels, and disease categories. Each row should point to one trait/study analysis unit.

## Example Inputs

```json
{
  "RunFastenloc.GWASManifest": "gwas_manifest.tsv",
  "RunFastenloc.GTF": "gencode.v44.annotation.gtf.gz",
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
  "RunFastenloc.harmonized_output_prefix": "harmonized_coloc",
  "RunFastenloc.summary_gene_threshold": "any",
  "RunFastenloc.raw_coloc_summary_prefix": "raw_coloc_summary",
  "RunFastenloc.gwas_units_per_shard": 10
}
```

## Workflow Behavior

The workflow first validates the GWAS manifest and QTL labels. It then runs a two-stage GWAS flow:

1. `LocalizeGWASData` materializes the row's `gwas_path` as a job-local gzip file.
2. `MergeCredibleSets` builds a localized manifest from those files and merges credible sets within each manifest `trait`.
3. `CreateGWASShards` groups manifest rows into shard index files using `gwas_units_per_shard`.
4. The workflow scatters over shards; each `RunColocShard` job runs fastENLOC and CLPP for several GWAS rows and every QTL layer, emitting raw outputs separated by `study_id` and `qtl_label`.
5. The workflow scatters over the resulting flat GWAS x QTL raw outputs and runs `HarmonizeColoc` once per pair with the consensus map.
6. Per-GWAS all-QTL harmonized outputs are aggregated.
7. Global all-GWAS/all-QTL outputs are aggregated.
8. `SummarizeRawColocByQTL` consumes the global raw fastENLOC and CLPP outputs to produce one summary table across QTL labels plus one summary file per QTL label.
9. `SummarizeColoc` consumes the global harmonized credible-set and gene outputs plus `GTF` to produce trait x layer and cross-layer union summaries.

The shard step reduces scheduler overhead for fastENLOC/CLPP while preserving the statistical boundary for harmonization. `HarmonizeColoc` still sees one GWAS analysis unit and one QTL layer per invocation, so Bayesian FDR thresholds are not pooled across studies, traits, or QTL layers.

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
| `raw_coloc_summary_out` | All-QTL summary of raw fastENLOC gene, enrich, MI, signal, SNP, and CLPP outputs grouped by `qtl_label`. |
| `per_qtl_raw_coloc_summary_out` | One raw output summary TSV per QTL label. Each file has one row per raw output family. |
| `harmonized_signal_out` | All-GWAS/all-QTL signal-level harmonized table. |
| `harmonized_credible_set_out` | All-GWAS/all-QTL credible-set-level harmonized table. |
| `harmonized_gene_out` | All-GWAS/all-QTL gene-level harmonized table. |
| `coloc_rate_by_trait_out` | Trait x layer and cross-layer union colocalization rates across de-duplicated consensus loci. |
| `gene_summary_by_trait_out` | Trait x layer and cross-layer union gene counts, genes per locus, protein-coding counts, and gene IDs. |
| `colocalizing_genes_long_out` | Long table of colocalizing genes by trait, layer, consensus locus, GTF gene type, and protein-coding status. |

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

These outputs are flat arrays, with one element per GWAS x QTL pair. Use `per_gwas_qtl_gwas_index`, `per_gwas_qtl_qtl_index`, `per_gwas_qtl_study_id`, and `per_gwas_qtl_qtl_label` to map each file back to the manifest row and QTL layer.

| Output | Description |
| --- | --- |
| `per_gwas_qtl_gwas_index` | Zero-based GWAS manifest row index for each flat GWAS/QTL output. |
| `per_gwas_qtl_qtl_index` | Zero-based QTL array index for each flat GWAS/QTL output. |
| `per_gwas_qtl_study_id` | Manifest `study_id` for each flat GWAS/QTL output. |
| `per_gwas_qtl_qtl_label` | QTL label for each flat GWAS/QTL output. |
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
