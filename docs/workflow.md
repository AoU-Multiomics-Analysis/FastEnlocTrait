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
| `aggregation.wdl` | `AggregateFiles`, `AggregateFilesWithQTLLabel`, `AggregateGzTsvFiles` |
| `clpp.wdl` | `CLPPFastEnloc` |
| `fastenloc.wdl` | `FastEnloc` |
| `harmonize.wdl` | `HarmonizeColoc` |
| `input_validation.wdl` | `ValidateQTLInputs` |
| `split.wdl` | `SplitFastenloc` |

## Inputs

| Input | Type | Description |
| --- | --- | --- |
| `FastEnlocTraitData` | `File` | Trait/GWAS fastENLOC input. Expected to be tab-delimited with six columns: chromosome, position, variant ID, reference allele, alternate allele, and an annotation/locus string. |
| `QTLData` | `Array[File]` | One or more QTL fastENLOC input files, each formatted for `fastenloc -eqtl`. |
| `QTLLabels` | `Array[String]` | Label for each QTL input, such as `eQTL`, `sQTL`, or `pQTL`. Must have the same length as `QTLData`, be unique, and match `[A-Za-z0-9._-]+`. |
| `NumberVariants` | `Int` | Total number of GWAS variants passed to `fastenloc -total_variants`. |
| `min_clpp` | `Float` | Minimum CLPP value to report. Defaults to `0.01`. |
| `clpp_output_prefix` | `String` | Prefix for CLPP per-chunk and combined output files. Defaults to `clpp`. |
| `harmonized_fdr_level` | `Float` | Bayesian FDR level used for harmonized pass flags. Defaults to `0.05`. |
| `harmonized_output_prefix` | `String` | Prefix for harmonized signal, credible-set, and gene outputs. Defaults to `harmonized_coloc`. |

`NumberVariants` is shared across all QTL layers because it depends on the GWAS/trait variant universe.

`SplitFastenloc` also supports `traits_per_chunk`, which defaults to `25` inside the task.

## Example Inputs

```json
{
  "RunFastenloc.FastEnlocTraitData": "MVP.all.fastenloc.vcf.gz",
  "RunFastenloc.QTLData": [
    "eqtl.fastenloc.vcf.gz",
    "sqtl.fastenloc.vcf.gz",
    "pqtl.fastenloc.vcf.gz"
  ],
  "RunFastenloc.QTLLabels": ["eQTL", "sQTL", "pQTL"],
  "RunFastenloc.NumberVariants": 1000000,
  "RunFastenloc.min_clpp": 0.01,
  "RunFastenloc.clpp_output_prefix": "clpp",
  "RunFastenloc.harmonized_fdr_level": 0.05,
  "RunFastenloc.harmonized_output_prefix": "harmonized_coloc"
}
```

## Workflow Behavior

The workflow splits `FastEnlocTraitData` once, then scatters over `QTLData` and `QTLLabels`. For each QTL layer, it runs fastENLOC and CLPP across all trait chunks, aggregates per-QTL outputs, and harmonizes the per-QTL combined results.

After per-QTL processing, the workflow also builds all-QTL combined files:

- Raw all-QTL fastENLOC and CLPP outputs get a leading `qtl_label` column.
- Harmonized all-QTL outputs keep the QTL label in the existing `layer` column.

## Primary Outputs

| Output | Description |
| --- | --- |
| `combined_gene_out` | All-QTL combined `*.enloc.gene.out` results. |
| `combined_enrich_out` | All-QTL combined `*.enloc.enrich.out` results. |
| `combined_mi_out` | All-QTL combined `*.enloc.mi.out` results. |
| `combined_sig_out` | All-QTL combined `*.enloc.sig.out` results. |
| `combined_snp_out` | All-QTL combined `*.enloc.snp.out` results. |
| `combined_clpp_out` | All-QTL combined CLPP results. Default filename is `clpp.combined.tsv`. |
| `harmonized_signal_out` | All-QTL signal-level harmonized table. |
| `harmonized_credible_set_out` | All-QTL credible-set-level harmonized table. |
| `harmonized_gene_out` | All-QTL gene-level harmonized table. |

## Per-QTL Outputs

| Output | Description |
| --- | --- |
| `per_qtl_combined_gene_out` | Per-QTL `*.enloc.gene.out` combined files. |
| `per_qtl_combined_enrich_out` | Per-QTL `*.enloc.enrich.out` combined files. |
| `per_qtl_combined_mi_out` | Per-QTL `*.enloc.mi.out` combined files. |
| `per_qtl_combined_sig_out` | Per-QTL `*.enloc.sig.out` combined files. |
| `per_qtl_combined_snp_out` | Per-QTL `*.enloc.snp.out` combined files. |
| `per_qtl_combined_clpp_out` | Per-QTL CLPP combined files. |
| `per_qtl_harmonized_signal_out` | Per-QTL harmonized signal-level files. |
| `per_qtl_harmonized_credible_set_out` | Per-QTL harmonized credible-set-level files. |
| `per_qtl_harmonized_gene_out` | Per-QTL harmonized gene-level files. |

## Input Format Notes

The core fastENLOC files are expected to use six tab-delimited columns:

```text
chromosome  position  variant_id  ref  alt  annotation
```

For trait data, the workflow splits the sixth column on `;` and uses the first field as the trait name. For QTL data, annotations can contain multiple records joined by `|`.
