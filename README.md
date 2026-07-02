# FastEnlocTrait

Workflow and helper scripts for running [fastENLOC](https://github.com/xqwen/fastenloc) colocalization across multiple GWAS traits and QTL fine-mapping datasets. The main workflow splits a multi-trait GWAS/trait input into smaller chunks, runs `fastenloc` for each trait against one or more QTL inputs, computes CLPP colocalization directly from the same fastENLOC input files, aggregates the standard `.enloc.*.out` and CLPP result files, and harmonizes fastENLOC and CLPP metrics into tidy downstream tables.

## Repository layout

```text
.
├── envs/
│   └── Dockerfile
├── scripts/
│   ├── PrepMVPFineMapping.R
│   ├── PrepQTLFinemapping.R
│   ├── SplitTraitData.R
│   ├── clpp_fastenloc.R
│   └── harmonize_coloc.R
└── workflows/
    └── RunFastEnlocTraits.wdl
```

## Main workflow

The WDL entry point is:

```text
workflows/RunFastEnlocTraits.wdl
```

Workflow name:

```text
RunFastenloc
```

### Inputs

| Input | Type | Description |
| --- | --- | --- |
| `FastEnlocTraitData` | `File` | Trait/GWAS fastENLOC input. Expected to be tab-delimited with six columns: chromosome, position, variant ID, reference allele, alternate allele, and an annotation/locus string. |
| `QTLData` | `Array[File]` | One or more QTL fastENLOC input files, each formatted for `fastenloc -eqtl`. |
| `QTLLabels` | `Array[String]` | Label for each QTL input, such as `eQTL`, `sQTL`, or `pQTL`. Must have the same length as `QTLData`, be unique, and match `[A-Za-z0-9._-]+`. |
| `NumberVariants` | `Int` | Total number of variants passed to `fastenloc -total_variants`. |
| `min_clpp` | `Float` | Minimum CLPP value to report. Defaults to `0.01`. |
| `clpp_output_prefix` | `String` | Prefix for CLPP per-chunk and combined output files. Defaults to `clpp`. |
| `harmonized_fdr_level` | `Float` | Bayesian FDR level used for harmonized pass flags. Defaults to `0.05`. |
| `harmonized_output_prefix` | `String` | Prefix for harmonized signal, credible-set, and gene outputs. Defaults to `harmonized_coloc`. |

`SplitFastenloc` also supports `traits_per_chunk`, which defaults to `25` inside the task.

For a single QTL layer, provide one-element `QTLData` and `QTLLabels` arrays.

### Outputs

The workflow aggregates per-trait fastENLOC outputs into:

| Output | Description |
| --- | --- |
| `combined_gene_out` | Combined `*.enloc.gene.out` results. |
| `combined_enrich_out` | Combined `*.enloc.enrich.out` results. |
| `combined_mi_out` | Combined `*.enloc.mi.out` results. |
| `combined_sig_out` | Combined `*.enloc.sig.out` results. |
| `combined_snp_out` | Combined `*.enloc.snp.out` results. |
| `combined_clpp_out` | Combined CLPP results from `clpp_fastenloc.R`. Default filename is `clpp.combined.tsv`. |
| `harmonized_signal_out` | Signal-level table joining fastENLOC RCP/LCP, CLPP, and gene-level GRCP/GLCP. |
| `harmonized_credible_set_out` | Credible-set-level rollup with colocalization flags and per-method gene lists. |
| `harmonized_gene_out` | Gene-level rollup with best signal metrics and gene-native fastENLOC metrics. |

The all-QTL raw combined fastENLOC and CLPP outputs include a leading `qtl_label` column. Each per-trait fastENLOC output is also annotated with an added leading `trait` column before aggregation. Harmonized outputs store the QTL label in the existing `layer` column.

The workflow also emits per-QTL output arrays for debugging and downstream layer-specific analysis:

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

## Docker image

The workflow runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

The Dockerfile installs R dependencies, builds `fastenloc` from `xqwen/fastenloc`, copies the scripts into `/home/mambauser`, and makes the `fastenloc` binary available on `PATH`.

To build locally:

```bash
docker build -t fastenloctrait:local -f envs/Dockerfile .
```

## Helper scripts

### Split trait data

`scripts/SplitTraitData.R` splits a multi-trait fastENLOC trait file into chunk files under `chunks/` and writes `chunk_manifest.txt`.

```bash
Rscript scripts/SplitTraitData.R \
  --input traits.fastenloc.txt \
  --traits-per-chunk 25
```

### Prepare QTL fine-mapping data

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

### Compute CLPP colocalization

`scripts/clpp_fastenloc.R` computes CLPP directly from fastENLOC-format GWAS/trait and QTL files:

```bash
Rscript scripts/clpp_fastenloc.R \
  --gwas MVP.all.fastenloc.vcf.gz \
  --qtl qtl.fastenloc.vcf.gz \
  --out clpp_pairs.tsv.gz \
  --min_clpp 0.01
```

The WDL runs this script on each split trait chunk and QTL input, then aggregates the resulting plain TSV files per QTL label and across all QTL labels. With the default `clpp_output_prefix`, per-layer combined outputs are named like `eQTL.clpp.combined.tsv`, and the all-layer combined output is `clpp.combined.tsv`.

### Harmonize fastENLOC and CLPP outputs

`scripts/harmonize_coloc.R` joins combined fastENLOC signal and gene outputs with the combined CLPP pairs file:

```bash
Rscript scripts/harmonize_coloc.R \
  --sig combined.enloc.sig.out \
  --gene combined.enloc.gene.out \
  --clpp clpp.combined.tsv \
  --gwas MVP.all.fastenloc.vcf.gz \
  --out harmonized_coloc.signal.tsv.gz \
  --cs_out harmonized_coloc.cs.tsv.gz \
  --gene_out harmonized_coloc.gene.tsv.gz \
  --fdr_level 0.05
```

The WDL runs this after per-QTL aggregation and emits signal-level, credible-set-level, and gene-level harmonized tables. The QTL label is passed as `--layer`, so multi-QTL harmonized outputs can be grouped by QTL layer.

### Prepare MVP fine-mapping data

`scripts/PrepMVPFineMapping.R` converts an MVP fine-mapping Excel file into fastENLOC trait format.

```bash
Rscript scripts/PrepMVPFineMapping.R \
  --TraitData mvp_finemapping.xlsx \
  --OutputFile MVP.all.fastenloc.vcf.gz
```

The current script filters the MVP input to `Trait == "WBC_Mean_INT"` before writing output.

## Example WDL inputs

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

## Output format notes

The core fastENLOC files are expected to use six tab-delimited columns:

```text
chromosome  position  variant_id  ref  alt  annotation
```

For trait data, the workflow splits the sixth column on `;` and uses the first field as the trait name. For QTL data, annotations can contain multiple records joined by `|`.

## Validation

The R scripts can be syntax-checked with:

```bash
Rscript -e 'for (f in list.files("scripts", pattern = "[.]R$", full.names = TRUE)) { parse(f); cat("parse OK:", f, "\n") }'
```

This repository does not currently include small example inputs or automated tests, so full end-to-end validation requires representative trait and QTL fine-mapping files plus a WDL runner such as Cromwell, miniwdl, or Dockstore.
