# FastEnlocTrait

FastEnlocTrait is a WDL workflow for running colocalization for one GWAS analysis unit against one or more QTL layers, such as eQTL, sQTL, and pQTL. It runs fastENLOC, computes CLPP from the same fastENLOC-format inputs, and harmonizes both metrics into tidy downstream tables.

## What the Workflow Does

1. Splits a GWAS fastENLOC input into manageable chunks.
2. Runs fastENLOC for each trait chunk against each QTL input layer.
3. Computes CLPP for the same GWAS/QTL input pairs.
4. Aggregates per-trait and per-QTL outputs.
5. Harmonizes fastENLOC and CLPP results at signal, credible-set, and gene levels.

The main workflow is:

```text
workflows/RunFastEnlocTraits.wdl
```

Workflow name:

```text
RunFastenloc
```

## Quick Start

Use parallel arrays for QTL files and labels:

```json
{
  "RunFastenloc.GWASData": "MVP.all.fastenloc.vcf.gz",
  "RunFastenloc.QTLData": [
    "eqtl.fastenloc.vcf.gz",
    "sqtl.fastenloc.vcf.gz",
    "pqtl.fastenloc.vcf.gz"
  ],
  "RunFastenloc.QTLLabels": ["eQTL", "sQTL", "pQTL"],
  "RunFastenloc.NumberVariants": 1000000
}
```

For a single QTL layer, provide one-element `QTLData` and `QTLLabels` arrays.

For GWAS sources with different `NumberVariants` values, run this workflow once per GWAS analysis unit. A future manifest layer can scatter over those units.

## Primary Outputs

The workflow emits all-QTL combined outputs plus per-QTL output arrays. The most useful downstream outputs are:

| Output | Description |
| --- | --- |
| `harmonized_signal_out` | Signal-level table joining fastENLOC RCP/LCP, CLPP, and gene-level GRCP/GLCP. |
| `harmonized_credible_set_out` | Credible-set-level rollup with colocalization flags and per-method gene lists. |
| `harmonized_gene_out` | Gene-level rollup with best signal metrics and gene-native fastENLOC metrics. |

Raw combined fastENLOC and CLPP outputs include a leading `qtl_label` column. Harmonized outputs store the QTL label in the existing `layer` column.

## Documentation

- [Workflow reference](docs/workflow.md): WDL inputs, outputs, multi-QTL behavior, and output naming.
- [Helper scripts](docs/helper-scripts.md): R script usage for preparing inputs, CLPP, and harmonization.
- [Runtime and validation](docs/runtime-validation.md): Docker image, local build, and validation commands.

## Repository Layout

```text
.
├── docs/
├── envs/
│   └── Dockerfile
├── scripts/
│   ├── PrepMVPFineMapping.R
│   ├── PrepQTLFinemapping.R
│   ├── SplitTraitData.R
│   ├── clpp_fastenloc.R
│   └── harmonize_coloc.R
└── workflows/
    ├── tasks/
    │   ├── aggregation.wdl
    │   ├── clpp.wdl
    │   ├── fastenloc.wdl
    │   ├── harmonize.wdl
    │   ├── input_validation.wdl
    │   └── split.wdl
    └── RunFastEnlocTraits.wdl
```

## Runtime Image

The WDL runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

After adding or changing files under `scripts/`, rebuild and publish the image so the workflow runtime contains the latest helper scripts.
