# FastEnlocTrait

FastEnlocTrait is a WDL workflow for running colocalization across one or more GWAS analysis units and one or more QTL layers, such as eQTL, sQTL, and pQTL. It runs fastENLOC, computes CLPP from the same fastENLOC-format inputs, merges GWAS credible sets into trait-level consensus loci, and harmonizes the coloc metrics into tidy downstream tables.

## What the Workflow Does

1. Validates a GWAS manifest and QTL layer labels.
2. Localizes each manifest `gwas_path` inside the job.
3. Merges GWAS credible sets across studies of the same trait into consensus loci.
4. Runs fastENLOC and CLPP for each GWAS/QTL pair.
5. Aggregates per-GWAS, per-QTL, and all-GWAS outputs.
6. Harmonizes fastENLOC and CLPP results at signal, credible-set, and gene levels, with consensus locus IDs on credible-set rollups.

The main workflow is:

```text
workflows/RunFastEnlocTraits.wdl
```

Workflow name:

```text
RunFastenloc
```

## Quick Start

Provide a GWAS manifest plus parallel arrays for QTL files and labels:

```json
{
  "RunFastenloc.GWASManifest": "gwas_manifest.tsv",
  "RunFastenloc.QTLData": [
    "eqtl.fastenloc.vcf.gz",
    "sqtl.fastenloc.vcf.gz",
    "pqtl.fastenloc.vcf.gz"
  ],
  "RunFastenloc.QTLLabels": ["eQTL", "sQTL", "pQTL"],
  "RunFastenloc.consensus_jaccard": 0.90
}
```

The manifest is tab-delimited:

```text
study_id	trait	n_variants	gwas_path	trait_category	n_credible_sets
GCST90027158	Alzheimer disease	1000000	gs://bucket/alzheimers.fastenloc.vcf.gz	neuro	245
FINNGEN_R12_G6_MS	Multiple sclerosis	950000	gs://bucket/ms.fastenloc.vcf.gz	immune	180
```

For a single QTL layer, provide one-element `QTLData` and `QTLLabels` arrays.

`gwas_path` is localized inside each GWAS job. It can point to `gs://`, HTTP(S), or a path already accessible inside the task runtime. The Docker image includes `gsutil` for Google Cloud Storage paths. Each manifest row is expected to point to one trait/study analysis unit.

## Primary Outputs

The workflow emits all-GWAS/all-QTL combined outputs plus per-GWAS and per-GWAS-by-QTL output arrays. The most useful downstream outputs are:

| Output | Description |
| --- | --- |
| `consensus_loci_out` | Trait-level consensus map from each original GWAS credible set to `consensus_locus_id`. |
| `harmonized_signal_out` | Signal-level table joining fastENLOC RCP/LCP, CLPP, and gene-level GRCP/GLCP. |
| `harmonized_credible_set_out` | Credible-set-level rollup with colocalization flags, per-method gene lists, and `consensus_locus_id` for de-duplicated trait coverage. |
| `harmonized_gene_out` | Gene-level rollup with best signal metrics and gene-native fastENLOC metrics. |

Raw combined fastENLOC and CLPP outputs include leading GWAS metadata columns plus `qtl_label`; the manifest trait is named `gwas_trait` there to avoid colliding with fastENLOC's own `trait` column. Harmonized outputs store the QTL label in `layer` and include `study`, `trait`, `trait_category`, `n_variants`, and `n_credible_sets`.

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
│   ├── harmonize_coloc.R
│   └── merge_credible_sets.R
└── workflows/
    ├── tasks/
    │   ├── aggregation.wdl
    │   ├── clpp.wdl
    │   ├── consensus.wdl
    │   ├── fastenloc.wdl
    │   ├── harmonize.wdl
    │   ├── input_validation.wdl
    │   ├── localize.wdl
    │   └── split.wdl
    └── RunFastEnlocTraits.wdl
```

## Runtime Image

The WDL runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

After adding or changing files under `scripts/`, rebuild and publish the image so the workflow runtime contains the latest helper scripts.
