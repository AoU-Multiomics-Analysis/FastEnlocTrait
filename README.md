# FastEnlocTrait

FastEnlocTrait is a WDL workflow for running colocalization across one or more GWAS analysis units and one or more QTL layers, such as eQTL, sQTL, and pQTL. It runs fastENLOC, computes CLPP from the same fastENLOC-format inputs, merges GWAS credible sets into trait-level consensus loci, harmonizes the coloc metrics, and summarizes trait-level coverage and gene results.

## What the Workflow Does

1. Validates a GWAS manifest and QTL layer labels.
2. Localizes each manifest `gwas_path` inside the job.
3. Merges GWAS credible sets across studies of the same trait into consensus loci.
4. Runs fastENLOC and CLPP for each GWAS/QTL pair in configurable GWAS shards.
5. Aggregates per-GWAS, per-QTL, and all-GWAS outputs.
6. Summarizes raw fastENLOC and CLPP output families by QTL layer.
7. Harmonizes fastENLOC and CLPP results at signal, credible-set, and gene levels, with consensus locus IDs on credible-set rollups.
8. Summarizes de-duplicated colocalization rates and colocalizing genes by trait, QTL layer, and cross-layer union.
9. Packages the main analysis-ready outputs into one archive for easier browsing.

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
  "RunFastenloc.GTF": "gencode.v44.annotation.gtf.gz",
  "RunFastenloc.QTLData": [
    "eqtl.fastenloc.vcf.gz",
    "sqtl.fastenloc.vcf.gz",
    "pqtl.fastenloc.vcf.gz"
  ],
  "RunFastenloc.QTLLabels": ["eQTL", "sQTL", "pQTL"],
  "RunFastenloc.gwas_units_per_shard": 10,
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

`gwas_units_per_shard` controls how many manifest rows are processed by each raw fastENLOC/CLPP shard job. Harmonization still runs separately for each GWAS x QTL pair, so Bayesian FDR thresholds remain per GWAS analysis unit and QTL layer.

## Primary Outputs

The workflow emits all-GWAS/all-QTL combined outputs plus per-GWAS and per-GWAS-by-QTL output arrays. The most useful downstream outputs are:

| Output | Description |
| --- | --- |
| `consensus_loci_out` | Trait-level consensus map from each original GWAS credible set to `consensus_locus_id`. |
| `raw_coloc_summary_out` | QTL-layer summary of each raw fastENLOC output family plus CLPP. |
| `per_qtl_raw_coloc_summary_out` | One raw fastENLOC/CLPP summary TSV per QTL label. |
| `harmonized_signal_out` | Signal-level table joining fastENLOC RCP/LCP, CLPP, and gene-level GRCP/GLCP. |
| `harmonized_credible_set_out` | Credible-set-level rollup with colocalization flags, per-method gene lists, and `consensus_locus_id` for de-duplicated trait coverage. |
| `harmonized_gene_out` | Gene-level rollup retaining all fastENLOC gene-level rows, with best signal/CLPP metrics joined when available. |
| `coloc_rate_by_trait_out` | Trait x layer and cross-layer union colocalization rates across consensus loci. |
| `gene_summary_by_trait_out` | Consensus-locus-based trait x layer and cross-layer union gene counts, genes per locus, protein-coding counts, and gene IDs. |
| `colocalizing_genes_long_out` | Long table of colocalizing genes by trait, layer, and consensus locus. |
| `high_level_outputs_archive` | Tarball containing the main manifests, consensus files, raw combined outputs, harmonized outputs, and final summary tables. |

Raw combined fastENLOC and CLPP outputs include leading GWAS metadata columns plus `qtl_label`; the manifest trait is named `gwas_trait` there to avoid colliding with fastENLOC's own `trait` column. Harmonized outputs store the QTL label in `layer` and include `study`, `trait`, `trait_category`, `n_variants`, and `n_credible_sets`.

Use `harmonized_gene_out` for complete fastENLOC gene-level evidence, including high-GLCP/GRCP genes without matched signal-level coloc rows. Use `gene_summary_by_trait_out` for downstream coverage-style gene summaries tied to consensus loci.

## Documentation

- [Workflow reference](docs/workflow.md): WDL inputs, outputs, multi-QTL behavior, and output naming.
- [Helper scripts](docs/helper-scripts.md): R script usage for preparing inputs, CLPP, harmonization, and summaries.
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
│   ├── merge_credible_sets.R
│   ├── summarize_coloc.R
│   └── summarize_raw_coloc.R
└── workflows/
    ├── tasks/
    │   ├── aggregation.wdl
    │   ├── clpp.wdl
    │   ├── consensus.wdl
    │   ├── fastenloc.wdl
    │   ├── harmonize.wdl
    │   ├── input_validation.wdl
    │   ├── localize.wdl
    │   ├── shard_coloc.wdl
    │   ├── split.wdl
    │   └── summarize.wdl
    └── RunFastEnlocTraits.wdl
```

## Runtime Image

The WDL runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

After adding or changing files under `scripts/`, rebuild and publish the image so the workflow runtime contains the latest helper scripts.
