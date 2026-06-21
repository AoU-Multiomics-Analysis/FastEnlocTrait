# FastEnlocTrait

Workflow and helper scripts for running [fastENLOC](https://github.com/xqwen/fastenloc) colocalization across multiple GWAS traits and QTL fine-mapping datasets. The main workflow splits a multi-trait GWAS/trait input into smaller chunks, runs `fastenloc` for each trait against a QTL input, and aggregates the standard `.enloc.*.out` result files.

## Repository layout

```text
.
├── envs/
│   └── Dockerfile
├── scripts/
│   ├── PrepMVPFineMapping.R
│   ├── PrepQTLFinemapping.R
│   └── SplitTraitData.R
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
| `QTLData` | `File` | QTL fastENLOC input file, formatted for `fastenloc -eqtl`. |
| `NumberVariants` | `Int` | Total number of variants passed to `fastenloc -total_variants`. |

`SplitFastenloc` also supports `traits_per_chunk`, which defaults to `25` inside the task.

### Outputs

The workflow aggregates per-trait fastENLOC outputs into:

| Output | Description |
| --- | --- |
| `combined_gene_out` | Combined `*.enloc.gene.out` results. |
| `combined_enrich_out` | Combined `*.enloc.enrich.out` results. |
| `combined_mi_out` | Combined `*.enloc.mi.out` results. |
| `combined_sig_out` | Combined `*.enloc.sig.out` results. |
| `combined_snp_out` | Combined `*.enloc.snp.out` results. |

Each per-trait output is annotated with an added leading `trait` column before aggregation.

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
  "RunFastenloc.QTLData": "qtl.fastenloc.vcf.gz",
  "RunFastenloc.NumberVariants": 1000000
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
