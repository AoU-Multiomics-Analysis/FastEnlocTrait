# Runtime and Validation

## Docker Image

The workflow runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

The Dockerfile installs R dependencies, builds `fastenloc` from `xqwen/fastenloc`, copies scripts into `/home/mambauser`, and makes the `fastenloc` binary available on `PATH`.

It also installs `google-cloud-sdk`, `curl`, and `gzip` so `LocalizeGWASData` can materialize manifest `gwas_path` values from `gs://`, HTTP(S), or in-runtime local paths inside each GWAS job.

Build locally with:

```bash
docker build -t fastenloctrait:local -f envs/Dockerfile .
```

After changing files under `scripts/`, `envs/`, or runtime task requirements, rebuild and publish the runtime image before running the WDL in Dockstore/Cromwell so the container includes the latest helper scripts and localization tools.

## Local Validation

Check WDL syntax and static issues:

```bash
miniwdl check workflows/RunFastEnlocTraits.wdl
```

Syntax-check all R scripts:

```bash
Rscript -e 'for (f in list.files("scripts", pattern = "[.]R$", full.names = TRUE)) { parse(f); cat("parse OK:", f, "\n") }'
```

Run the final summary-plot unit test:

```bash
Rscript tests/test_plot_coloc_summary.R
Rscript tests/test_summarize_raw_coloc.R
```

Run a representative local data smoke test without a WDL engine:

```bash
scripts/run_local_smoke_test.sh \
  "/path/to/Open targets fine-mapping" \
  "/path/to/susie_files" \
  test_runs/opentargets_smoke \
  /path/to/fastenloc \
  /path/to/gencode.annotation.gtf.gz
```

The runner selects the study with the most credible sets in each trait category, runs all three QTL layers, and preserves intermediate/raw data, logs, harmonized outputs, summaries, figures, and a QC table. It reuses completed pair outputs when rerun against the same output directory.

To run one representative study for every distinct trait, set
`SELECTION_MODE=trait_max`. Use `SELECTION_MODE=all` to retain every manifest
study. For each study, the eQTL, sQTL, and pQTL analyses run concurrently.

```bash
SELECTION_MODE=trait_max scripts/run_local_smoke_test.sh \
  "/path/to/Open targets fine-mapping" \
  "/path/to/susie_files" \
  test_runs/opentargets_all_traits \
  /path/to/fastenloc
```

The runner writes plot variants requiring at least 1, 5, and 10 distinct
colocalizing genes. These thresholds affect only the figures and their audit
TSVs, not the underlying analysis or complete summaries.

Run the fastENLOC parsing and aggregation regression tests:

```bash
bash tests/test_fastenloc_output_normalization.sh
Rscript tests/test_fastenloc_parsing.R
Rscript tests/test_gzip_aggregation.R
```

Check for whitespace issues before committing:

```bash
git diff --check
```

## End-To-End Testing

This repository does not currently include small example GWAS/QTL fixtures. Full end-to-end validation requires representative fastENLOC-format trait and QTL files plus a WDL runner such as Cromwell, miniwdl, or Dockstore.
