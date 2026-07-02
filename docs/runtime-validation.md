# Runtime and Validation

## Docker Image

The workflow runtime uses:

```text
ghcr.io/aou-multiomics-analysis/fastenloctrait:main
```

The Dockerfile installs R dependencies, builds `fastenloc` from `xqwen/fastenloc`, copies scripts into `/home/mambauser`, and makes the `fastenloc` binary available on `PATH`.

Build locally with:

```bash
docker build -t fastenloctrait:local -f envs/Dockerfile .
```

After changing files under `scripts/`, rebuild and publish the runtime image before running the WDL in Dockstore/Cromwell so the container includes the latest helper scripts.

## Local Validation

Check WDL syntax and static issues:

```bash
miniwdl check workflows/RunFastEnlocTraits.wdl
```

Syntax-check all R scripts:

```bash
Rscript -e 'for (f in list.files("scripts", pattern = "[.]R$", full.names = TRUE)) { parse(f); cat("parse OK:", f, "\n") }'
```

Check for whitespace issues before committing:

```bash
git diff --check
```

## End-To-End Testing

This repository does not currently include small example GWAS/QTL fixtures. Full end-to-end validation requires representative fastENLOC-format trait and QTL files plus a WDL runner such as Cromwell, miniwdl, or Dockstore.
