# Re-pull Open Targets fine-mapping data

This directory provides a reproducible, standard-library-only downloader and
converter for the Open Targets fine-mapping studies used by FastEnlocTrait. It
re-pulls every study listed in the path-redacted manifest, retrieves complete
credible-set and nested-locus pagination, validates the downloaded data, writes
deterministic fastENLOC GWAS inputs, and creates a pipeline-ready manifest.

The implementation uses the documented
[Open Targets Platform GraphQL API](https://platform-docs.opentargets.org/data-access/graphql-api).
Open Targets recommends its bulk data downloads for unconstrained,
platform-wide extraction; this utility is intentionally limited to a curated
study snapshot and batches and caches its API requests.

The checked-in manifest is an analysis snapshot, not an automatic discovery
query. This distinction is intentional: rerunning it retrieves the same set of
approved studies instead of silently adding new phenotypes released later by
Open Targets.

## Included files

- `pull_finemapping.py`: downloader, cache/resume logic, QC, conversion, and
  pipeline-manifest writer.
- `manifests/opentargets_studies.no_paths.tsv`: the approved 401-study
  configuration. It contains no local or cloud file paths. MTAG and explicit
  multi-trait/pleiotropic studies are excluded.
- `queries/credible_sets.graphql`: outer credible-set query with the first page
  of each nested locus.
- `queries/credible_set_locus_page.graphql`: additional nested-locus pages.
- `tests/test_pull_finemapping.py`: offline conversion and manifest tests.

## Requirements

- Python 3.10 or newer.
- HTTPS access to `https://api.platform.opentargets.org/api/v4/graphql`.
- No third-party Python packages are required.

## Validate the configuration without downloading

From the FastEnlocTrait repository root:

```bash
python3 opentargets_finemapping/pull_finemapping.py --dry-run
```

## Re-pull every configured study

```bash
python3 opentargets_finemapping/pull_finemapping.py \
  --output-dir opentargets_pull
```

For a current-release refresh where credible-set counts may legitimately differ
from the checked-in snapshot, add `--allow-count-mismatch`. The generated QC
table retains both the snapshot count and the newly written count:

```bash
python3 opentargets_finemapping/pull_finemapping.py \
  --allow-count-mismatch \
  --output-dir opentargets_pull
```

The run writes:

- `opentargets_pull/gwas/*.fastenloc.gwas.vcf.gz`
- `opentargets_pull/gwas_manifest.tsv`
- `opentargets_pull/retrieval_qc.tsv`
- `opentargets_pull/provenance/credible_sets/*.json`
- `opentargets_pull/provenance/locus_pages/*.json`

The generated `gwas_manifest.tsv` contains absolute local `gwas_path` values
and can be supplied directly as `RunFastenloc.GWASManifest` for a local WDL
run.

Every emitted credible-set identifier contains its chromosome, observed
variant interval, and Open Targets `studyLocusId`, for example
`GCST010571_chr1.2561226.2581666_L<studyLocusId>`. This prevents multiple
fine-mapped sets on the same chromosome from collapsing into one consensus
bucket.

## Test a subset first

```bash
python3 opentargets_finemapping/pull_finemapping.py \
  --study-id GCST011956 \
  --output-dir opentargets_pull_sle
```

Alternatively, use `--max-studies 3` to run the first three manifest rows.

## Resume, refresh, and failures

Completed GraphQL requests are cached under `provenance/`. Re-running the same
command resumes from those files. Use `--refresh` to replace the cache.

By default, a count mismatch, incomplete pagination, duplicate locus variant,
empty 95% credible set, or unexpected PIP sum stops the run. This protects the
pipeline from silently accepting API or schema changes. Useful controls:

- `--keep-going`: record failed studies in `retrieval_qc.tsv` and continue.
- `--allow-count-mismatch`: accept a current Open Targets credible-set count
  that differs from the checked-in snapshot.
- `--pip-sum-min`: change the default minimum 95% credible-set PIP sum of
  `0.90`.
- `--credible-set-page-size`, `--locus-page-size`, and `--batch-size`: tune
  GraphQL pagination.

A live schema/conversion check on 2026-07-17 found one example of normal release
drift: `GCST004030` had 10 current SuSiE-inf credible sets versus 11 in the
snapshot. Strict mode detected the difference; `--allow-count-mismatch`
successfully wrote the current 10-set input and recorded both values in QC.

If Open Targets changes its GraphQL schema, update the query documents and the
small response-parsing sections in `pull_finemapping.py`; cached requests
include their exact query and variables so stale responses cannot be reused
accidentally.
