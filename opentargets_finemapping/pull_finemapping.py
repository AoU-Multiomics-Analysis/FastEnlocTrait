#!/usr/bin/env python3
"""Re-pull configured Open Targets fine-mapping data and build fastENLOC inputs.

The input manifest intentionally contains no file paths. Each row identifies a
study, the fine-mapping method to retain, and metadata to carry into the output
pipeline manifest. GraphQL responses are cached so interrupted downloads can
resume without repeating completed requests.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


DEFAULT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "manifests" / "opentargets_studies.no_paths.tsv"
QUERY_DIR = ROOT / "queries"
REQUIRED_COLUMNS = {
    "study_id",
    "trait",
    "trait_category",
    "finemap_method",
    "n_credible_sets",
    "n_variants",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download all configured Open Targets credible sets, retrieve complete "
            "nested loci, build deterministic fastENLOC GWAS inputs, and emit a "
            "pipeline-ready manifest."
        )
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=Path("opentargets_pull"))
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    parser.add_argument("--study-id", action="append", default=[])
    parser.add_argument("--max-studies", type=int)
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--credible-set-page-size", type=int, default=100)
    parser.add_argument("--locus-page-size", type=int, default=500)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--retry-wait", type=float, default=2.0)
    parser.add_argument("--pip-sum-min", type=float, default=0.90)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--allow-count-mismatch", action="store_true")
    parser.add_argument("--keep-going", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and summarize the path-redacted manifest without network calls.",
    )
    args = parser.parse_args(argv)
    for name in (
        "batch_size",
        "credible_set_page_size",
        "locus_page_size",
        "timeout",
        "retries",
    ):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.max_studies is not None and args.max_studies <= 0:
        parser.error("--max-studies must be positive")
    return args


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"{path}: manifest has no header")
        headers = list(reader.fieldnames)
        rows = list(reader)
    missing = sorted(REQUIRED_COLUMNS - set(headers))
    if missing:
        raise RuntimeError(f"{path}: missing required columns: {', '.join(missing)}")
    if "gwas_path" in headers and any(row.get("gwas_path", "").strip() for row in rows):
        raise RuntimeError(
            f"{path}: retrieval manifest must not contain populated gwas_path values"
        )
    if not rows:
        raise RuntimeError(f"{path}: manifest has no study rows")
    ids = [row["study_id"].strip() for row in rows]
    if any(not value for value in ids):
        raise RuntimeError(f"{path}: study_id cannot be empty")
    duplicates = sorted({value for value in ids if ids.count(value) > 1})
    if duplicates:
        raise RuntimeError(f"{path}: duplicate study IDs: {', '.join(duplicates)}")
    for row in rows:
        for field in ("n_credible_sets", "n_variants"):
            try:
                value = int(row[field])
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: {row['study_id']} has invalid {field}={row[field]!r}"
                ) from exc
            if value < 0 or (field == "n_variants" and value == 0):
                raise RuntimeError(
                    f"{path}: {row['study_id']} has invalid {field}={value}"
                )
    return headers, rows


def write_tsv(path: Path, headers: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=headers,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def chunks(values: list[Any], size: int) -> Iterable[list[Any]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def normalize_method(value: str | None) -> str:
    if value == "SuSie":
        return "SuSiE"
    return (value or "unknown").strip()


def configured_methods(value: str) -> set[str]:
    return {
        normalize_method(item.strip())
        for item in re.split(r"[;,|]", value)
        if item.strip()
    }


def slugify(value: str, max_length: int = 120) -> str:
    clean = value.replace("'", "")
    clean = re.sub(r"[^A-Za-z0-9]+", "_", clean).strip("_")
    return clean[:max_length].rstrip("_") or "trait"


def chromosome_key(value: Any) -> tuple[int, str]:
    clean = str(value).removeprefix("chr")
    aliases = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    try:
        return int(clean), clean
    except ValueError:
        return aliases.get(clean.upper(), 99), clean


class GraphQLClient:
    def __init__(
        self,
        endpoint: str,
        cache_dir: Path,
        *,
        timeout: int,
        retries: int,
        retry_wait: float,
        refresh: bool,
    ) -> None:
        self.endpoint = endpoint
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.retries = retries
        self.retry_wait = retry_wait
        self.refresh = refresh

    def execute(
        self,
        query: str,
        variables: dict[str, Any],
        cache_name: str,
    ) -> dict[str, Any]:
        cache_path = self.cache_dir / cache_name
        request_payload = {"query": query, "variables": variables}
        if cache_path.exists() and not self.refresh:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("request") != request_payload:
                raise RuntimeError(
                    f"Cached request mismatch at {cache_path}; use --refresh or a new output directory"
                )
            return self._response_data(cached["response"], cache_path)

        encoded = json.dumps(request_payload).encode("utf-8")
        request = urllib.request.Request(
            self.endpoint,
            data=encoded,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "FastEnlocTrait-OpenTargets-retriever/1.0",
            },
            method="POST",
        )
        response_payload: dict[str, Any] | None = None
        for attempt in range(1, self.retries + 1):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    response_payload = json.loads(response.read().decode("utf-8"))
                break
            except urllib.error.HTTPError as exc:
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if not retryable or attempt == self.retries:
                    body = exc.read().decode("utf-8", errors="replace")[:1000]
                    raise RuntimeError(
                        f"Open Targets HTTP {exc.code} for {cache_name}: {body}"
                    ) from exc
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                if attempt == self.retries:
                    raise RuntimeError(
                        f"Open Targets request failed for {cache_name}: {exc}"
                    ) from exc
            time.sleep(self.retry_wait * (2 ** (attempt - 1)))

        if response_payload is None:
            raise RuntimeError(f"No response received for {cache_name}")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(
            json.dumps(
                {"request": request_payload, "response": response_payload},
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return self._response_data(response_payload, cache_path)

    @staticmethod
    def _response_data(response: dict[str, Any], cache_path: Path) -> dict[str, Any]:
        if response.get("errors"):
            raise RuntimeError(
                f"Open Targets GraphQL error in {cache_path}: "
                f"{json.dumps(response['errors'])[:2000]}"
            )
        data = response.get("data")
        if not isinstance(data, dict):
            raise RuntimeError(f"Missing GraphQL data object in {cache_path}")
        return data


def query_fingerprint(study_ids: list[str]) -> str:
    return hashlib.sha256("\n".join(study_ids).encode("utf-8")).hexdigest()[:12]


def fetch_credible_sets(
    client: GraphQLClient,
    query: str,
    study_ids: list[str],
    *,
    batch_size: int,
    outer_page_size: int,
    locus_page_size: int,
) -> list[dict[str, Any]]:
    all_rows: list[dict[str, Any]] = []
    for batch in chunks(study_ids, batch_size):
        fingerprint = query_fingerprint(batch)
        page_index = 0
        batch_rows: list[dict[str, Any]] = []
        expected: int | None = None
        while True:
            variables = {
                "studyIds": batch,
                "page": {"index": page_index, "size": outer_page_size},
                "locusPage": {"index": 0, "size": locus_page_size},
            }
            data = client.execute(
                query,
                variables,
                f"credible_sets/{fingerprint}.page_{page_index:05d}.json",
            )
            block = data.get("credibleSets")
            if not isinstance(block, dict):
                raise RuntimeError("credibleSets query returned no credibleSets block")
            count = int(block["count"])
            expected = count if expected is None else expected
            if count != expected:
                raise RuntimeError(f"Inconsistent credible-set count for batch {batch}")
            page_rows = block.get("rows") or []
            batch_rows.extend(page_rows)
            if len(batch_rows) >= expected:
                break
            if not page_rows:
                raise RuntimeError(
                    f"Credible-set pagination stopped at {len(batch_rows)}/{expected} for {batch}"
                )
            page_index += 1
        if len(batch_rows) != expected:
            raise RuntimeError(
                f"Retrieved {len(batch_rows)} credible sets, expected {expected}, for {batch}"
            )
        all_rows.extend(batch_rows)

    locus_ids = [row["studyLocusId"] for row in all_rows]
    if len(locus_ids) != len(set(locus_ids)):
        raise RuntimeError("Duplicate studyLocusId values across credible-set pages")
    returned_studies = {row["studyId"] for row in all_rows}
    unexpected = returned_studies - set(study_ids)
    if unexpected:
        raise RuntimeError(f"API returned unexpected studies: {sorted(unexpected)}")
    return all_rows


def complete_locus_pages(
    client: GraphQLClient,
    query: str,
    credible_sets: list[dict[str, Any]],
    *,
    locus_page_size: int,
) -> None:
    for credible_set in credible_sets:
        locus = credible_set.get("locus") or {}
        expected = int(locus.get("count") or 0)
        rows = list(locus.get("rows") or [])
        page_index = 1
        while len(rows) < expected:
            variables = {
                "studyLocusId": credible_set["studyLocusId"],
                "locusPage": {"index": page_index, "size": locus_page_size},
            }
            data = client.execute(
                query,
                variables,
                (
                    "locus_pages/"
                    f"{credible_set['studyLocusId']}.page_{page_index:05d}.json"
                ),
            )
            returned = data.get("credibleSet")
            if not isinstance(returned, dict):
                raise RuntimeError(
                    f"{credible_set['studyLocusId']}: missing credibleSet response"
                )
            block = returned.get("locus") or {}
            if int(block.get("count") or 0) != expected:
                raise RuntimeError(
                    f"{credible_set['studyLocusId']}: locus count changed during pagination"
                )
            page_rows = block.get("rows") or []
            if not page_rows:
                raise RuntimeError(
                    f"{credible_set['studyLocusId']}: incomplete locus pagination "
                    f"({len(rows)}/{expected})"
                )
            rows.extend(page_rows)
            page_index += 1
        if len(rows) != expected:
            raise RuntimeError(
                f"{credible_set['studyLocusId']}: retrieved {len(rows)}/{expected} locus rows"
            )
        variant_ids = [row["variant"]["id"] for row in rows]
        if len(variant_ids) != len(set(variant_ids)):
            raise RuntimeError(
                f"{credible_set['studyLocusId']}: duplicate variants across locus pages"
            )
        credible_set["locus"]["rows"] = rows


def deterministic_gzip(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                text.writelines(lines)
    temp.replace(path)


def build_study_file(
    config: dict[str, str],
    credible_sets: list[dict[str, Any]],
    output_dir: Path,
    *,
    pip_sum_min: float,
    allow_count_mismatch: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    study_id = config["study_id"]
    methods = configured_methods(config["finemap_method"])
    selected = [
        row
        for row in credible_sets
        if row["studyId"] == study_id
        and normalize_method(row.get("finemappingMethod")) in methods
    ]
    selected.sort(
        key=lambda row: (
            chromosome_key(row.get("chromosome")),
            row.get("position") or 0,
            row["studyLocusId"],
        )
    )
    expected_count = int(config["n_credible_sets"])
    if len(selected) != expected_count and not allow_count_mismatch:
        available = sorted(
            {
                normalize_method(row.get("finemappingMethod"))
                for row in credible_sets
                if row["studyId"] == study_id
            }
        )
        raise RuntimeError(
            f"{study_id}: selected {len(selected)} credible sets for methods "
            f"{sorted(methods)}, expected {expected_count}; available methods={available}"
        )
    if not selected:
        raise RuntimeError(f"{study_id}: no credible sets matched configured methods {methods}")

    trait_label = slugify(config["trait"])
    filename = f"{trait_label}__{study_id}.fastenloc.gwas.vcf.gz"
    output_path = (output_dir / "gwas" / filename).resolve()
    lines: list[str] = []
    unique_variants: set[str] = set()
    pip_sums: list[float] = []
    seen_sets: set[str] = set()

    for credible_set in selected:
        study_locus_id = credible_set["studyLocusId"]
        if study_locus_id in seen_sets:
            raise RuntimeError(f"{study_id}: duplicate studyLocusId {study_locus_id}")
        seen_sets.add(study_locus_id)
        locus_rows = [
            row
            for row in credible_set["locus"]["rows"]
            if row.get("is95CredibleSet") is True
        ]
        if not locus_rows:
            raise RuntimeError(f"{study_locus_id}: empty 95% credible set")
        locus_rows.sort(
            key=lambda row: (
                chromosome_key(row["variant"]["chromosome"]),
                int(row["variant"]["position"]),
                row["variant"]["id"],
            )
        )
        pip_sum = sum(float(row["posteriorProbability"]) for row in locus_rows)
        if pip_sum < pip_sum_min or pip_sum > 1.000001:
            raise RuntimeError(
                f"{study_locus_id}: unexpected 95% credible-set PIP sum {pip_sum:.8f}"
            )
        pip_sums.append(pip_sum)
        locus_chromosomes = {
            str(row["variant"]["chromosome"]).removeprefix("chr")
            for row in locus_rows
        }
        if len(locus_chromosomes) != 1:
            raise RuntimeError(
                f"{study_locus_id}: credible-set variants span chromosomes "
                f"{sorted(locus_chromosomes)}"
            )
        locus_chromosome = next(iter(locus_chromosomes))
        locus_start = min(int(row["variant"]["position"]) for row in locus_rows)
        locus_end = max(int(row["variant"]["position"]) for row in locus_rows)
        locus_token = slugify(study_locus_id)
        locus_suffix = (
            locus_token if locus_token.lower().startswith("l")
            else f"L{locus_token}"
        )
        set_id = (
            f"{study_id}_chr{locus_chromosome}."
            f"{locus_start}.{locus_end}_{locus_suffix}"
        )
        set_size = len(locus_rows)
        for row in locus_rows:
            variant = row["variant"]
            chrom = str(variant["chromosome"]).removeprefix("chr")
            chrom_label = f"chr{chrom}"
            position = int(variant["position"])
            reference = variant["referenceAllele"]
            alternate = variant["alternateAllele"]
            variant_id = f"{chrom_label}_{position}_{reference}_{alternate}"
            unique_variants.add(variant["id"])
            annotation = (
                f"{trait_label};{set_id}={float(row['posteriorProbability']):.4e}"
                f"[{pip_sum:.3e}:{set_size}]"
            )
            lines.append(
                "\t".join(
                    [
                        chrom_label,
                        str(position),
                        variant_id,
                        reference,
                        alternate,
                        annotation,
                    ]
                )
                + "\n"
            )
    deterministic_gzip(output_path, lines)

    manifest_row = dict(config)
    manifest_row["n_credible_sets"] = str(len(selected))
    manifest_row["gwas_path"] = str(output_path)
    qc_row = {
        "study_id": study_id,
        "status": "complete",
        "configured_methods": ";".join(sorted(methods)),
        "n_credible_sets_expected": expected_count,
        "n_credible_sets_written": len(selected),
        "n_fastenloc_rows": len(lines),
        "n_unique_variants": len(unique_variants),
        "min_credible_set_pip_sum": f"{min(pip_sums):.8f}",
        "max_credible_set_pip_sum": f"{max(pip_sums):.8f}",
        "gwas_path": str(output_path),
        "error": "",
    }
    return manifest_row, qc_row


def output_headers(input_headers: list[str]) -> list[str]:
    headers = [header for header in input_headers if header != "gwas_path"]
    insert_after = "n_samples" if "n_samples" in headers else "n_variants"
    headers.insert(headers.index(insert_after) + 1, "gwas_path")
    return headers


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    headers, configs = read_tsv(args.manifest)
    if args.study_id:
        requested = set(args.study_id)
        known = {row["study_id"] for row in configs}
        missing = sorted(requested - known)
        if missing:
            raise RuntimeError(f"Requested study IDs are absent from manifest: {missing}")
        configs = [row for row in configs if row["study_id"] in requested]
    if args.max_studies is not None:
        configs = configs[: args.max_studies]

    method_counts: dict[str, int] = defaultdict(int)
    for row in configs:
        method_counts[row["finemap_method"]] += 1
    print(
        json.dumps(
            {
                "manifest": str(args.manifest.resolve()),
                "studies": len(configs),
                "fine_mapping_methods": dict(sorted(method_counts.items())),
                "output_dir": str(args.output_dir.resolve()),
                "dry_run": args.dry_run,
            },
            indent=2,
        )
    )
    if args.dry_run:
        return 0

    credible_query = (QUERY_DIR / "credible_sets.graphql").read_text(encoding="utf-8")
    locus_query = (QUERY_DIR / "credible_set_locus_page.graphql").read_text(
        encoding="utf-8"
    )
    output_dir = args.output_dir.resolve()
    client = GraphQLClient(
        args.endpoint,
        output_dir / "provenance",
        timeout=args.timeout,
        retries=args.retries,
        retry_wait=args.retry_wait,
        refresh=args.refresh,
    )
    study_ids = [row["study_id"] for row in configs]
    credible_sets = fetch_credible_sets(
        client,
        credible_query,
        study_ids,
        batch_size=args.batch_size,
        outer_page_size=args.credible_set_page_size,
        locus_page_size=args.locus_page_size,
    )
    complete_locus_pages(
        client,
        locus_query,
        credible_sets,
        locus_page_size=args.locus_page_size,
    )

    by_study: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for credible_set in credible_sets:
        by_study[credible_set["studyId"]].append(credible_set)

    manifest_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for index, config in enumerate(configs, start=1):
        study_id = config["study_id"]
        print(f"[{index}/{len(configs)}] building {study_id}", file=sys.stderr)
        try:
            manifest_row, qc_row = build_study_file(
                config,
                by_study.get(study_id, []),
                output_dir,
                pip_sum_min=args.pip_sum_min,
                allow_count_mismatch=args.allow_count_mismatch,
            )
        except Exception as exc:  # noqa: BLE001
            if not args.keep_going:
                raise
            qc_rows.append(
                {
                    "study_id": study_id,
                    "status": "error",
                    "configured_methods": config["finemap_method"],
                    "n_credible_sets_expected": config["n_credible_sets"],
                    "n_credible_sets_written": 0,
                    "n_fastenloc_rows": 0,
                    "n_unique_variants": 0,
                    "min_credible_set_pip_sum": "",
                    "max_credible_set_pip_sum": "",
                    "gwas_path": "",
                    "error": str(exc),
                }
            )
            print(f"ERROR {study_id}: {exc}", file=sys.stderr)
            continue
        manifest_rows.append(manifest_row)
        qc_rows.append(qc_row)

    manifest_path = output_dir / "gwas_manifest.tsv"
    qc_path = output_dir / "retrieval_qc.tsv"
    write_tsv(manifest_path, output_headers(headers), manifest_rows)
    qc_headers = list(qc_rows[0]) if qc_rows else [
        "study_id",
        "status",
        "error",
    ]
    write_tsv(qc_path, qc_headers, qc_rows)
    failures = sum(row["status"] != "complete" for row in qc_rows)
    print(
        json.dumps(
            {
                "credible_sets_downloaded": len(credible_sets),
                "studies_complete": len(manifest_rows),
                "studies_failed": failures,
                "pipeline_manifest": str(manifest_path),
                "qc": str(qc_path),
                "provenance": str(output_dir / "provenance"),
            },
            indent=2,
        )
    )
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
