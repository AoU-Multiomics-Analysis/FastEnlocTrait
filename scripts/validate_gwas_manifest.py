#!/usr/bin/env python3
"""Validate every fastENLOC GWAS file referenced by a pipeline manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from validate_fastenloc_gwas import validate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"{args.manifest}: manifest contains no studies")

    results: list[dict[str, object]] = []
    failures = 0
    for row in rows:
        study_id = row["study_id"]
        try:
            result = validate(
                Path(row["gwas_path"]),
                study_id,
                int(row["n_credible_sets"]),
            )
            result["error"] = ""
        except (OSError, ValueError) as exc:
            failures += 1
            result = {
                "study_id": study_id,
                "status": "FAIL",
                "n_fastenloc_rows": "",
                "expected_credible_sets": row["n_credible_sets"],
                "observed_credible_sets": "",
                "min_set_size": "",
                "max_set_size": "",
                "min_credible_set_pip_sum": "",
                "max_credible_set_pip_sum": "",
                "error": str(exc).replace("\n", " | "),
            }
        results.append(result)

    headers = [
        "study_id",
        "status",
        "n_fastenloc_rows",
        "expected_credible_sets",
        "observed_credible_sets",
        "min_set_size",
        "max_set_size",
        "min_credible_set_pip_sum",
        "max_credible_set_pip_sum",
        "error",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=headers,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(results)
    print(
        f"Validated {len(results)} studies: "
        f"{len(results) - failures} passed, {failures} failed"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
