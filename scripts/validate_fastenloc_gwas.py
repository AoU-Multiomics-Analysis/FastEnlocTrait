#!/usr/bin/env python3
"""Validate a localized fastENLOC GWAS file before colocalization."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path


ANNOTATION = re.compile(
    r"^(?P<set_id>[^=]+)=(?P<pip>[0-9.eE+-]+)"
    r"\[(?P<cs_pip>[0-9.eE+-]+):(?P<set_size>[0-9]+)\]$"
)
PLACEHOLDER_SET_ID = re.compile(
    r"(?:_L(?:none|na|nan)$|chr(?:[0-9]+|X|Y|M|MT)\.0\.0_)",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas", required=True, type=Path)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--expected-credible-sets", required=True, type=int)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def validate(
    path: Path, study_id: str, expected_credible_sets: int
) -> dict[str, object]:
    if expected_credible_sets < 1:
        raise ValueError("expected credible-set count must be positive")

    row_counts: dict[str, int] = defaultdict(int)
    declared_sizes: dict[str, int] = {}
    declared_cs_pip: dict[str, float] = {}
    pip_sums: dict[str, float] = defaultdict(float)
    chromosomes: dict[str, set[str]] = defaultdict(set)
    variants: dict[str, set[str]] = defaultdict(set)
    errors: list[str] = []
    n_rows = 0

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, start=1):
            if not fields or all(not value for value in fields):
                continue
            n_rows += 1
            if len(fields) != 6:
                errors.append(
                    f"line {line_number}: expected 6 tab-separated fields, got {len(fields)}"
                )
                continue
            chrom, position, variant_id, reference, alternate, annotation = fields
            if not chrom or not position or not variant_id or not reference or not alternate:
                errors.append(f"line {line_number}: required variant field is empty")
                continue
            try:
                if int(position) < 1:
                    raise ValueError
            except ValueError:
                errors.append(f"line {line_number}: invalid position {position!r}")
                continue

            match = ANNOTATION.fullmatch(annotation)
            if not match:
                errors.append(
                    f"line {line_number}: malformed fastENLOC annotation {annotation!r}"
                )
                continue
            set_id = match.group("set_id")
            if PLACEHOLDER_SET_ID.search(set_id):
                errors.append(
                    f"line {line_number}: placeholder credible-set ID {set_id!r}"
                )
            try:
                pip = float(match.group("pip"))
                cs_pip = float(match.group("cs_pip"))
                set_size = int(match.group("set_size"))
            except ValueError:
                errors.append(f"line {line_number}: non-numeric annotation value")
                continue
            if not math.isfinite(pip) or not 0 <= pip <= 1:
                errors.append(f"line {line_number}: invalid variant PIP {pip}")
            if not math.isfinite(cs_pip) or not 0 < cs_pip <= 1.000001:
                errors.append(f"line {line_number}: invalid credible-set PIP {cs_pip}")
            if set_size < 1:
                errors.append(f"line {line_number}: invalid declared set size {set_size}")

            if set_id in declared_sizes and declared_sizes[set_id] != set_size:
                errors.append(
                    f"line {line_number}: inconsistent size for {set_id}: "
                    f"{declared_sizes[set_id]} vs {set_size}"
                )
            if (
                set_id in declared_cs_pip
                and abs(declared_cs_pip[set_id] - cs_pip) > 1e-8
            ):
                errors.append(
                    f"line {line_number}: inconsistent cumulative PIP for {set_id}"
                )
            if variant_id in variants[set_id]:
                errors.append(
                    f"line {line_number}: duplicate variant {variant_id!r} in {set_id}"
                )

            declared_sizes[set_id] = set_size
            declared_cs_pip[set_id] = cs_pip
            row_counts[set_id] += 1
            pip_sums[set_id] += pip
            chromosomes[set_id].add(chrom.removeprefix("chr"))
            variants[set_id].add(variant_id)

    if n_rows == 0:
        errors.append("GWAS file has no data rows")
    observed = len(row_counts)
    if observed != expected_credible_sets:
        errors.append(
            f"{study_id}: observed {observed} distinct credible-set IDs, "
            f"expected {expected_credible_sets} from the manifest"
        )
    for set_id in sorted(row_counts):
        if row_counts[set_id] != declared_sizes[set_id]:
            errors.append(
                f"{set_id}: observed {row_counts[set_id]} rows but annotation "
                f"declares {declared_sizes[set_id]}"
            )
        if len(chromosomes[set_id]) != 1:
            errors.append(
                f"{set_id}: variants span chromosomes "
                f"{','.join(sorted(chromosomes[set_id]))}"
            )
        tolerance = max(5e-4, 0.005 * declared_cs_pip[set_id])
        if abs(pip_sums[set_id] - declared_cs_pip[set_id]) > tolerance:
            errors.append(
                f"{set_id}: variant PIPs sum to {pip_sums[set_id]:.8g}, "
                f"annotation declares {declared_cs_pip[set_id]:.8g}"
            )

    if errors:
        preview = "\n  - ".join(errors[:25])
        remainder = len(errors) - min(len(errors), 25)
        suffix = f"\n  ... and {remainder} more" if remainder else ""
        raise ValueError(
            f"Invalid fastENLOC GWAS input for {study_id}:\n  - {preview}{suffix}"
        )

    return {
        "study_id": study_id,
        "status": "PASS",
        "n_fastenloc_rows": n_rows,
        "expected_credible_sets": expected_credible_sets,
        "observed_credible_sets": observed,
        "min_set_size": min(row_counts.values()),
        "max_set_size": max(row_counts.values()),
        "min_credible_set_pip_sum": f"{min(pip_sums.values()):.8g}",
        "max_credible_set_pip_sum": f"{max(pip_sums.values()):.8g}",
    }


def main() -> None:
    args = parse_args()
    result = validate(args.gwas, args.study_id, args.expected_credible_sets)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(result), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(result)


if __name__ == "__main__":
    main()
