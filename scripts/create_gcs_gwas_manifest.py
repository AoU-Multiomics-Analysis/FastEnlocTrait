#!/usr/bin/env python3
"""Create a cloud-ready GWAS manifest from an organized local release."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--release-root", required=True, type=Path)
    parser.add_argument("--gcs-prefix", required=True)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    release_root = args.release_root.resolve()
    prefix = args.gcs_prefix.rstrip("/")
    if not prefix.startswith("gs://"):
        raise RuntimeError("--gcs-prefix must start with gs://")

    with args.manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"{args.manifest}: missing header")
        headers = list(reader.fieldnames)
        rows = list(reader)
    if "gwas_path" not in headers:
        raise RuntimeError(f"{args.manifest}: missing gwas_path column")
    if not rows:
        raise RuntimeError(f"{args.manifest}: contains no studies")

    for row in rows:
        local_path = Path(row["gwas_path"]).resolve()
        if not local_path.is_file():
            raise RuntimeError(
                f"{row['study_id']}: local GWAS file does not exist: {local_path}"
            )
        try:
            relative = local_path.relative_to(release_root)
        except ValueError as exc:
            raise RuntimeError(
                f"{row['study_id']}: GWAS file is outside the release root"
            ) from exc
        row["gwas_path"] = f"{prefix}/{relative.as_posix()}"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(args.out.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=headers,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.out)
    print(f"Wrote {len(rows)} cloud paths to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
