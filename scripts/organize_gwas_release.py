#!/usr/bin/env python3
"""Reorganize a GWAS release using paths from a prior manifest."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path, PurePosixPath


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--layout-manifest", required=True, type=Path)
    parser.add_argument("--release-root", required=True, type=Path)
    parser.add_argument("--retrieval-qc", type=Path)
    parser.add_argument(
        "--layout-marker",
        default="GWASColocDataV2/",
        help="Keep the portion of each prior gwas_path after this marker.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"{path}: missing header")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
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
    temporary.replace(path)


def relative_layout(path: str, marker: str) -> Path:
    if marker not in path:
        raise RuntimeError(f"Layout path does not contain {marker!r}: {path}")
    relative = PurePosixPath(path.split(marker, 1)[1])
    if relative.is_absolute() or ".." in relative.parts or len(relative.parts) < 2:
        raise RuntimeError(f"Unsafe or unorganized relative layout: {relative}")
    return Path(*relative.parts)


def main() -> int:
    args = parse_args()
    headers, rows = read_tsv(args.manifest)
    _, layout_rows = read_tsv(args.layout_manifest)
    layouts = {
        row["study_id"]: relative_layout(row["gwas_path"], args.layout_marker)
        for row in layout_rows
    }
    if len(layouts) != len(layout_rows):
        raise RuntimeError(f"{args.layout_manifest}: duplicate study IDs")

    release_root = args.release_root.resolve()
    moves: list[tuple[Path, Path, dict[str, str]]] = []
    for row in rows:
        study_id = row["study_id"]
        if study_id not in layouts:
            raise RuntimeError(f"{study_id}: absent from layout manifest")
        source = Path(row["gwas_path"]).resolve()
        destination = (release_root / layouts[study_id]).resolve()
        if release_root not in destination.parents:
            raise RuntimeError(f"{study_id}: destination escapes release root")
        if not source.is_file():
            raise RuntimeError(f"{study_id}: source file does not exist: {source}")
        if destination.exists() and destination != source:
            raise RuntimeError(f"{study_id}: destination already exists: {destination}")
        moves.append((source, destination, row))

    updated_paths: dict[str, str] = {}
    for source, destination, row in moves:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source != destination:
            shutil.move(source, destination)
        row["gwas_path"] = str(destination)
        updated_paths[row["study_id"]] = str(destination)
    write_tsv(args.manifest, headers, rows)

    if args.retrieval_qc:
        qc_headers, qc_rows = read_tsv(args.retrieval_qc)
        for row in qc_rows:
            if row["study_id"] in updated_paths:
                row["gwas_path"] = updated_paths[row["study_id"]]
        write_tsv(args.retrieval_qc, qc_headers, qc_rows)

    old_gwas_dir = release_root / "gwas"
    if old_gwas_dir.is_dir() and not any(old_gwas_dir.iterdir()):
        old_gwas_dir.rmdir()
    print(f"Organized {len(rows)} studies under {release_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
