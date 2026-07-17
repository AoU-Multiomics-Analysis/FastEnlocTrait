from __future__ import annotations

import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "pull_finemapping", ROOT / "pull_finemapping.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def locus_row(variant_id: str, position: int, pip: float) -> dict:
    return {
        "posteriorProbability": pip,
        "is95CredibleSet": True,
        "is99CredibleSet": True,
        "variant": {
            "id": variant_id,
            "chromosome": "1",
            "position": position,
            "referenceAllele": "A",
            "alternateAllele": "G",
        },
    }


class PullFineMappingTests(unittest.TestCase):
    def test_checked_in_manifest_is_redacted_and_filtered(self) -> None:
        manifest = ROOT / "manifests" / "opentargets_studies.no_paths.tsv"
        with manifest.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            headers = reader.fieldnames or []
        self.assertEqual(len(rows), 401)
        self.assertNotIn("gwas_path", headers)
        self.assertEqual(len({row["study_id"] for row in rows}), 401)
        flagged = [
            row["trait"]
            for row in rows
            if "(mtag)" in row["trait"].lower()
            or "pleiotrop" in row["trait"].lower()
            or " or " in row["trait"].lower()
        ]
        self.assertEqual(flagged, [])

    def test_build_study_file_filters_method_and_writes_fastenloc(self) -> None:
        config = {
            "study_id": "TEST001",
            "trait": "Example disease",
            "trait_category": "test",
            "finemap_method": "SuSiE",
            "n_credible_sets": "1",
            "n_variants": "1000",
        }
        credible_sets = [
            {
                "studyId": "TEST001",
                "studyLocusId": "L1",
                "chromosome": "1",
                "position": 100,
                "finemappingMethod": "SuSie",
                "locus": {
                    "count": 2,
                    "rows": [
                        locus_row("1_100_A_G", 100, 0.60),
                        locus_row("1_101_A_G", 101, 0.35),
                    ],
                },
            },
            {
                "studyId": "TEST001",
                "studyLocusId": "PICS1",
                "chromosome": "1",
                "position": 200,
                "finemappingMethod": "PICS",
                "locus": {
                    "count": 2,
                    "rows": [
                        locus_row("1_200_A_G", 200, 0.70),
                        locus_row("1_201_A_G", 201, 0.25),
                    ],
                },
            },
        ]
        with tempfile.TemporaryDirectory() as temp:
            manifest_row, qc = MODULE.build_study_file(
                config,
                credible_sets,
                Path(temp),
                pip_sum_min=0.90,
                allow_count_mismatch=False,
            )
            output = Path(manifest_row["gwas_path"])
            self.assertTrue(output.exists())
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                lines = handle.readlines()
            self.assertEqual(len(lines), 2)
            self.assertTrue(lines[0].startswith("chr1\t100\tchr1_100_A_G\tA\tG\t"))
            self.assertIn("Example_disease;TEST001_L1=", lines[0])
            self.assertEqual(qc["n_credible_sets_written"], 1)
            self.assertEqual(qc["n_fastenloc_rows"], 2)

    def test_count_mismatch_fails_closed(self) -> None:
        config = {
            "study_id": "TEST002",
            "trait": "Example disease",
            "trait_category": "test",
            "finemap_method": "SuSiE",
            "n_credible_sets": "2",
            "n_variants": "1000",
        }
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(RuntimeError, "expected 2"):
                MODULE.build_study_file(
                    config,
                    [],
                    Path(temp),
                    pip_sum_min=0.90,
                    allow_count_mismatch=False,
                )


if __name__ == "__main__":
    unittest.main()
