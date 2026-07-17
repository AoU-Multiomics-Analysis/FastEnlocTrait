from __future__ import annotations

import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "validate_fastenloc_gwas", ROOT / "scripts" / "validate_fastenloc_gwas.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_gwas(path: Path, lines: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.writelines(line + "\n" for line in lines)


class ValidateFastenlocGWASTests(unittest.TestCase):
    def test_valid_file_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "good.tsv.gz"
            write_gwas(
                path,
                [
                    "chr1\t100\tchr1_100_A_G\tA\tG\ttrait;S1_chr1.100.101_LL1=6.0000e-01[9.500e-01:2]",
                    "chr1\t101\tchr1_101_A_G\tA\tG\ttrait;S1_chr1.100.101_LL1=3.5000e-01[9.500e-01:2]",
                ],
            )
            result = MODULE.validate(path, "S1", 1)
            self.assertEqual(result["observed_credible_sets"], 1)
            self.assertEqual(result["status"], "PASS")

    def test_placeholder_lnone_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bad.tsv.gz"
            write_gwas(
                path,
                [
                    "chr1\t100\tchr1_100_A_G\tA\tG\ttrait;S1_chr1.0.0_LNone=9.5000e-01[9.500e-01:1]"
                ],
            )
            with self.assertRaisesRegex(ValueError, "placeholder credible-set ID"):
                MODULE.validate(path, "S1", 1)

    def test_manifest_count_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bad.tsv.gz"
            write_gwas(
                path,
                [
                    "chr1\t100\tchr1_100_A_G\tA\tG\ttrait;S1_chr1.100.100_LL1=9.5000e-01[9.500e-01:1]"
                ],
            )
            with self.assertRaisesRegex(ValueError, "observed 1.*expected 2"):
                MODULE.validate(path, "S1", 2)

    def test_declared_set_size_mismatch_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bad.tsv.gz"
            write_gwas(
                path,
                [
                    "chr1\t100\tchr1_100_A_G\tA\tG\ttrait;S1_chr1.100.100_LL1=9.5000e-01[9.500e-01:2]"
                ],
            )
            with self.assertRaisesRegex(ValueError, "observed 1 rows"):
                MODULE.validate(path, "S1", 1)


if __name__ == "__main__":
    unittest.main()
