import argparse
import tempfile
import unittest
from pathlib import Path

import build_variant_parquet as bvp
import score_blocks


def make_boundaries_file(path: Path, pairs: list[tuple[int, int]]) -> None:
    lines = ["START\tEND\n"] + [f"{s}\t{e}\n" for s, e in pairs]
    path.write_text("".join(lines))


class LoadAndVerifyPilotBlocksTest(unittest.TestCase):
    def test_returns_sorted_blocks_when_all_pairs_present(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "boundaries.tsv"
            # shuffled order, plus unrelated extra rows -- should still work
            pairs = list(reversed(bvp.PILOT_BLOCKS_CHR22)) + [(1, 2), (3, 4)]
            make_boundaries_file(path, pairs)

            blocks = bvp.load_and_verify_pilot_blocks(path)

            self.assertEqual(len(blocks), 12)
            starts = [b.start for b in blocks]
            self.assertEqual(starts, sorted(starts))
            first = blocks[0]
            self.assertEqual(first.chrom, "chr22")
            self.assertEqual(first.block_id, f"chr22_{first.start}-{first.end}")

    def test_raises_when_a_pilot_pair_is_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "boundaries.tsv"
            make_boundaries_file(path, bvp.PILOT_BLOCKS_CHR22[1:])  # drop the first pair

            with self.assertRaisesRegex(ValueError, "not found in upstream boundaries"):
                bvp.load_and_verify_pilot_blocks(path)


class ParseBcftoolsRowTest(unittest.TestCase):
    def test_parses_a_well_formed_snv_row(self):
        variant = bvp.parse_bcftools_row("22\t23507496\tC\tT\t110\t5096\t0.02\tSNP\n")
        self.assertEqual(variant.chrom, "chr22")
        self.assertEqual(variant.pos, 23507496)
        self.assertEqual(variant.ref, "C")
        self.assertEqual(variant.alt, "T")
        self.assertEqual(variant.ac, 110)
        self.assertEqual(variant.an, 5096)
        self.assertAlmostEqual(variant.af, 0.02)
        self.assertEqual(variant.vt, "SNP")

    def test_rejects_af_outside_zero_one(self):
        with self.assertRaisesRegex(ValueError, "AF=1.5 outside"):
            bvp.parse_bcftools_row("22\t100\tC\tT\t1\t2\t1.5\tSNP\n")

    def test_rejects_multiallelic_alt(self):
        with self.assertRaisesRegex(ValueError, "multiallelic"):
            bvp.parse_bcftools_row("22\t100\tC\tT,G\t1,1\t2\t0.5,0.5\tSNP\n")

    def test_rejects_ac_greater_than_an(self):
        with self.assertRaisesRegex(ValueError, "implausible AC"):
            bvp.parse_bcftools_row("22\t100\tC\tT\t10\t2\t0.5\tSNP\n")


class AssignAndDedupeTest(unittest.TestCase):
    def _variant(self, pos, ref="A", alt="G", vt="SNP"):
        return bvp.RawVariant("chr22", pos, ref, alt, ac=1, an=2, af=0.5, vt=vt)

    def test_assigns_boundary_position_to_the_next_block_not_the_previous(self):
        # Half-open [start, end): position == end belongs to the *next* block.
        blocks = [
            score_blocks.Block("chr22", 100, 200, "chr22_100-200"),
            score_blocks.Block("chr22", 200, 300, "chr22_200-300"),
        ]
        variants = [
            self._variant(pos=200),  # 1-based VCF pos 200 -> 0-based 199 -> first block
            self._variant(pos=201),  # 0-based 200 -> second block, exactly at its start
        ]

        rows, duplicates = bvp.assign_and_dedupe(variants, blocks)

        self.assertEqual(duplicates, 0)
        self.assertEqual(rows[0]["block_id"], "chr22_100-200")
        self.assertEqual(rows[1]["block_id"], "chr22_200-300")

    def test_every_variant_lands_in_exactly_one_block(self):
        blocks = [
            score_blocks.Block("chr22", 0, 10, "chr22_0-10"),
            score_blocks.Block("chr22", 20, 30, "chr22_20-30"),  # gap between blocks
        ]
        variants = [self._variant(pos=p) for p in (1, 5, 10, 21, 29, 30)]
        # pos 10 (0-based 9) -> first block; pos 30 (0-based 29) -> second block

        rows, _ = bvp.assign_and_dedupe(variants, blocks)

        self.assertEqual(len(rows), len(variants))
        for row in rows:
            self.assertIn(row["block_id"], {"chr22_0-10", "chr22_20-30"})

    def test_variant_outside_every_block_raises(self):
        blocks = [score_blocks.Block("chr22", 0, 10, "chr22_0-10")]
        with self.assertRaisesRegex(ValueError, "does not fall inside any pilot block"):
            bvp.assign_and_dedupe([self._variant(pos=50)], blocks)

    def test_exact_duplicate_is_dropped_but_same_position_different_allele_is_not(self):
        blocks = [score_blocks.Block("chr22", 0, 100, "chr22_0-100")]
        variants = [
            self._variant(pos=5, ref="A", alt="G"),
            self._variant(pos=5, ref="A", alt="G"),  # exact duplicate
            self._variant(pos=5, ref="A", alt="C"),  # same position, different allele
        ]

        rows, duplicates = bvp.assign_and_dedupe(variants, blocks)

        self.assertEqual(duplicates, 1)
        self.assertEqual(len(rows), 2)
        alts = sorted(row["alt"] for row in rows)
        self.assertEqual(alts, ["C", "G"])


class AviJoinTest(unittest.TestCase):
    def test_join_matches_by_chrom_pos_alt_and_counts_unmatched(self):
        avi_rows = [
            "chr22\t100\tA\tG\t0.9\t35\n",
            "chr22\t200\tA\tT\t0.1\t5\n",
        ]
        index = bvp.build_avi_index(
            avi_rows,
            chrom_column=0,
            position_column=1,
            alt_column=3,
            raw_column=4,
            phred_column=5,
        )
        variant_rows = [
            {"chrom": "chr22", "pos": 100, "ref": "A", "alt": "G"},
            {"chrom": "chr22", "pos": 300, "ref": "A", "alt": "G"},  # no AVI row
        ]

        unmatched = bvp.join_avi(variant_rows, index)

        self.assertEqual(unmatched, 1)
        self.assertAlmostEqual(variant_rows[0]["avi_raw"], 0.9)
        self.assertAlmostEqual(variant_rows[0]["avi_phred"], 35.0)
        self.assertIsNone(variant_rows[1]["avi_raw"])
        self.assertIsNone(variant_rows[1]["avi_phred"])

    def test_join_ignores_ref_when_matching(self):
        # AVI rows carry no REF column; a REF discrepancy must not block the join.
        avi_rows = ["chr22\t100\tA\tG\t0.9\t35\n"]
        index = bvp.build_avi_index(
            avi_rows,
            chrom_column=0,
            position_column=1,
            alt_column=3,
            raw_column=4,
            phred_column=5,
        )
        variant_rows = [{"chrom": "chr22", "pos": 100, "ref": "C", "alt": "G"}]

        unmatched = bvp.join_avi(variant_rows, index)

        self.assertEqual(unmatched, 0)
        self.assertAlmostEqual(variant_rows[0]["avi_raw"], 0.9)


class ResolveAviColumnsTest(unittest.TestCase):
    def test_resolves_named_header_columns(self):
        args = argparse.Namespace(
            avi_chrom_column=None,
            avi_position_column=None,
            avi_alt_column=None,
            avi_raw_column=None,
            avi_phred_column=None,
        )
        header = ["chrom", "position", "ref", "alt", "avi_raw", "avi_phred"]
        self.assertEqual(bvp.resolve_avi_columns(args, header), (0, 1, 3, 4, 5))

    def test_headerless_input_requires_explicit_alt_column(self):
        args = argparse.Namespace(
            avi_chrom_column=None,
            avi_position_column=None,
            avi_alt_column=None,
            avi_raw_column="4",
            avi_phred_column="5",
        )
        with self.assertRaisesRegex(ValueError, "could not infer the avi-alt column"):
            bvp.resolve_avi_columns(args, None)


if __name__ == "__main__":
    unittest.main()
