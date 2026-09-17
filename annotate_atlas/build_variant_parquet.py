#!/usr/bin/env python3
"""Build a per-SNV Parquet table for the chr22 pilot haploblocks.

Extracts observed SNVs from the phased 1000G VCF inside the 12 chr22 pilot
blocks (https://data.haploblocks.org/chr22_small_dataset/), assigns each to
exactly one haploblock, and joins the AlphaGenome Atlas AVI score by
(chromosome, position, ALT). Indels in the same regions are written to a
separate Parquet file with no AVI columns -- the static AVI download is
SNVs only.

Boundary convention (verified, not assumed): chr22_haploblock_boundaries_chr22.tsv
gives 0-based, half-open [start, end) intervals. Evidence:
  - Adjacent blocks share exactly one coordinate (block[i].end ==
    block[i+1].start); a 1-based-inclusive reading would double-assign that
    position to both blocks, which contradicts "assign each SNV to exactly
    one haploblock".
  - block_stats.tsv (upstream QC) reports block_length == end - start
    exactly for every block, e.g. chr22_44303006-44307622 -> block_length
    4616 == 44307622 - 44303006.
  - annotate_atlas/README already documents "half-open block interval", and
    score_blocks.py / build_blocks_bed.py already implement blocks this way.

Usage:
  python3 build_variant_parquet.py VCF BOUNDARIES_TSV AVI_PATH [--output OUT.parquet]

AVI_PATH: local path or HTTP(S) URL to the tabix-indexed AVI SNV scores
(needs a matching .tbi). Column names are resolved from the real file
header, or via --avi-*-column if the file is headerless -- never guessed.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

from score_blocks import (
    Block,
    canonical_chrom,
    detect_contig_style,
    parse_column_spec,
    read_tabix_contigs,
    read_tabix_header,
    stream_tabix_rows,
    write_tabix_regions,
)

# The 12 chr22 pilot blocks published at
# https://data.haploblocks.org/chr22_small_dataset/. Verified against
# chr22_haploblock_boundaries_chr22.tsv: every pair below is an exact
# (START, END) row in that file (see load_and_verify_pilot_blocks).
PILOT_BLOCKS_CHR22: list[tuple[int, int]] = [
    (23507389, 23529638),
    (26245115, 26270093),
    (27020234, 27031341),
    (27396539, 27435114),
    (30383333, 30440683),
    (33827023, 33873503),
    (36238986, 36271022),
    (37075250, 37101553),
    (44303006, 44307622),
    (44791710, 44821141),
    (46418842, 46438932),
    (49264062, 49281800),
]

PILOT_CHROM = "chr22"  # canonical / output chrom
DEFAULT_OUTPUT = Path("/dcai/projects/iu_0142/team11/data/parquet/chr22_small_variants_avi.parquet")

BCFTOOLS_QUERY_FORMAT = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/VT\n"

SNV_COLUMNS = [
    "chrom", "pos", "ref", "alt", "variant_id",
    "AC", "AN", "AF",
    "block_id", "block_start", "block_end",
    "avi_raw", "avi_phred",
]
INDEL_COLUMNS = [c for c in SNV_COLUMNS if c not in ("avi_raw", "avi_phred")]

AVI_CHROM_ALIASES = ("chrom", "chr", "chromosome", "seqname")
AVI_POSITION_ALIASES = ("pos", "position", "bp")
AVI_ALT_ALIASES = ("alt", "allele", "alternate", "alt_allele")
AVI_RAW_ALIASES = ("avi_raw", "avi_score", "raw", "raw_score", "score")
AVI_PHRED_ALIASES = ("avi_phred", "phred", "phred_score", "avi_phred_score")


@dataclass(frozen=True)
class RawVariant:
    chrom: str  # canonical, e.g. "chr22"
    pos: int  # 1-based VCF position
    ref: str
    alt: str
    ac: int
    an: int
    af: float
    vt: str


def pilot_blocks() -> list[Block]:
    return [
        Block(PILOT_CHROM, start, end, f"{PILOT_CHROM}_{start}-{end}")
        for start, end in sorted(PILOT_BLOCKS_CHR22)
    ]


def load_and_verify_pilot_blocks(boundaries_tsv: Path) -> list[Block]:
    """Cross-check the hardcoded pilot coordinates against the upstream boundaries file."""
    with boundaries_tsv.open() as stream:
        header = stream.readline().rstrip("\n").split("\t")
        if header != ["START", "END"]:
            raise ValueError(f"{boundaries_tsv}: unexpected header {header!r}")
        known: set[tuple[int, int]] = set()
        for line_number, line in enumerate(stream, 2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{boundaries_tsv}:{line_number}: expected 2 columns")
            known.add((int(fields[0]), int(fields[1])))

    missing = [pair for pair in PILOT_BLOCKS_CHR22 if pair not in known]
    if missing:
        raise ValueError(
            f"{boundaries_tsv}: pilot blocks not found in upstream boundaries: {missing}"
        )
    return pilot_blocks()


def parse_bcftools_row(line: str) -> RawVariant:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 8:
        raise ValueError(f"expected 8 bcftools query columns, got {len(fields)}: {line[:200]!r}")
    chrom_text, pos_text, ref, alt, ac_text, an_text, af_text, vt = fields

    for label, text in (("ALT", alt), ("AC", ac_text), ("AF", af_text)):
        if "," in text:
            raise ValueError(
                f"{label}={text!r} at {chrom_text}:{pos_text} looks multiallelic; "
                "run `bcftools norm -m -any` on the VCF first"
            )

    try:
        pos = int(pos_text)
        ac = int(ac_text)
        an = int(an_text)
        af = float(af_text)
    except ValueError as error:
        raise ValueError(
            f"non-numeric POS/AC/AN/AF at {chrom_text}:{pos_text}: {line[:200]!r}"
        ) from error

    if not (0.0 <= af <= 1.0):
        raise ValueError(f"AF={af} outside [0, 1] at {chrom_text}:{pos_text}")
    if an <= 0 or ac < 0 or ac > an:
        raise ValueError(f"implausible AC={ac}/AN={an} at {chrom_text}:{pos_text}")

    return RawVariant(canonical_chrom(chrom_text), pos, ref, alt, ac, an, af, vt)


def assign_and_dedupe(
    raw_variants: Iterable[RawVariant], blocks: Sequence[Block]
) -> tuple[list[dict], int]:
    """Assign each variant to its one pilot block; drop exact (chrom,pos,ref,alt) duplicates.

    `raw_variants` must be position-sorted per chrom (true of bcftools query
    output) and every row must fall inside one of `blocks`, since the VCF
    query is itself region-restricted to those exact half-open intervals.
    """
    cursor = 0
    seen: set[tuple[str, int, str, str]] = set()
    rows: list[dict] = []
    duplicates = 0

    for variant in raw_variants:
        key = (variant.chrom, variant.pos, variant.ref, variant.alt)
        if key in seen:
            duplicates += 1
            continue
        seen.add(key)

        position_zero_based = variant.pos - 1
        while cursor < len(blocks) and blocks[cursor].end <= position_zero_based:
            cursor += 1
        if cursor >= len(blocks) or not (
            blocks[cursor].start <= position_zero_based < blocks[cursor].end
        ):
            raise ValueError(
                f"{variant.chrom}:{variant.pos} does not fall inside any pilot block; "
                "the VCF query region and block table have diverged"
            )
        block = blocks[cursor]

        rows.append(
            {
                "chrom": variant.chrom,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "variant_id": f"{variant.chrom}-{variant.pos}-{variant.ref}-{variant.alt}",
                "AC": variant.ac,
                "AN": variant.an,
                "AF": variant.af,
                "block_id": block.block_id,
                "block_start": block.start,
                "block_end": block.end,
            }
        )

    return rows, duplicates


def stream_bcftools_query(bcftools: str, vcf_path: Path, regions_path: Path) -> Iterator[str]:
    try:
        process = subprocess.Popen(
            [bcftools, "query", "-R", str(regions_path), "-f", BCFTOOLS_QUERY_FORMAT, str(vcf_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(f"bcftools executable not found: {bcftools!r}") from None

    assert process.stdout is not None
    for line in process.stdout:
        if line.strip():
            yield line
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bcftools query failed ({return_code}): {stderr.strip()}")


def resolve_avi_columns(
    args: argparse.Namespace, header: Sequence[str] | None
) -> tuple[int, int, int, int, int]:
    chrom_column = parse_column_spec(
        args.avi_chrom_column, header, AVI_CHROM_ALIASES, "avi-chrom", fallback=0
    )
    position_column = parse_column_spec(
        args.avi_position_column, header, AVI_POSITION_ALIASES, "avi-position", fallback=1
    )
    alt_column = parse_column_spec(args.avi_alt_column, header, AVI_ALT_ALIASES, "avi-alt")
    raw_column = parse_column_spec(args.avi_raw_column, header, AVI_RAW_ALIASES, "avi-raw")
    phred_column = parse_column_spec(args.avi_phred_column, header, AVI_PHRED_ALIASES, "avi-phred")
    return chrom_column, position_column, alt_column, raw_column, phred_column


def parse_avi_row(
    line: str,
    *,
    chrom_column: int,
    position_column: int,
    alt_column: int,
    raw_column: int,
    phred_column: int,
) -> tuple[tuple[str, int, str], float | None, float | None]:
    fields = line.split("\t")
    required = max(chrom_column, position_column, alt_column, raw_column, phred_column)
    if len(fields) <= required:
        raise ValueError(
            f"AVI row has {len(fields)} columns; need column {required + 1}: {line[:200]!r}"
        )
    chrom = canonical_chrom(fields[chrom_column])
    position = int(fields[position_column])
    alt = fields[alt_column].strip()

    def parse_score(text: str) -> float | None:
        text = text.strip()
        return float(text) if text not in {"", ".", "NA", "NaN"} else None

    avi_raw = parse_score(fields[raw_column])
    avi_phred = parse_score(fields[phred_column])
    return (chrom, position, alt), avi_raw, avi_phred


def build_avi_index(
    rows: Iterable[str], **column_kwargs: int
) -> dict[tuple[str, int, str], tuple[float | None, float | None]]:
    index: dict[tuple[str, int, str], tuple[float | None, float | None]] = {}
    for line in rows:
        key, avi_raw, avi_phred = parse_avi_row(line, **column_kwargs)
        index[key] = (avi_raw, avi_phred)
    return index


def join_avi(
    variant_rows: list[dict],
    avi_index: dict[tuple[str, int, str], tuple[float | None, float | None]],
) -> int:
    """Join avi_raw/avi_phred onto variant_rows in place; return the unmatched count."""
    unmatched = 0
    for row in variant_rows:
        hit = avi_index.get((row["chrom"], row["pos"], row["alt"]))
        if hit is None:
            unmatched += 1
            row["avi_raw"], row["avi_phred"] = None, None
        else:
            row["avi_raw"], row["avi_phred"] = hit
    return unmatched


def write_parquet(rows: list[dict], columns: Sequence[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table({column: [row[column] for row in rows] for column in columns})
    pq.write_table(table, path, compression="zstd")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf_path", type=Path)
    parser.add_argument("boundaries_tsv", type=Path)
    parser.add_argument("avi_path", help="local path or HTTP(S) URL to AVI SNV .tsv.gz")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--indels-output",
        type=Path,
        default=None,
        help="default: --output with _avi.parquet replaced by _indels.parquet",
    )
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--avi-chrom-column", help="AVI chrom column name or 1-based number")
    parser.add_argument("--avi-position-column", help="AVI position column name or 1-based number")
    parser.add_argument("--avi-alt-column", help="AVI ALT column name or 1-based number")
    parser.add_argument("--avi-raw-column", help="AVI raw-score column name or 1-based number")
    parser.add_argument("--avi-phred-column", help="AVI PHRED column name or 1-based number")
    parser.add_argument(
        "--contig-style",
        choices=("auto", "chr", "bare"),
        default="auto",
        help="contig naming in the AVI index (default: detect with tabix -l)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    os.umask(0o002)  # keep files group-writable (cabale/iu_0142 project share)
    args = build_parser().parse_args(argv)
    indels_output = args.indels_output or Path(
        str(args.output).replace("_avi.parquet", "_indels.parquet")
    )

    blocks = load_and_verify_pilot_blocks(args.boundaries_tsv)
    vcf_regions_path = write_tabix_regions(blocks, "bare")  # this VCF's contigs are "22", not "chr22"
    try:
        raw_snvs: list[RawVariant] = []
        raw_indels: list[RawVariant] = []
        for line in stream_bcftools_query(args.bcftools, args.vcf_path, vcf_regions_path):
            variant = parse_bcftools_row(line)
            if variant.vt == "SNP":
                if len(variant.ref) != 1 or len(variant.alt) != 1:
                    raise ValueError(
                        f"VT=SNP but REF/ALT not single-base at "
                        f"{variant.chrom}:{variant.pos} ({variant.ref}>{variant.alt})"
                    )
                raw_snvs.append(variant)
            elif variant.vt == "INDEL":
                raw_indels.append(variant)
            else:
                raise ValueError(f"unhandled VT={variant.vt!r} at {variant.chrom}:{variant.pos}")

        snv_rows, snv_duplicates = assign_and_dedupe(raw_snvs, blocks)
        indel_rows, indel_duplicates = assign_and_dedupe(raw_indels, blocks)
    finally:
        vcf_regions_path.unlink(missing_ok=True)

    avi_header = read_tabix_header(args.tabix, args.avi_path)
    contig_style = args.contig_style
    if contig_style == "auto":
        contig_style = detect_contig_style(read_tabix_contigs(args.tabix, args.avi_path))
    chrom_column, position_column, alt_column, raw_column, phred_column = resolve_avi_columns(
        args, avi_header
    )
    avi_regions_path = write_tabix_regions(blocks, contig_style)
    try:
        avi_index = build_avi_index(
            stream_tabix_rows(args.tabix, args.avi_path, avi_regions_path),
            chrom_column=chrom_column,
            position_column=position_column,
            alt_column=alt_column,
            raw_column=raw_column,
            phred_column=phred_column,
        )
    finally:
        avi_regions_path.unlink(missing_ok=True)

    unmatched = join_avi(snv_rows, avi_index)

    write_parquet(snv_rows, SNV_COLUMNS, args.output)
    write_parquet(indel_rows, INDEL_COLUMNS, indels_output)

    per_block = Counter(row["block_id"] for row in snv_rows)
    print(f"wrote {args.output}: {len(snv_rows)} SNV rows", file=sys.stderr)
    print(f"wrote {indels_output}: {len(indel_rows)} indel rows", file=sys.stderr)
    print(f"unmatched AVI: {unmatched}", file=sys.stderr)
    print(f"duplicate SNV rows dropped: {snv_duplicates}", file=sys.stderr)
    print(f"duplicate indel rows dropped: {indel_duplicates}", file=sys.stderr)
    print("SNVs per block:", file=sys.stderr)
    for block in blocks:
        print(f"  {block.block_id}: {per_block.get(block.block_id, 0)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
