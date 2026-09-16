#!/usr/bin/env python3
"""Annotate haploblocks with AlphaGenome Atlas AVI scores via tabix.

Usage:
  python3 score_blocks.py BLOCKS_BED AVI_PATH_OR_URL [OUT_TSV] [--inspect]

BLOCKS_BED: output of build_blocks_bed.py (chrom, start, end, name; 0-based
half-open).
AVI_PATH_OR_URL: local path or http(s) URL to the tabix-indexed avi.tsv.gz
(needs a matching .tbi alongside it, local or remote). tabix handles both
local files and HTTP(S) range-queryable URLs the same way.
--inspect: print the first block's raw tabix rows and exit, instead of
scoring everything. Run this first against real data to confirm the two
TODO column indices below.

Rules baked in (see hackathon task spec):
- AVI is unsigned PHRED: pooled directly, no sign handling.
- PHRED 20 == genome-wide top 1% of predictions (per AlphaGenome Atlas
  paper), so avi_top1 / avi_dens filter on PHRED >= 20 rather than
  re-deriving a per-block rank cutoff.
- SNVs only: the Atlas static AVI download is genome-wide SNVs only as of
  this writing (no indels), so every tabix row is already an SNV -- no
  extra REF/ALT length filtering needed.
"""
import subprocess
import sys
import time
from pathlib import Path

# TODO(markus): confirm real avi.tsv.gz column layout once you've unzipped
# avi_scores_snvs_tabix.zip on the cluster VM (avi.tsv.gz + avi.tsv.gz.tbi).
# The Atlas paper (Fig 1B) says AVI is exposed as Raw / Percentile / PHRED,
# but the exact column count/order of the static tabix file isn't
# documented anywhere I could find (paper, readthedocs, downloads page).
# Run with --inspect first and fix these two indices before trusting any
# real output.
COL_POS = 1  # 0-indexed column holding variant position
COL_PHRED = -1  # 0-indexed column holding unsigned PHRED-scaled AVI score

PHRED_THRESHOLD = 20  # paper: PHRED 20 == genome-wide top 1% of predictions


def tabix_query(avi_path: str, chrom: str, start: int, end: int) -> list[str]:
    # BED is 0-based half-open; tabix regions are 1-based inclusive.
    region = f"{chrom}:{start + 1}-{end}"
    result = subprocess.run(
        ["tabix", avi_path, region],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def score_block(avi_path: str, chrom: str, start: int, end: int) -> dict:
    rows = tabix_query(avi_path, chrom, start, end)
    n_snv = len(rows)
    length_kb = (end - start) / 1000.0

    hits = []
    for line in rows:
        fields = line.split("\t")
        try:
            phred = float(fields[COL_PHRED])
        except (IndexError, ValueError):
            continue
        if phred >= PHRED_THRESHOLD:
            hits.append(phred)

    avi_top1 = sum(hits) / len(hits) if hits else ""
    avi_dens = len(hits) / length_kb if length_kb > 0 else ""

    return {
        "n_snv": n_snv,
        "avi_top1": avi_top1,
        "avi_dens": avi_dens,
        "length_kb": length_kb,
    }


def inspect(avi_path: str, chrom: str, start: int, end: int) -> None:
    rows = tabix_query(avi_path, chrom, start, end)
    print(f"{len(rows)} rows for {chrom}:{start}-{end}", file=sys.stderr)
    for line in rows[:5]:
        fields = line.split("\t")
        print(f"  {len(fields)} cols: {fields}", file=sys.stderr)
    print(
        f"\nusing COL_POS={COL_POS} COL_PHRED={COL_PHRED}"
        " -- fix at top of script if these look wrong",
        file=sys.stderr,
    )


def load_blocks(blocks_bed: Path) -> list[tuple[str, int, int, str]]:
    blocks = []
    with blocks_bed.open() as f:
        for line in f:
            chrom, start, end, name = line.rstrip("\n").split("\t")
            blocks.append((chrom, int(start), int(end), name))
    return blocks


def main():
    if len(sys.argv) < 3:
        print(
            f"usage: {sys.argv[0]} BLOCKS_BED AVI_PATH_OR_URL [OUT_TSV] [--inspect]",
            file=sys.stderr,
        )
        sys.exit(1)

    blocks_bed = Path(sys.argv[1])
    avi_path = sys.argv[2]
    do_inspect = "--inspect" in sys.argv
    positional = [a for a in sys.argv[3:] if a != "--inspect"]
    out_path = Path(positional[0]) if positional else Path("block_scores.tsv")

    blocks = load_blocks(blocks_bed)

    if do_inspect:
        chrom, start, end, _name = blocks[0]
        inspect(avi_path, chrom, start, end)
        return

    t0 = time.time()
    with out_path.open("w") as out:
        out.write(
            "block_id\tchrom\tstart\tend\tlength_kb\tn_snv\tavi_top1\tavi_dens\n"
        )
        for i, (chrom, start, end, name) in enumerate(blocks, 1):
            s = score_block(avi_path, chrom, start, end)
            out.write(
                f"{name}\t{chrom}\t{start}\t{end}\t{s['length_kb']:.3f}\t"
                f"{s['n_snv']}\t{s['avi_top1']}\t{s['avi_dens']}\n"
            )
            if i % 500 == 0:
                print(f"{i}/{len(blocks)} blocks, {time.time() - t0:.1f}s", file=sys.stderr)

    print(f"wrote {out_path}: {len(blocks)} blocks in {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
