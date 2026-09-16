#!/usr/bin/env python3
"""Scrape haploblock interval filenames from data.haploblocks.org and write blocks.bed (hg38).

Usage: python3 build_blocks_bed.py [out.bed]
"""
import re
import sys
import urllib.request
from pathlib import Path

ROOT = "https://data.haploblocks.org/haploblock_hashes/1000G"
CHROMS = [str(i) for i in range(1, 23)] + ["X"]
FILE_RE = re.compile(r'href="(chr\w+_cluster_hashes_(\d+)-(\d+)\.tsv)"')
EXPECTED_TOTAL = 39141


def fetch(url: str) -> str:
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode("utf-8", errors="replace")


def scrape_chrom(chrom: str) -> list[tuple[str, int, int, str]]:
    html = fetch(f"{ROOT}/chr{chrom}/")
    rows = []
    for _fname, start, end in FILE_RE.findall(html):
        start, end = int(start), int(end)
        name = f"chr{chrom}_{start}_{end}"
        rows.append((f"chr{chrom}", start, end, name))
    return rows


def check_overlaps(rows: list[tuple[str, int, int, str]]) -> list[str]:
    problems = []
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for chrom, start, end, name in rows:
        by_chrom.setdefault(chrom, []).append((start, end, name))
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        for (s1, e1, n1), (s2, e2, n2) in zip(intervals, intervals[1:]):
            if s2 < e1:
                problems.append(f"{chrom}: {n1} ({s1}-{e1}) overlaps {n2} ({s2}-{e2})")
    return problems


def main():
    out_path = Path(sys.argv[1] if len(sys.argv) > 1 else "blocks.bed")

    all_rows = []
    for chrom in CHROMS:
        rows = scrape_chrom(chrom)
        print(f"chr{chrom}: {len(rows)} blocks", file=sys.stderr)
        all_rows.extend(rows)

    chrom_order = {f"chr{c}": i for i, c in enumerate(CHROMS)}
    all_rows.sort(key=lambda r: (chrom_order[r[0]], r[1]))

    with out_path.open("w") as f:
        for chrom, start, end, name in all_rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(f"\ntotal blocks: {len(all_rows)} (expected ~{EXPECTED_TOTAL})", file=sys.stderr)
    if len(all_rows) != EXPECTED_TOTAL:
        print("WARNING: count mismatch vs expected", file=sys.stderr)

    problems = check_overlaps(all_rows)
    if problems:
        print(f"WARNING: {len(problems)} overlapping block pairs:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  {p}", file=sys.stderr)
    else:
        print("no overlaps", file=sys.stderr)

    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
