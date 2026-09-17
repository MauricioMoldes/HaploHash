# How thousands of SNVs become one number per block

**In plain words:**

```
top-1% fraction  =  dangerous mutations found in this block
                     ---------------------------------------
                     all mutations checked in this block
```

"Dangerous" = scored PHRED >= 20 = worse than 99% of all possible mutations genome-wide.

A haploblock is a region, e.g. `chr1:1,583,825-1,958,353`.

The AVI file has one row per possible mutation in the genome (every position, every possible letter swap), each with an impact score (PHRED).

**Step 1 — collect:** grab every AVI row whose position falls inside that block's start-end range. Everything outside the range is ignored.

**Step 2 — count:** for the rows you collected —

| SNV position | inside block? | PHRED score | counts as "top hit"? (>=20) |
|---|---|---|---|
| 1,583,900 | yes | 25.1 | yes |
| 1,584,010 | yes | 8.3 | no |
| 1,584,010 | yes | 41.7 | yes |
| 2,000,500 | no (outside block) | - | ignored |
| 1,600,200 | yes | 12.0 | no |

In this toy example: 4 rows landed inside the block, 2 of them scored >=20.

**Step 3 — divide:**

```
avi_top1_fraction = (rows inside block scoring >= 20) / (all scored rows inside block)
                   = 2 / 4 = 0.50  (50%)
```

Do this once per block (39,074 times), and you get one `avi_top1_fraction` number per block — that's the y-axis value plotted in the charts.

Source: [`annotate_atlas/score_blocks.py`](../annotate_atlas/score_blocks.py), `aggregate_rows()` + `BlockSummary`.
