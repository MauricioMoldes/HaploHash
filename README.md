# HaploHash

Haploblock-level functional annotation + privacy-preserving genomic hashes. 16–18 Sep 2026.

```
39,141 haploblocks (hg38)   +   AlphaGenome Atlas | gnomAD
            \                          /
             \___ aggregate to block __/   top-1% statistic
                         |
               annotated haploblocks
                 /                \
           haplograph          64-bit hash
   nodes / edges / islands   S|CHROM|HAPLO|CLUSTER|VAR
      lift-weighted          re-identification attack
```

![workflow](docs/image.png)

## Features per block

| Feature | Source | Note |
|---|---|---|
| `avi_top1` | Atlas AVI, Tabix | unsigned, permissive |
| `rnaseq_abs` | Atlas raw, API | signed — keep sum *and* abs |
| `gnocchi_mean` | gnomAD non-coding | positional, 1 kb |
| `loeuf_min` | gnomAD constraint | gene-joined; no gene = NA, not 0 |
| covariates | length, SNV density, coding fraction, genes |

Circular, excluded: PhastCons, Cactus, allele frequency (AVI inputs); recombination rate (defines blocks).

## Tracks

| Who | Track |
|---|---|
| Mauricio Moldes | block filtering, gene counts |
| Alejandra Caballero, Mina | annotation databases, pruning |
| Markus Marandi | AlphaGenome atlas |
| Robert Campbell, David Bonet | effect sizes (burden / SKAT / ACAT) |
| Aditya Kumar Karna | encoding, hashing, hash attack |

## Data

1000G — 2,548 individuals, 26 populations, `data.haploblocks.org`. hg38, SNVs only. UKB pending.
