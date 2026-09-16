# HaploHash

Haploblock-level functional annotation + privacy-preserving genomic hashes.
Hackathon, 16–18 Sep 2026.

## Flow

```
data.haploblocks.org              alphagenome.google/downloads
39,141 haploblocks (hg38)         AVI (Tabix) | raw features (API)
            \                              /
             \___ aggregate to block level /
                          |
                annotated haploblocks
                  /                \
            haplograph          64-bit hash
    nodes / edges / islands   S|CHROM|HAPLO|CLUSTER|VAR
       lift-weighted          re-identification attack
```

Team whiteboard the flow came from:

![workflow](docs/image.png)

## What each stage is

**Haploblocks** — genome cut into recombination-defined regions, so each block is a unit that is actually inherited together.
- 39,141 blocks, GRCh38, chr1–22 + chrX
- per block, haplotypes clustered with MMseqs2; each individual gets a cluster ID
- block coordinates live in the filenames: `chr21_cluster_hashes_14215892-14284114.tsv`

**AlphaGenome Atlas** — precomputed effect predictions for every possible SNV in hg38.
- AVI score: Tabix, region-indexed, permissive license — `tabix avi.hg38.tsv.gz chr21:14215892-14284114`
- raw features (RNA-seq, DNase, ATAC, ChIP, CAGE): API only, non-commercial, rate-limited
- SNVs only for now; indels after final publication

**Aggregation** — the actual scientific claim. Lift variant-level scores to block level.
- candidate rules: mean, max, top-k percentile over the block
- paper's strongest result used top 1% of AVI (burden test), so top-k is the first thing to try
- swaps the paper's gene-anchored window for a recombination-defined block

**Haplograph** — how annotated blocks relate to each other.
- co-occurrence graph, already published: `nodes.csv.gz`, `edges.csv.gz`, `islands.csv.gz`
- edges weighted by lift; `islands.csv.gz` gives ready-made dense subgraphs to test
- annotated blocks join on as node attributes

**Hashing + attack** — privacy-preserving sharing and lookup.
- 64-bit packed field: `S(4) | CHROM(10) | HAPLO(20) | CLUSTER(20) | VAR(10)`
- not a one-way function — HASH is the binary cluster index, so this is an identifier, not a digest
- privacy comes from shipping cluster IDs instead of sequence
- so the attack is re-identification, not preimage: how few block IDs uniquely pin one person? rare and singleton clusters are the leak

## Tracks

| Who | Track |
|---|---|
| Mauricio Moldes | block filtering, genes per block |
| Alejandra Caballero, Mina | annotation databases, parameter pruning |
| Markus Marandi | AlphaGenome atlas: conservation, expression |
| Robert Campbell, David Bonet | variant → block effect sizes (burden / SKAT / ACAT) |
| Aditya Kumar Karna | encoding, hashing, hash attack |

Lead: Mauricio · Writer: Aditya

## Data

- 1000G — 2,548 individuals, 26 populations, from `data.haploblocks.org`
- hg38 on both sides, no liftover needed
- UKB pending access; published haploblocks are 1000G only, so UKB haplotypes must be assigned to existing clusters rather than reclustered

## Open

- aggregation rule not settled
- no UKB haplotype→cluster assignment step confirmed
- raw-feature API throughput unknown; scope expression work to one chromosome
