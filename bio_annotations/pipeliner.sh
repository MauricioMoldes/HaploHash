#!/usr/bin/env bash
############
##Create env
conda create -n haploblock_mvp -c conda-forge -c bioconda \
    minimap2 \
    samtools \
    bedtools 

##Genome download
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38.knownGene.gtf.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
############

##Pipeline
set -euo pipefail

#0. Input files 
##reference genome (fasta/gtf)
REFERENCE="/dcai/projects/iu_0142/team11/ref_genome/hg38.fna"
GTF="/dcai/projects/iu_0142/team11/ref_genome/hg38.knownGene.gtf"

##haploblocks of interest, maybe a dir makes more sense
HAPLOBLOCKS_DIR="/dcai/projects/iu_0142/team11/data/haploblock_sequences/1000G/chr22/clusters"

OUTDIR="ch22_genes"
NAME='ch22'
mkdir -p "${OUTDIR}"

#1. Check present tools
##a. minimap2
command -v minimap2 >/dev/null 2>&1 || {
    echo "ERROR: minimap2 not found"
    exit 1
}

##b. samtools
command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found"
    exit 1
}

##c. bedtools
command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools not found"
    exit 1
}

#2. Index reference genome
echo ">>> Indexing reference genome..."

if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx "${REFERENCE}"
fi

if [ ! -f "${REFERENCE}.mmi" ]; then
    minimap2 -d "${REFERENCE}.mmi" "${REFERENCE}"
fi

#3. Extract genes and CDS from GTF 
echo ">>> Preparing GTF annotations..."

GENES_GTF="${OUTDIR}/genes.gtf"
CDS_GTF="${OUTDIR}/cds.gtf"

if [ ! -f "$GENES_GTF" ]; then
    awk '$3 == "gene"' "${GTF}" > "$GENES_GTF"
fi

if [ ! -f "$CDS_GTF" ]; then
    awk '$3 == "CDS"' "${GTF}" > "$CDS_GTF"
fi

#4. Process each haploblock
for HAPLOBLOCKS in "$HAPLOBLOCKS_DIR"/*_all_seqs.fasta; do
    [[ -e "$HAPLOBLOCKS" ]] || continue

    #Name change
    NAME=$(basename "$HAPLOBLOCKS" _all_seqs.fasta)

    if [[ -f "${OUTDIR}/${NAME}_haploblock_cds.tsv" ]]; then
        echo ">>> ${NAME} already exists, skipping..."
        continue
    fi
    
    #4.1 Prepare FASTA for mapping
    echo ">>> Checking FASTA..."
    CLEAN_FASTA="${OUTDIR}/${NAME}.clean.fasta"

    awk '
    /^>/ {
        if (seen && len > 0) {
            print header
            print seq
        }
        header=$0
        seq=""
        len=0
        seen=1
        next
    }
    {
        seq=seq $0
        len+=length($0)
    }
    END {
        if (seen && len > 0) {
            print header
            print seq
        }
    }
    ' "$HAPLOBLOCKS" > "$CLEAN_FASTA"

    # Skip files with no valid sequences
    if [ ! -s "$CLEAN_FASTA" ] || ! grep -q '^>' "$CLEAN_FASTA"; then
        echo "WARNING: No valid sequences found in $NAME. Skipping."
        rm -f "$CLEAN_FASTA"
        continue
    fi

    #4.2 Map haploblocks to GRCh38
    echo ">>> Mapping haploblocks..."

    minimap2 \
        -ax sr \
        "${REFERENCE}.mmi" \
        "${CLEAN_FASTA}" \
        > "${OUTDIR}/${NAME}.sam"

    rm "$CLEAN_FASTA"

    #4.3 SAM == sorted BAM
    echo ">>> Converting SAM to sorted BAM..."

    samtools view \
        -b \
        "${OUTDIR}/${NAME}.sam" |
        samtools sort \
        -o "${OUTDIR}/${NAME}.bam"

    samtools index "${OUTDIR}/${NAME}.bam"

    #4.4 Keep primary alignments
    echo ">>> Extracting primary alignments..."

    samtools view \
        -b \
        -F 256 \
        "${OUTDIR}/${NAME}.bam" \
        > "${OUTDIR}/${NAME}.primary.bam"

    samtools index "${OUTDIR}/${NAME}.primary.bam"

    #4.5 BAM == BED
    echo ">>> Creating BED..."

    bedtools bamtobed \
        -i "${OUTDIR}/${NAME}.primary.bam" \
        > "${OUTDIR}/${NAME}.bed"

    #4.6 Haploblock == genes
    echo ">>> Annotating haploblocks with genes..."

    bedtools intersect \
        -a "${OUTDIR}/${NAME}.bed" \
        -b "${GENES_GTF}" \
        -wa \
        -wb \
        > "${OUTDIR}/${NAME}_haploblock_genes.tsv"

    #4.7 Haploblock == CDS
    echo ">>> Annotating haploblocks with CDS..."

    bedtools intersect \
        -a "${OUTDIR}/${NAME}.bed" \
        -b "${CDS_GTF}" \
        -wao \
        > "${OUTDIR}/${NAME}_haploblock_cds.tsv"

    #4.8 Mapping statistics
    echo ">>> Calculating mapping statistics..."

    samtools flagstat \
        "${OUTDIR}/${NAME}.primary.bam" \
        > "${OUTDIR}/${NAME}_mapping_stats.txt"

    echo ""
    echo ">>> Finished: $NAME"
    echo ""
    echo "Outputs:"
    echo "  BED:    ${OUTDIR}/${NAME}.bed"
    echo "  Genes:  ${OUTDIR}/${NAME}_haploblock_genes.tsv"
    echo "  CDS:    ${OUTDIR}/${NAME}_haploblock_cds.tsv"
    echo "  Stats:  ${OUTDIR}/${NAME}_mapping_stats.txt"

done

echo ""
echo "===="
echo "DONE"
echo "===="
