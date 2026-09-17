#!/bin/bash
#SBATCH --job-name=avi_score_blocks
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/avi_score_blocks_%j.log

set -euo pipefail

REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
BIN=$REPO/annotate_atlas/bin
SRC_TSV=/dcai/projects/iu_0142/team11/data/avi/alphagenome_variant_impact_score_snvs.tsv.gz
NWORKERS=32
LOCAL=/tmp/avi_${SLURM_JOB_ID}
OUT=$REPO/block_scores.tsv

mkdir -p "$LOCAL"
trap 'rm -rf "$LOCAL"' EXIT

echo "[$(date)] staging AVI file to local disk on $(hostname)"
cp "$SRC_TSV" "$SRC_TSV.tbi" "$LOCAL/"
LOCAL_TSV="$LOCAL/$(basename "$SRC_TSV")"

echo "[$(date)] splitting blocks.bed into $NWORKERS chunks"
split -n l/$NWORKERS -d -a 2 "$REPO/annotate_atlas/blocks.bed" "$LOCAL/chunk_"

echo "[$(date)] launching $NWORKERS parallel workers"
for f in "$LOCAL"/chunk_*; do
  idx=$(basename "$f")
  python3 "$REPO/annotate_atlas/score_blocks.py" "$f" "$LOCAL_TSV" \
    "$LOCAL/out_${idx}.tsv" --tabix "$BIN/tabix" \
    > "$LOCAL/log_${idx}.txt" 2>&1 &
done
wait

echo "[$(date)] all workers finished, checking logs for errors"
for f in "$LOCAL"/log_*.txt; do
  if grep -qiE 'error|traceback' "$f"; then
    echo "FAILED: $f"
    cat "$f"
    exit 1
  fi
done

echo "[$(date)] merging chunk outputs"
head -n1 "$LOCAL"/out_chunk_00.tsv > "$OUT"
for f in "$LOCAL"/out_chunk_*.tsv; do
  tail -n +2 "$f" >> "$OUT"
done

echo "[$(date)] done. rows: $(wc -l < "$OUT")"
