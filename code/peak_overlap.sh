#!/bin/bash -l
#
# Identify unique and shared H3K27ac peaks across the four MCF7 conditions
# using bedtools, and report the peak counts used to build the Venn diagram
# (Figure 1A). Run interactively (the intersections are lightweight):
#
#   srun --partition=amdsmall --cpus-per-task=1 --mem=30G --time=2:00:00 --pty bash -i
#   conda activate chipseq      # or: module load bedtools
#   bash peak_overlap.sh
#
# Inputs : the MACS2 *_summits.bed file for each condition (see results/).
# Outputs: pairwise / triple / all-four shared BEDs, per-condition unique BEDs,
#          and a printed peak-count summary.

set -euo pipefail

# ---- Input summit files (one per condition) ----
CONTROL="Control_shRNA_summits.bed"
GATA3="GATA3_shRNA_summits.bed"
GATA3_JUN="GATA3_shRNA_JUN_OE_summits.bed"
JUN="JUN_OE_DOX_summits.bed"

OUTDIR="peak_overlap"
mkdir -p "$OUTDIR"

# ---- Pairwise intersections ----
bedtools intersect -a "$CONTROL"   -b "$GATA3"      > "$OUTDIR/shared_Control_GATA3.bed"
bedtools intersect -a "$CONTROL"   -b "$GATA3_JUN"  > "$OUTDIR/shared_Control_GATA3JUN.bed"
bedtools intersect -a "$CONTROL"   -b "$JUN"        > "$OUTDIR/shared_Control_JUN.bed"
bedtools intersect -a "$GATA3"     -b "$GATA3_JUN"  > "$OUTDIR/shared_GATA3_GATA3JUN.bed"
bedtools intersect -a "$GATA3"     -b "$JUN"        > "$OUTDIR/shared_GATA3_JUN.bed"
bedtools intersect -a "$GATA3_JUN" -b "$JUN"        > "$OUTDIR/shared_GATA3JUN_JUN.bed"

# ---- Triple intersections ----
bedtools intersect -a "$CONTROL" -b "$GATA3" "$JUN"       > "$OUTDIR/shared_Control_GATA3_JUN.bed"
bedtools intersect -a "$CONTROL" -b "$GATA3" "$GATA3_JUN" > "$OUTDIR/shared_Control_GATA3_GATA3JUN.bed"
bedtools intersect -a "$CONTROL" -b "$JUN"   "$GATA3_JUN" > "$OUTDIR/shared_Control_JUN_GATA3JUN.bed"
bedtools intersect -a "$GATA3"   -b "$JUN"   "$GATA3_JUN" > "$OUTDIR/shared_GATA3_JUN_GATA3JUN.bed"

# ---- Shared by all four ----
bedtools intersect -a "$CONTROL" -b "$GATA3" "$GATA3_JUN" "$JUN" > "$OUTDIR/shared_all4.bed"

# ---- Unique to each condition (-v = no overlap with any of the others) ----
bedtools intersect -a "$CONTROL"   -b "$GATA3" "$GATA3_JUN" "$JUN" -v > "$OUTDIR/unique_Control.bed"
bedtools intersect -a "$GATA3"     -b "$CONTROL" "$GATA3_JUN" "$JUN" -v > "$OUTDIR/unique_GATA3.bed"
bedtools intersect -a "$GATA3_JUN" -b "$CONTROL" "$GATA3" "$JUN" -v     > "$OUTDIR/unique_GATA3JUN.bed"
bedtools intersect -a "$JUN"       -b "$CONTROL" "$GATA3" "$GATA3_JUN" -v > "$OUTDIR/unique_JUN.bed"

# ---- Peak-count summary (feeds the Venn diagram) ----
echo "=== Peak counts ==="
wc -l "$OUTDIR"/*.bed
