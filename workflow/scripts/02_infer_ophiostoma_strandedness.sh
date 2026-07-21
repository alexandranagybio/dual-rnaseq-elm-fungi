#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

REPO_ROOT="${1:-$HOME/dual-rnaseq-elm-fungi}"
BED12="${2:-$REPO_ROOT/audit/strandedness/ophiostoma.bed12}"
BAM_DIR="${3:-$HOME/rnaseq/03_alignment/ophiostoma_alignments}"
OUT_DIR="$REPO_ROOT/audit/strandedness"

mkdir -p "$OUT_DIR"

if [[ ! -s "$BED12" ]]; then
    echo "ERROR: BED12 file is missing or empty: $BED12" >&2
    exit 1
fi

command -v infer_experiment.py >/dev/null 2>&1 || {
    echo "ERROR: infer_experiment.py is not available." >&2
    exit 1
}

bam_count=0

for bam in "$BAM_DIR"/*_ophiostoma.sorted.bam
do
    [[ -e "$bam" ]] || continue

    sample=$(basename "$bam" _ophiostoma.sorted.bam)
    output="$OUT_DIR/${sample}_infer_experiment.txt"

    printf 'Inferring strandedness: sample %s\n' "$sample"

    infer_experiment.py \
        -r "$BED12" \
        -i "$bam" \
        > "$output"

    bam_count=$((bam_count + 1))
done

if [[ "$bam_count" -ne 9 ]]; then
    echo "ERROR: expected 9 BAM files but processed $bam_count." >&2
    exit 1
fi

printf 'PASS: strandedness inferred for %d BAM files\n' "$bam_count"
