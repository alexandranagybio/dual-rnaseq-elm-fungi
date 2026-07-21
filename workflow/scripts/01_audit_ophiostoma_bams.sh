#!/usr/bin/env bash
set -euo pipefail

BAM_DIR="${1:-/home/alex/rnaseq/03_alignment/ophiostoma_alignments}"
OUT_DIR="${2:-audit/bam_qc}"

mkdir -p "$OUT_DIR"

printf "sample\tbam\tquickcheck_status\treference_header_sha256\n" \
  > "$OUT_DIR/ophiostoma_bam_audit.tsv"

for bam in "$BAM_DIR"/*_ophiostoma.sorted.bam
do
    sample=$(basename "$bam" | cut -d_ -f1)

    if samtools quickcheck -v "$bam" \
        > "$OUT_DIR/${sample}_quickcheck.txt" 2>&1
    then
        quickcheck_status="pass"
    else
        quickcheck_status="fail"
    fi

    header_hash=$(
        samtools view -H "$bam" |
        awk '$1=="@SQ"{print $2"\t"$3}' |
        sha256sum |
        cut -d' ' -f1
    )

    samtools view -H "$bam" \
        > "$OUT_DIR/${sample}_header.sam"

    samtools flagstat -@ 4 "$bam" \
        > "$OUT_DIR/${sample}_flagstat.txt"

    samtools stats -@ 4 "$bam" \
        > "$OUT_DIR/${sample}_samtools_stats.txt"

    printf "%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$bam" \
        "$quickcheck_status" \
        "$header_hash" \
        >> "$OUT_DIR/ophiostoma_bam_audit.tsv"
done
