#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# 04_count_ophiostoma_genes.sh
#
# Purpose:
#   Generate a stranded, gene-level featureCounts matrix for the nine
#   Ophiostoma novo-ulmi libraries.
#
# Validated biological/technical decisions:
#   - Library orientation: forward stranded
#   - featureCounts strandedness: -s 1
#   - Expected annotation universe: 8,640 genes
#   - Samples: 149-157
#
# Repository/data layout:
#   This script is intended to live in:
#     ~/dual-rnaseq-elm-fungi/workflow/scripts/
#
#   The BAM files remain in:
#     ~/rnaseq/03_alignment/ophiostoma_alignments/
#
# Usage:
#   chmod +x workflow/scripts/04_count_ophiostoma_genes.sh
#   bash workflow/scripts/04_count_ophiostoma_genes.sh
#
# Optional overrides:
#   RNASEQ_ROOT=/different/path \
#   OPN_GTF=/path/to/validated.annotation.gtf \
#   THREADS=8 \
#   bash workflow/scripts/04_count_ophiostoma_genes.sh
# ---------------------------------------------------------------------------

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

readonly RNASEQ_ROOT="${RNASEQ_ROOT:-${HOME}/rnaseq}"
readonly BAM_DIR="${BAM_DIR:-${RNASEQ_ROOT}/03_alignment/ophiostoma_alignments}"
readonly OUT_DIR="${OUT_DIR:-${REPO_ROOT}/results/ophiostoma/gene_counts}"
readonly LOG_DIR="${LOG_DIR:-${REPO_ROOT}/logs}"
readonly THREADS="${THREADS:-4}"
readonly EXPECTED_GENES="${EXPECTED_GENES:-8640}"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

readonly LOG_FILE="${LOG_DIR}/04_count_ophiostoma_genes.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

die() {
    printf '[%s] ERROR: %s\n' "$(timestamp)" "$*" >&2
    exit 1
}

log() {
    printf '[%s] %s\n' "$(timestamp)" "$*"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

find_validated_gtf() {
    local candidates=()

    if [[ -n "${OPN_GTF:-}" ]]; then
        printf '%s\n' "${OPN_GTF}"
        return 0
    fi

    # Prefer audit-repository annotation products.
    while IFS= read -r -d '' f; do
        candidates+=("$f")
    done < <(
        find "${REPO_ROOT}" -type f \
            \( -iname '*ophiostoma*.gtf' -o -iname '*ophnu1*.gtf' -o -iname '*annotation*.gtf' \) \
            -print0 2>/dev/null
    )

    # Fall back to the RNA-seq project only if exactly one plausible GTF exists.
    while IFS= read -r -d '' f; do
        candidates+=("$f")
    done < <(
        find "${RNASEQ_ROOT}/02_reference" -type f \
            \( -iname '*ophiostoma*.gtf' -o -iname '*ophnu1*.gtf' -o -iname '*annotation*.gtf' \) \
            -print0 2>/dev/null
    )

    if (( ${#candidates[@]} == 1 )); then
        printf '%s\n' "${candidates[0]}"
        return 0
    fi

    if (( ${#candidates[@]} == 0 )); then
        return 1
    fi

    printf 'Multiple candidate GTF files were found:\n' >&2
    printf '  %s\n' "${candidates[@]}" >&2
    printf 'Set OPN_GTF explicitly to the validated 8,640-gene GTF.\n' >&2
    return 1
}

require_command featureCounts
require_command samtools
require_command awk
require_command python3
require_command sha256sum

GTF="$(find_validated_gtf)" || die \
    "Could not uniquely identify the validated GTF. Run with OPN_GTF=/absolute/path/to/validated.gtf"

require_file "${GTF}"

SAMPLE_IDS=(149 150 151 152 153 154 155 156 157)
BAMS=()

for sample_id in "${SAMPLE_IDS[@]}"; do
    bam="${BAM_DIR}/${sample_id}_ophiostoma.sorted.bam"
    bai="${bam}.bai"
    require_file "${bam}"
    require_file "${bai}"

    samtools quickcheck -v "${bam}" \
        || die "samtools quickcheck failed for ${bam}"

    BAMS+=("${bam}")
done

readonly RAW_COUNTS="${OUT_DIR}/ophiostoma_gene_featurecounts.txt"
readonly SUMMARY="${RAW_COUNTS}.summary"
readonly RUN_INFO="${OUT_DIR}/04_count_ophiostoma_genes.run_info.txt"
readonly CHECKSUMS="${OUT_DIR}/04_count_ophiostoma_genes.sha256"

log "Repository root: ${REPO_ROOT}"
log "RNA-seq data root: ${RNASEQ_ROOT}"
log "BAM directory: ${BAM_DIR}"
log "Validated GTF: ${GTF}"
log "Output directory: ${OUT_DIR}"
log "Expected genes: ${EXPECTED_GENES}"
log "Strandedness: forward (-s 1)"
log "Starting featureCounts."

# The validated GTF produced during the audit must contain gene_id attributes.
# We count exon features and aggregate them by gene_id, yielding one row per gene.
featureCounts \
    -T "${THREADS}" \
    -p \
    --countReadPairs \
    -B \
    -C \
    -s 1 \
    -t exon \
    -g gene_id \
    -a "${GTF}" \
    -o "${RAW_COUNTS}" \
    "${BAMS[@]}"

require_file "${RAW_COUNTS}"
require_file "${SUMMARY}"

observed_genes="$(
    awk 'BEGIN{n=0} !/^#/ && $1!="Geneid" {n++} END{print n}' "${RAW_COUNTS}"
)"

[[ "${observed_genes}" =~ ^[0-9]+$ ]] \
    || die "Could not determine the number of rows in ${RAW_COUNTS}"

if [[ "${observed_genes}" -ne "${EXPECTED_GENES}" ]]; then
    die "Expected ${EXPECTED_GENES} gene rows, but featureCounts produced ${observed_genes}. Do not continue to DESeq2."
fi

{
    printf 'run_timestamp\t%s\n' "$(timestamp)"
    printf 'repository_root\t%s\n' "${REPO_ROOT}"
    printf 'rnaseq_root\t%s\n' "${RNASEQ_ROOT}"
    printf 'annotation_gtf\t%s\n' "${GTF}"
    printf 'bam_directory\t%s\n' "${BAM_DIR}"
    printf 'expected_gene_rows\t%s\n' "${EXPECTED_GENES}"
    printf 'observed_gene_rows\t%s\n' "${observed_genes}"
    printf 'feature_type\texon\n'
    printf 'grouping_attribute\tgene_id\n'
    printf 'strandedness\t1\n'
    printf 'paired_end\ttrue\n'
    printf 'require_both_ends_mapped\ttrue\n'
    printf 'exclude_chimeric_fragments\ttrue\n'
    printf 'threads\t%s\n' "${THREADS}"
    printf '\n[featureCounts version]\n'
    featureCounts -v 2>&1 || true
    printf '\n[samtools version]\n'
    samtools --version | head -n 2
    printf '\n[python version]\n'
    python3 --version
    printf '\n[input BAM files]\n'
    printf '%s\n' "${BAMS[@]}"
} > "${RUN_INFO}"

sha256sum \
    "${GTF}" \
    "${BAMS[@]}" \
    "${RAW_COUNTS}" \
    "${SUMMARY}" \
    "${RUN_INFO}" \
    > "${CHECKSUMS}"

log "SUCCESS: generated ${observed_genes}-gene featureCounts table."
log "Counts: ${RAW_COUNTS}"
log "Assignment summary: ${SUMMARY}"
log "Run information: ${RUN_INFO}"
log "Checksums: ${CHECKSUMS}"
