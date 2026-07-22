#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# 04_count_ophiostoma_genes.sh
#
# Purpose:
#   Generate a forward-stranded, gene-level featureCounts matrix for the nine
#   Ophiostoma novo-ulmi libraries listed in config/samples.tsv.
#
# Validated decisions:
#   - Library orientation: forward stranded
#   - featureCounts strandedness: -s 1
#   - Feature type: exon
#   - Grouping attribute: gene_id
#   - Expected annotation universe: 8,640 genes
#   - Paired-end fragments counted with -B and -C
#
# Usage:
#   bash workflow/scripts/04_count_ophiostoma_genes.sh
#
# Optional overrides:
#   THREADS=8 bash workflow/scripts/04_count_ophiostoma_genes.sh
# ---------------------------------------------------------------------------

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

readonly SAMPLE_TABLE="${SAMPLE_TABLE:-${REPO_ROOT}/config/samples.tsv}"
readonly GTF="${OPN_GTF:-${REPO_ROOT}/data/annotation/ophiostoma/annotation.gtf}"
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
    command -v "$1" >/dev/null 2>&1 ||
        die "Required command not found: $1"
}

require_file() {
    [[ -s "$1" ]] ||
        die "Missing or empty file: $1"
}

require_command featureCounts
require_command samtools
require_command awk
require_command sha256sum

require_file "${SAMPLE_TABLE}"
require_file "${GTF}"

mapfile -t BAM_RECORDS < <(
    awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "sample_id") sample_col = i
            if ($i == "ophiostoma_bam") bam_col = i
            if ($i == "include_ophiostoma") include_col = i
        }

        if (!sample_col || !bam_col || !include_col) {
            print "ERROR: required columns are missing from sample table" \
                > "/dev/stderr"
            exit 1
        }

        next
    }

    $include_col == "yes" {
        print $sample_col "\t" $bam_col
    }
    ' "${SAMPLE_TABLE}"
)

[[ "${#BAM_RECORDS[@]}" -eq 9 ]] ||
    die "Expected 9 Ophiostoma samples, found ${#BAM_RECORDS[@]} in ${SAMPLE_TABLE}"

SAMPLE_IDS=()
BAMS=()

for record in "${BAM_RECORDS[@]}"; do
    IFS=$'\t' read -r sample_id bam <<< "${record}"

    [[ -n "${sample_id}" ]] ||
        die "Encountered an empty sample ID in ${SAMPLE_TABLE}"

    [[ -n "${bam}" && "${bam}" != "NA" ]] ||
        die "Missing BAM path for Ophiostoma sample ${sample_id}"

    if [[ "${bam}" != /* ]]; then
        bam="${REPO_ROOT}/${bam}"
    fi

    require_file "${bam}"

    bai="${bam}.bai"

    if [[ ! -s "${bai}" ]]; then
        alternate_bai="${bam%.bam}.bai"

        if [[ -s "${alternate_bai}" ]]; then
            bai="${alternate_bai}"
        else
            die "Missing BAM index for ${bam}"
        fi
    fi

    samtools quickcheck -v "${bam}" ||
        die "samtools quickcheck failed for ${bam}"

    SAMPLE_IDS+=("${sample_id}")
    BAMS+=("${bam}")
done

readonly RAW_COUNTS="${OUT_DIR}/ophiostoma_gene_featurecounts.txt"
readonly SUMMARY="${RAW_COUNTS}.summary"
readonly RUN_INFO="${OUT_DIR}/04_count_ophiostoma_genes.run_info.txt"
readonly CHECKSUMS="${OUT_DIR}/04_count_ophiostoma_genes.sha256"

log "Repository root: ${REPO_ROOT}"
log "Sample table: ${SAMPLE_TABLE}"
log "Validated GTF: ${GTF}"
log "Output directory: ${OUT_DIR}"
log "Expected genes: ${EXPECTED_GENES}"
log "Ophiostoma samples: ${SAMPLE_IDS[*]}"
log "Strandedness: forward (-s 1)"
log "Starting featureCounts."

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
    awk '
    BEGIN {
        n = 0
    }

    !/^#/ && $1 != "Geneid" {
        n++
    }

    END {
        print n
    }
    ' "${RAW_COUNTS}"
)"

[[ "${observed_genes}" =~ ^[0-9]+$ ]] ||
    die "Could not determine the number of gene rows in ${RAW_COUNTS}"

if [[ "${observed_genes}" -ne "${EXPECTED_GENES}" ]]; then
    die "Expected ${EXPECTED_GENES} gene rows, but featureCounts produced ${observed_genes}. Do not continue to DESeq2."
fi

{
    printf 'run_timestamp\t%s\n' "$(timestamp)"
    printf 'repository_root\t%s\n' "${REPO_ROOT}"
    printf 'sample_table\t%s\n' "${SAMPLE_TABLE}"
    printf 'annotation_gtf\t%s\n' "${GTF}"
    printf 'expected_gene_rows\t%s\n' "${EXPECTED_GENES}"
    printf 'observed_gene_rows\t%s\n' "${observed_genes}"
    printf 'feature_type\texon\n'
    printf 'grouping_attribute\tgene_id\n'
    printf 'strandedness\t1\n'
    printf 'paired_end\ttrue\n'
    printf 'require_both_ends_mapped\ttrue\n'
    printf 'exclude_chimeric_fragments\ttrue\n'
    printf 'threads\t%s\n' "${THREADS}"

    printf '\n[samples]\n'
    for i in "${!SAMPLE_IDS[@]}"; do
        printf '%s\t%s\n' "${SAMPLE_IDS[$i]}" "${BAMS[$i]}"
    done

    printf '\n[featureCounts version]\n'
    featureCounts -v 2>&1 || true

    printf '\n[samtools version]\n'
    samtools --version | head -n 2
} > "${RUN_INFO}"

sha256sum \
    "${SAMPLE_TABLE}" \
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
