#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C

REPO_ROOT="${1:-$HOME/dual-rnaseq-elm-fungi}"
REPORT_DIR="$REPO_ROOT/audit/strandedness"
OUTPUT="$REPORT_DIR/ophiostoma_strandedness_summary.tsv"

printf '%s\t%s\t%s\t%s\t%s\n' \
    sample_id \
    failed_fraction \
    forward_fraction \
    reverse_fraction \
    inferred_featurecounts_setting \
    > "$OUTPUT"

report_count=0

for report in "$REPORT_DIR"/[0-9]*_infer_experiment.txt
do
    [[ -e "$report" ]] || continue

    sample=$(basename "$report" _infer_experiment.txt)

    failed=$(
        awk -F': ' \
            '/Fraction of reads failed to determine:/ {print $2}' \
            "$report"
    )

    forward=$(
        awk -F': ' \
            '/1\+\+,1--,2\+-,2-\+/{print $2}' \
            "$report"
    )

    reverse=$(
        awk -F': ' \
            '/1\+-,1-\+,2\+\+,2--/{print $2}' \
            "$report"
    )

    if [[ -z "$failed" || -z "$forward" || -z "$reverse" ]]; then
        echo "ERROR: could not parse $report" >&2
        exit 1
    fi

    setting=$(
        awk -v f="$forward" -v r="$reverse" '
        BEGIN {
            if (f >= 0.80 && f > r)
                print "-s 1"
            else if (r >= 0.80 && r > f)
                print "-s 2"
            else
                print "-s 0 or ambiguous"
        }'
    )

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$sample" \
        "$failed" \
        "$forward" \
        "$reverse" \
        "$setting" \
        >> "$OUTPUT"

    report_count=$((report_count + 1))
done

if [[ "$report_count" -ne 9 ]]; then
    echo "ERROR: expected 9 reports but parsed $report_count." >&2
    exit 1
fi

echo "Created $OUTPUT"
