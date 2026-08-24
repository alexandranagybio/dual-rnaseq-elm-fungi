#!/usr/bin/env python3
"""
35_build_extracellular_summary_table.py

Build the canonical publication summary of extracellular transcriptional
responses in Fusarium cf. salinense and Ophiostoma novo-ulmi.

No new statistical test is performed.

Definitions
-----------
Significant DE gene:
    adjusted P value < 0.05

Secreted gene:
    confident SignalP-positive gene

High-confidence CAZyme:
    at least one gene-associated protein supported by >= 2 dbCAN tools

Secreted CAZyme:
    intersection of significant, confident SignalP-positive, and
    high-confidence dbCAN-positive genes

Contrasts
---------
Fusarium:
    interaction versus self/control

Ophiostoma:
    interaction versus self

Outputs
-------
results/publication/extracellular_response/
    extracellular_summary.tsv
    extracellular_summary_long.tsv
    extracellular_overlap_long.tsv
    extracellular_audit.tsv
    run_info.tsv
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path.cwd()

FUSARIUM_INPUT = (
    ROOT
    / "results/fusarium/publication_annotation/tables/"
      "fusarium_publication_annotation_full.tsv"
)

OPHIOSTOMA_INPUT = (
    ROOT
    / "results/ophiostoma/cazyme_analysis/tables/"
      "interaction_vs_onu_all_genes_with_dbcan.tsv"
)

OUTDIR = (
    ROOT
    / "results/publication/extracellular_response"
)

SUMMARY_OUT = OUTDIR / "extracellular_summary.tsv"
SUMMARY_LONG_OUT = OUTDIR / "extracellular_summary_long.tsv"
OVERLAP_LONG_OUT = OUTDIR / "extracellular_overlap_long.tsv"
AUDIT_OUT = OUTDIR / "extracellular_audit.tsv"
RUN_INFO_OUT = OUTDIR / "run_info.tsv"

ALPHA = 0.05

EXPECTED_SIGNIFICANT = {
    "Fusarium": 2973,
    "Ophiostoma": 4938,
}

EXPECTED_SIGNIFICANT_CAZYMES = {
    "Fusarium": 142,
    "Ophiostoma": 206,
}

EMPTY_VALUES = {
    "",
    "NA",
    "NaN",
    "nan",
    "NULL",
    "null",
    "None",
    "none",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def read_tsv(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    require_file(path)

    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            fail(f"no header found in {path}")

        return list(reader.fieldnames), list(reader)


def write_tsv(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, object]],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def parse_float(
    value: str,
    field: str,
    gene_id: str,
) -> float:
    if value.strip() in EMPTY_VALUES:
        return float("nan")

    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid {field} for {gene_id}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return math.isnan(value)


def parse_bool(
    value: str,
    field: str,
    gene_id: str,
) -> bool:
    normalized = value.strip().upper()

    if normalized in {"TRUE", "T", "1", "YES", "Y"}:
        return True

    if normalized in {"FALSE", "F", "0", "NO", "N", ""}:
        return False

    fail(
        f"invalid Boolean in {field} for {gene_id}: {value}"
    )

    return False


def require_columns(
    path: Path,
    observed: list[str],
    required: set[str],
) -> None:
    missing = required - set(observed)

    if missing:
        fail(
            f"{path} is missing columns: {sorted(missing)}"
        )


def process_fusarium() -> tuple[dict[str, object], set[str]]:
    fields, rows = read_tsv(FUSARIUM_INPUT)

    required = {
        "gene_id",
        "in_deseq2_dataset",
        "raw_padj",
        "signalp_confident_positive",
        "dbcan_high_confidence",
    }

    require_columns(
        FUSARIUM_INPUT,
        fields,
        required,
    )

    seen: set[str] = set()

    significant: set[str] = set()
    significant_secreted: set[str] = set()
    significant_cazymes: set[str] = set()
    significant_secreted_cazymes: set[str] = set()

    for row in rows:
        gene_id = row["gene_id"]

        if gene_id in seen:
            fail(
                f"duplicated Fusarium gene ID: {gene_id}"
            )

        seen.add(gene_id)

        in_dataset = parse_bool(
            row["in_deseq2_dataset"],
            "in_deseq2_dataset",
            gene_id,
        )

        if not in_dataset:
            continue

        padj = parse_float(
            row["raw_padj"],
            "raw_padj",
            gene_id,
        )

        is_significant = (
            not is_nan(padj)
            and padj < ALPHA
        )

        if not is_significant:
            continue

        is_secreted = parse_bool(
            row["signalp_confident_positive"],
            "signalp_confident_positive",
            gene_id,
        )

        is_cazyme = parse_bool(
            row["dbcan_high_confidence"],
            "dbcan_high_confidence",
            gene_id,
        )

        significant.add(gene_id)

        if is_secreted:
            significant_secreted.add(gene_id)

        if is_cazyme:
            significant_cazymes.add(gene_id)

        if is_secreted and is_cazyme:
            significant_secreted_cazymes.add(gene_id)

    summary = summarize_species(
        species="Fusarium",
        significant=significant,
        significant_secreted=significant_secreted,
        significant_cazymes=significant_cazymes,
        significant_secreted_cazymes=(
            significant_secreted_cazymes
        ),
        input_file=FUSARIUM_INPUT,
    )

    return summary, seen


def process_ophiostoma() -> tuple[dict[str, object], set[str]]:
    fields, rows = read_tsv(OPHIOSTOMA_INPUT)

    required = {
        "gene_id",
        "contrast",
        "padj",
        "signalp_is_sp",
        "dbcan_high_confidence",
    }

    require_columns(
        OPHIOSTOMA_INPUT,
        fields,
        required,
    )

    seen: set[str] = set()

    significant: set[str] = set()
    significant_secreted: set[str] = set()
    significant_cazymes: set[str] = set()
    significant_secreted_cazymes: set[str] = set()

    expected_contrast = "interaction_vs_onu"

    for row in rows:
        gene_id = row["gene_id"]

        if gene_id in seen:
            fail(
                f"duplicated Ophiostoma gene ID: {gene_id}"
            )

        seen.add(gene_id)

        if row["contrast"] != expected_contrast:
            fail(
                f"unexpected Ophiostoma contrast for {gene_id}: "
                f"{row['contrast']}"
            )

        padj = parse_float(
            row["padj"],
            "padj",
            gene_id,
        )

        is_significant = (
            not is_nan(padj)
            and padj < ALPHA
        )

        if not is_significant:
            continue

        is_secreted = parse_bool(
            row["signalp_is_sp"],
            "signalp_is_sp",
            gene_id,
        )

        is_cazyme = parse_bool(
            row["dbcan_high_confidence"],
            "dbcan_high_confidence",
            gene_id,
        )

        significant.add(gene_id)

        if is_secreted:
            significant_secreted.add(gene_id)

        if is_cazyme:
            significant_cazymes.add(gene_id)

        if is_secreted and is_cazyme:
            significant_secreted_cazymes.add(gene_id)

    summary = summarize_species(
        species="Ophiostoma",
        significant=significant,
        significant_secreted=significant_secreted,
        significant_cazymes=significant_cazymes,
        significant_secreted_cazymes=(
            significant_secreted_cazymes
        ),
        input_file=OPHIOSTOMA_INPUT,
    )

    return summary, seen


def summarize_species(
    species: str,
    significant: set[str],
    significant_secreted: set[str],
    significant_cazymes: set[str],
    significant_secreted_cazymes: set[str],
    input_file: Path,
) -> dict[str, object]:
    nonsecreted_cazymes = (
        significant_cazymes
        - significant_secreted_cazymes
    )

    secreted_non_cazymes = (
        significant_secreted
        - significant_secreted_cazymes
    )

    significant_nonsecreted_non_cazymes = (
        significant
        - significant_secreted
        - significant_cazymes
    )

    significant_count = len(significant)
    secreted_count = len(significant_secreted)
    cazyme_count = len(significant_cazymes)
    secreted_cazyme_count = len(
        significant_secreted_cazymes
    )

    return {
        "species": species,
        "significant_de_genes": significant_count,
        "significant_secreted_genes": secreted_count,
        "significant_nonsecreted_genes": (
            significant_count - secreted_count
        ),
        "significant_high_confidence_cazymes": cazyme_count,
        "significant_secreted_cazymes": secreted_cazyme_count,
        "significant_nonsecreted_cazymes": len(
            nonsecreted_cazymes
        ),
        "significant_secreted_non_cazymes": len(
            secreted_non_cazymes
        ),
        "significant_nonsecreted_non_cazymes": len(
            significant_nonsecreted_non_cazymes
        ),
        "secreted_fraction_of_significant_de": (
            secreted_count / significant_count
            if significant_count
            else 0.0
        ),
        "cazyme_fraction_of_significant_de": (
            cazyme_count / significant_count
            if significant_count
            else 0.0
        ),
        "secreted_cazyme_fraction_of_significant_de": (
            secreted_cazyme_count / significant_count
            if significant_count
            else 0.0
        ),
        "secreted_fraction_of_significant_cazymes": (
            secreted_cazyme_count / cazyme_count
            if cazyme_count
            else 0.0
        ),
        "cazyme_fraction_of_significant_secreted": (
            secreted_cazyme_count / secreted_count
            if secreted_count
            else 0.0
        ),
        "input_file": str(input_file),
        "_significant_set": significant,
        "_secreted_set": significant_secreted,
        "_cazyme_set": significant_cazymes,
        "_secreted_cazyme_set": (
            significant_secreted_cazymes
        ),
    }


def audit_species(
    summary: dict[str, object],
) -> list[dict[str, object]]:
    species = str(summary["species"])

    significant = int(summary["significant_de_genes"])
    secreted = int(summary["significant_secreted_genes"])
    cazymes = int(
        summary[
            "significant_high_confidence_cazymes"
        ]
    )
    secreted_cazymes = int(
        summary["significant_secreted_cazymes"]
    )
    nonsecreted_cazymes = int(
        summary["significant_nonsecreted_cazymes"]
    )

    checks = [
        (
            "significant_de_matches_expected",
            significant,
            EXPECTED_SIGNIFICANT[species],
            significant == EXPECTED_SIGNIFICANT[species],
        ),
        (
            "significant_cazymes_match_expected",
            cazymes,
            EXPECTED_SIGNIFICANT_CAZYMES[species],
            cazymes
            == EXPECTED_SIGNIFICANT_CAZYMES[species],
        ),
        (
            "secreted_subset_of_significant",
            secreted,
            f"<= {significant}",
            secreted <= significant,
        ),
        (
            "cazyme_subset_of_significant",
            cazymes,
            f"<= {significant}",
            cazymes <= significant,
        ),
        (
            "secreted_cazymes_subset_of_secreted",
            secreted_cazymes,
            f"<= {secreted}",
            secreted_cazymes <= secreted,
        ),
        (
            "secreted_cazymes_subset_of_cazymes",
            secreted_cazymes,
            f"<= {cazymes}",
            secreted_cazymes <= cazymes,
        ),
        (
            "cazyme_partition",
            secreted_cazymes + nonsecreted_cazymes,
            cazymes,
            (
                secreted_cazymes
                + nonsecreted_cazymes
                == cazymes
            ),
        ),
        (
            "nonsecreted_cazymes_nonnegative",
            nonsecreted_cazymes,
            ">= 0",
            nonsecreted_cazymes >= 0,
        ),
    ]

    return [
        {
            "species": species,
            "check": check,
            "value": value,
            "expected": expected,
            "status": (
                "PASS"
                if passed
                else "FAIL"
            ),
        }
        for check, value, expected, passed in checks
    ]


def public_summary(
    summary: dict[str, object],
) -> dict[str, object]:
    return {
        key: value
        for key, value in summary.items()
        if not key.startswith("_")
    }


def main() -> None:
    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fusarium, fusarium_seen = process_fusarium()
    ophiostoma, ophiostoma_seen = process_ophiostoma()

    summaries_internal = [
        fusarium,
        ophiostoma,
    ]

    summaries = [
        public_summary(summary)
        for summary in summaries_internal
    ]

    audit_rows = [
        row
        for summary in summaries_internal
        for row in audit_species(summary)
    ]

    failed = [
        row
        for row in audit_rows
        if row["status"] == "FAIL"
    ]

    if failed:
        for row in failed:
            print(
                "FAILED:",
                row,
                file=sys.stderr,
            )

        fail(
            "one or more extracellular summary "
            "validation checks failed"
        )

    write_tsv(
        SUMMARY_OUT,
        list(summaries[0]),
        summaries,
    )

    long_metrics = [
        (
            "Significant DE genes",
            "significant_de_genes",
            "all",
        ),
        (
            "Significant secreted genes",
            "significant_secreted_genes",
            "secreted",
        ),
        (
            "Significant high-confidence CAZymes",
            "significant_high_confidence_cazymes",
            "cazyme",
        ),
        (
            "Significant secreted CAZymes",
            "significant_secreted_cazymes",
            "secreted_cazyme",
        ),
        (
            "Significant non-secreted CAZymes",
            "significant_nonsecreted_cazymes",
            "nonsecreted_cazyme",
        ),
    ]

    summary_long: list[dict[str, object]] = []

    for summary in summaries:
        for metric, field, category in long_metrics:
            summary_long.append(
                {
                    "species": summary["species"],
                    "metric": metric,
                    "category": category,
                    "count": summary[field],
                }
            )

    write_tsv(
        SUMMARY_LONG_OUT,
        [
            "species",
            "metric",
            "category",
            "count",
        ],
        summary_long,
    )

    overlap_rows: list[dict[str, object]] = []

    for summary in summaries:
        overlap_rows.extend(
            [
                {
                    "species": summary["species"],
                    "secreted_status": "Secreted",
                    "cazyme_status": "CAZyme",
                    "count": summary[
                        "significant_secreted_cazymes"
                    ],
                },
                {
                    "species": summary["species"],
                    "secreted_status": "Secreted",
                    "cazyme_status": "Non-CAZyme",
                    "count": summary[
                        "significant_secreted_non_cazymes"
                    ],
                },
                {
                    "species": summary["species"],
                    "secreted_status": "Non-secreted",
                    "cazyme_status": "CAZyme",
                    "count": summary[
                        "significant_nonsecreted_cazymes"
                    ],
                },
                {
                    "species": summary["species"],
                    "secreted_status": "Non-secreted",
                    "cazyme_status": "Non-CAZyme",
                    "count": summary[
                        "significant_nonsecreted_non_cazymes"
                    ],
                },
            ]
        )

    write_tsv(
        OVERLAP_LONG_OUT,
        [
            "species",
            "secreted_status",
            "cazyme_status",
            "count",
        ],
        overlap_rows,
    )

    write_tsv(
        AUDIT_OUT,
        [
            "species",
            "check",
            "value",
            "expected",
            "status",
        ],
        audit_rows,
    )

    run_info = [
        {
            "field": "script",
            "value":
                "35_build_extracellular_summary_table.py",
        },
        {
            "field": "run_timestamp_utc",
            "value":
                datetime.now(timezone.utc).isoformat(),
        },
        {
            "field": "repository_root",
            "value": str(ROOT),
        },
        {
            "field": "significance_definition",
            "value": "adjusted P value < 0.05",
        },
        {
            "field": "secreted_definition",
            "value": (
                "confident SignalP-positive gene; "
                "Prediction == SP with validated "
                "gene-level annotation"
            ),
        },
        {
            "field": "cazyme_definition",
            "value": (
                "high-confidence dbCAN-positive gene; "
                "at least one associated protein "
                "supported by >= 2 dbCAN tools"
            ),
        },
        {
            "field": "secreted_cazyme_definition",
            "value": (
                "intersection of significant, "
                "confident SignalP-positive, and "
                "high-confidence dbCAN-positive genes"
            ),
        },
        {
            "field": "fusarium_contrast",
            "value": "interaction versus self/control",
        },
        {
            "field": "ophiostoma_contrast",
            "value": "interaction versus self",
        },
        {
            "field": "fusarium_input_file",
            "value": str(FUSARIUM_INPUT),
        },
        {
            "field": "fusarium_input_sha256",
            "value": sha256(FUSARIUM_INPUT),
        },
        {
            "field": "fusarium_input_rows",
            "value": len(fusarium_seen),
        },
        {
            "field": "ophiostoma_input_file",
            "value": str(OPHIOSTOMA_INPUT),
        },
        {
            "field": "ophiostoma_input_sha256",
            "value": sha256(OPHIOSTOMA_INPUT),
        },
        {
            "field": "ophiostoma_input_rows",
            "value": len(ophiostoma_seen),
        },
        {
            "field": "statistical_testing",
            "value": (
                "none; descriptive set intersections "
                "of validated gene-level annotations "
                "and DESeq2 results"
            ),
        },
    ]

    write_tsv(
        RUN_INFO_OUT,
        [
            "field",
            "value",
        ],
        run_info,
    )

    print()
    print("============================================================")
    print("EXTRACELLULAR RESPONSE SUMMARY COMPLETE")
    print("============================================================")
    print()

    for summary in summaries:
        species = summary["species"]

        print(species)
        print(
            "  Significant DE genes: "
            f"{summary['significant_de_genes']}"
        )
        print(
            "  Significant secreted genes: "
            f"{summary['significant_secreted_genes']}"
        )
        print(
            "  Significant high-confidence CAZymes: "
            f"{summary['significant_high_confidence_cazymes']}"
        )
        print(
            "  Significant secreted CAZymes: "
            f"{summary['significant_secreted_cazymes']}"
        )
        print(
            "  Significant non-secreted CAZymes: "
            f"{summary['significant_nonsecreted_cazymes']}"
        )
        print(
            "  Secreted fraction of significant CAZymes: "
            f"{100 * float(summary['secreted_fraction_of_significant_cazymes']):.1f}%"
        )
        print()

    print("Validation status: PASS")
    print()
    print(f"Summary: {SUMMARY_OUT}")
    print(f"Long:    {SUMMARY_LONG_OUT}")
    print(f"Overlap: {OVERLAP_LONG_OUT}")
    print(f"Audit:   {AUDIT_OUT}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
