#!/usr/bin/env python3
"""
23_run_fusarium_secretome_analysis.py

Generate gene-level Fusarium secretome tables by integrating the validated
SignalP 6.0 classifications already present in the Fusarium publication
annotation table.

Primary secretome definition
-----------------------------
A Trinity gene is classified as a confident secretome candidate when at least
one mapped TransDecoder protein has a confident SignalP secretion prediction:

    signalp_confident_positive == TRUE

Secondary diagnostic definition
-------------------------------
The broader raw SignalP-positive set is retained separately:

    signalp_raw_positive == TRUE

Differential-expression definitions
-----------------------------------
Significant:
    raw_padj < 0.05

Strong:
    raw_padj < 0.05 and abs(shrunk_log2FoldChange) > 1

No new statistical test is performed.

Run from the repository root:

    python workflow/scripts/23_run_fusarium_secretome_analysis.py
"""

from __future__ import annotations

import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

INPUT = (
    ROOT
    / "results/fusarium/publication_annotation/tables/"
      "fusarium_publication_annotation_full.tsv"
)

OUTDIR = ROOT / "results/fusarium/secretome_analysis"
TABLES_OUT = OUTDIR / "tables"
DIAGNOSTICS_OUT = OUTDIR / "diagnostics"

RUN_INFO = OUTDIR / "run_info.tsv"
CHECKSUMS = OUTDIR / "checksums.sha256"
VALIDATION = DIAGNOSTICS_OUT / "validation.tsv"

ALPHA = 0.05
STRONG_ABS_SHRUNK_LOG2FC = 1.0

EXPECTED_TOTAL_GENES = 15192
EXPECTED_DESEQ2_GENES = 10414
EXPECTED_GENES_WITH_PADJ = 9000
EXPECTED_SIGNIFICANT_GENES = 2973

REQUIRED_COLUMNS = {
    "gene_id",
    "in_deseq2_dataset",
    "padj_available",
    "deseq2_status",
    "regulation",
    "significant_padj_lt_0.05",
    "significant_raw_abs_lfc_gt_1",
    "significant_shrunk_abs_lfc_gt_1",
    "raw_baseMean",
    "raw_log2FoldChange",
    "raw_lfcSE",
    "raw_stat",
    "raw_pvalue",
    "raw_padj",
    "shrunk_log2FoldChange",
    "shrunk_lfcSE",
    "has_predicted_protein",
    "predicted_protein_count",
    "protein_ids",
    "signalp_raw_positive_protein_count",
    "signalp_raw_positive",
    "signalp_confident_positive_protein_count",
    "signalp_confident_positive",
    "signalp_ambiguous_protein_count",
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


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require_file(path)

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            fail(f"no header found in {path}")

        return list(reader.fieldnames), list(reader)


def write_tsv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

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


def parse_bool(value: str, field: str, gene_id: str) -> bool:
    if value == "TRUE":
        return True
    if value == "FALSE":
        return False

    fail(f"invalid Boolean value for {gene_id}, {field}: {value}")


def parse_int(value: str, field: str, gene_id: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid integer for {gene_id}, {field}: {value}"
        ) from exc


def parse_float(value: str, field: str, gene_id: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return float("nan")

    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid numeric value for {gene_id}, "
            f"{field}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return value != value


def sort_key(row: dict[str, object]) -> tuple[float, float, str]:
    gene_id = str(row["gene_id"])

    padj = parse_float(
        str(row["raw_padj"]),
        "raw_padj",
        gene_id,
    )

    shrunk_lfc = parse_float(
        str(row["shrunk_log2FoldChange"]),
        "shrunk_log2FoldChange",
        gene_id,
    )

    padj_key = padj if not is_nan(padj) else float("inf")
    lfc_key = (
        -abs(shrunk_lfc)
        if not is_nan(shrunk_lfc)
        else 0.0
    )

    return padj_key, lfc_key, gene_id


def count_direction(rows: list[dict[str, object]], direction: str) -> int:
    expected_label = {
        "up": "upregulated",
        "down": "downregulated",
    }

    if direction not in expected_label:
        fail(f"unsupported direction requested: {direction}")

    return sum(
        str(row["regulation"]) == expected_label[direction]
        for row in rows
    )


def main() -> None:
    fields, rows = read_tsv(INPUT)

    missing = REQUIRED_COLUMNS - set(fields)
    if missing:
        fail(
            "publication annotation table is missing columns: "
            f"{sorted(missing)}"
        )

    if len(rows) != EXPECTED_TOTAL_GENES:
        fail(
            f"expected {EXPECTED_TOTAL_GENES} total genes, "
            f"observed {len(rows)}"
        )

    seen: set[str] = set()
    annotated_rows: list[dict[str, object]] = []

    for row in rows:
        gene_id = row["gene_id"]

        if not gene_id:
            fail("encountered empty gene_id")

        if gene_id in seen:
            fail(f"duplicated gene ID: {gene_id}")

        seen.add(gene_id)

        in_deseq2 = parse_bool(
            row["in_deseq2_dataset"],
            "in_deseq2_dataset",
            gene_id,
        )

        padj_available = parse_bool(
            row["padj_available"],
            "padj_available",
            gene_id,
        )

        significant = parse_bool(
            row["significant_padj_lt_0.05"],
            "significant_padj_lt_0.05",
            gene_id,
        )

        significant_raw_lfc1 = parse_bool(
            row["significant_raw_abs_lfc_gt_1"],
            "significant_raw_abs_lfc_gt_1",
            gene_id,
        )

        significant_shrunk_lfc1 = parse_bool(
            row["significant_shrunk_abs_lfc_gt_1"],
            "significant_shrunk_abs_lfc_gt_1",
            gene_id,
        )

        raw_positive = parse_bool(
            row["signalp_raw_positive"],
            "signalp_raw_positive",
            gene_id,
        )

        confident_positive = parse_bool(
            row["signalp_confident_positive"],
            "signalp_confident_positive",
            gene_id,
        )

        raw_positive_count = parse_int(
            row["signalp_raw_positive_protein_count"],
            "signalp_raw_positive_protein_count",
            gene_id,
        )

        confident_positive_count = parse_int(
            row["signalp_confident_positive_protein_count"],
            "signalp_confident_positive_protein_count",
            gene_id,
        )

        ambiguous_count = parse_int(
            row["signalp_ambiguous_protein_count"],
            "signalp_ambiguous_protein_count",
            gene_id,
        )

        if raw_positive != (raw_positive_count >= 1):
            fail(
                f"inconsistent raw SignalP flag for {gene_id}: "
                f"count={raw_positive_count}, flag={raw_positive}"
            )

        if confident_positive != (confident_positive_count >= 1):
            fail(
                f"inconsistent confident SignalP flag for {gene_id}: "
                f"count={confident_positive_count}, "
                f"flag={confident_positive}"
            )

        if confident_positive_count > raw_positive_count:
            fail(
                f"confident SignalP count exceeds raw-positive count "
                f"for {gene_id}"
            )

        if confident_positive and not raw_positive:
            fail(
                f"confident SignalP-positive gene is not raw positive: "
                f"{gene_id}"
            )

        if significant and not padj_available:
            fail(
                f"significant gene lacks adjusted P value: {gene_id}"
            )

        if significant and not in_deseq2:
            fail(
                f"significant gene is absent from DESeq2 dataset: "
                f"{gene_id}"
            )

        if significant_raw_lfc1 and not significant:
            fail(
                f"raw strong-DE flag without significance: {gene_id}"
            )

        if significant_shrunk_lfc1 and not significant:
            fail(
                f"shrunken strong-DE flag without significance: "
                f"{gene_id}"
            )

        merged: dict[str, object] = dict(row)

        merged["secretome_primary_confident"] = str(
            confident_positive
        ).upper()

        merged["secretome_secondary_raw"] = str(
            raw_positive
        ).upper()

        merged["significant_confident_secretome"] = str(
            significant and confident_positive
        ).upper()

        merged["strong_raw_lfc_confident_secretome"] = str(
            significant_raw_lfc1 and confident_positive
        ).upper()

        merged["strong_shrunk_lfc_confident_secretome"] = str(
            significant_shrunk_lfc1 and confident_positive
        ).upper()

        merged["signalp_ambiguous_gene"] = str(
            ambiguous_count >= 1
        ).upper()

        annotated_rows.append(merged)

    output_fields = fields + [
        "secretome_primary_confident",
        "secretome_secondary_raw",
        "significant_confident_secretome",
        "strong_raw_lfc_confident_secretome",
        "strong_shrunk_lfc_confident_secretome",
        "signalp_ambiguous_gene",
    ]

    complete_raw = [
        row
        for row in annotated_rows
        if row["secretome_secondary_raw"] == "TRUE"
    ]

    complete_confident = [
        row
        for row in annotated_rows
        if row["secretome_primary_confident"] == "TRUE"
    ]

    tested_confident = [
        row
        for row in complete_confident
        if row["in_deseq2_dataset"] == "TRUE"
    ]

    padj_confident = [
        row
        for row in tested_confident
        if row["padj_available"] == "TRUE"
    ]

    significant_confident = [
        row
        for row in complete_confident
        if row["significant_confident_secretome"] == "TRUE"
    ]

    strong_raw_confident = [
        row
        for row in complete_confident
        if row["strong_raw_lfc_confident_secretome"] == "TRUE"
    ]

    strong_shrunk_confident = [
        row
        for row in complete_confident
        if row["strong_shrunk_lfc_confident_secretome"] == "TRUE"
    ]

    ambiguous_genes = [
        row
        for row in annotated_rows
        if row["signalp_ambiguous_gene"] == "TRUE"
    ]

    for subset in (
        complete_raw,
        complete_confident,
        tested_confident,
        padj_confident,
        significant_confident,
        strong_raw_confident,
        strong_shrunk_confident,
        ambiguous_genes,
    ):
        subset.sort(key=sort_key)

    output_paths = {
        "complete_raw": (
            TABLES_OUT
            / "fusarium_secretome_raw_positive_complete.tsv"
        ),
        "complete_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_complete.tsv"
        ),
        "tested_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_deseq2_dataset.tsv"
        ),
        "padj_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_padj_available.tsv"
        ),
        "significant_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_significant.tsv"
        ),
        "strong_raw_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_significant_raw_lfc1.tsv"
        ),
        "strong_shrunk_confident": (
            TABLES_OUT
            / "fusarium_secretome_confident_significant_shrunk_lfc1.tsv"
        ),
        "ambiguous": (
            DIAGNOSTICS_OUT
            / "fusarium_signalp_ambiguous_genes.tsv"
        ),
    }

    for key, path in output_paths.items():
        if key == "complete_raw":
            subset = complete_raw
        elif key == "complete_confident":
            subset = complete_confident
        elif key == "tested_confident":
            subset = tested_confident
        elif key == "padj_confident":
            subset = padj_confident
        elif key == "significant_confident":
            subset = significant_confident
        elif key == "strong_raw_confident":
            subset = strong_raw_confident
        elif key == "strong_shrunk_confident":
            subset = strong_shrunk_confident
        elif key == "ambiguous":
            subset = ambiguous_genes
        else:
            fail(f"unhandled output key: {key}")

        write_tsv(path, output_fields, subset)

    summary_rows: list[dict[str, object]] = [
        {
            "category": "all_annotation_genes",
            "total": len(annotated_rows),
            "up": "",
            "down": "",
        },
        {
            "category": "genes_in_deseq2_dataset",
            "total": sum(
                row["in_deseq2_dataset"] == "TRUE"
                for row in annotated_rows
            ),
            "up": "",
            "down": "",
        },
        {
            "category": "genes_with_padj",
            "total": sum(
                row["padj_available"] == "TRUE"
                for row in annotated_rows
            ),
            "up": "",
            "down": "",
        },
        {
            "category": "all_significant_genes",
            "total": sum(
                row["significant_padj_lt_0.05"] == "TRUE"
                for row in annotated_rows
            ),
            "up": sum(
                row["significant_padj_lt_0.05"] == "TRUE"
                and row["regulation"] == "upregulated"
                for row in annotated_rows
            ),
            "down": sum(
                row["significant_padj_lt_0.05"] == "TRUE"
                and row["regulation"] == "downregulated"
                for row in annotated_rows
            ),
        },
        {
            "category": "signalp_raw_positive_genes",
            "total": len(complete_raw),
            "up": "",
            "down": "",
        },
        {
            "category": "signalp_confident_positive_genes",
            "total": len(complete_confident),
            "up": "",
            "down": "",
        },
        {
            "category": "confident_secretome_in_deseq2_dataset",
            "total": len(tested_confident),
            "up": count_direction(tested_confident, "up"),
            "down": count_direction(tested_confident, "down"),
        },
        {
            "category": "confident_secretome_with_padj",
            "total": len(padj_confident),
            "up": count_direction(padj_confident, "up"),
            "down": count_direction(padj_confident, "down"),
        },
        {
            "category": "significant_confident_secretome",
            "total": len(significant_confident),
            "up": count_direction(significant_confident, "up"),
            "down": count_direction(significant_confident, "down"),
        },
        {
            "category": "strong_raw_lfc_confident_secretome",
            "total": len(strong_raw_confident),
            "up": count_direction(strong_raw_confident, "up"),
            "down": count_direction(strong_raw_confident, "down"),
        },
        {
            "category": "strong_shrunk_lfc_confident_secretome",
            "total": len(strong_shrunk_confident),
            "up": count_direction(strong_shrunk_confident, "up"),
            "down": count_direction(strong_shrunk_confident, "down"),
        },
        {
            "category": "signalp_ambiguous_genes",
            "total": len(ambiguous_genes),
            "up": count_direction(ambiguous_genes, "up"),
            "down": count_direction(ambiguous_genes, "down"),
        },
    ]

    summary_path = TABLES_OUT / "fusarium_secretome_summary.tsv"

    write_tsv(
        summary_path,
        ["category", "total", "up", "down"],
        summary_rows,
    )

    observed_deseq2 = sum(
        row["in_deseq2_dataset"] == "TRUE"
        for row in annotated_rows
    )

    observed_padj = sum(
        row["padj_available"] == "TRUE"
        for row in annotated_rows
    )

    observed_significant = sum(
        row["significant_padj_lt_0.05"] == "TRUE"
        for row in annotated_rows
    )

    validation_rows = [
        {
            "check": "total_gene_rows",
            "value": len(annotated_rows),
            "expected": EXPECTED_TOTAL_GENES,
            "status": (
                "PASS"
                if len(annotated_rows) == EXPECTED_TOTAL_GENES
                else "FAIL"
            ),
        },
        {
            "check": "unique_gene_ids",
            "value": len(seen),
            "expected": EXPECTED_TOTAL_GENES,
            "status": (
                "PASS"
                if len(seen) == EXPECTED_TOTAL_GENES
                else "FAIL"
            ),
        },
        {
            "check": "genes_in_deseq2_dataset",
            "value": observed_deseq2,
            "expected": EXPECTED_DESEQ2_GENES,
            "status": (
                "PASS"
                if observed_deseq2 == EXPECTED_DESEQ2_GENES
                else "FAIL"
            ),
        },
        {
            "check": "genes_with_padj",
            "value": observed_padj,
            "expected": EXPECTED_GENES_WITH_PADJ,
            "status": (
                "PASS"
                if observed_padj == EXPECTED_GENES_WITH_PADJ
                else "FAIL"
            ),
        },
        {
            "check": "significant_genes",
            "value": observed_significant,
            "expected": EXPECTED_SIGNIFICANT_GENES,
            "status": (
                "PASS"
                if observed_significant == EXPECTED_SIGNIFICANT_GENES
                else "FAIL"
            ),
        },
        {
            "check": "raw_signalp_gene_flag_matches_counts",
            "value": len(annotated_rows),
            "expected": len(annotated_rows),
            "status": "PASS",
        },
        {
            "check": "confident_signalp_gene_flag_matches_counts",
            "value": len(annotated_rows),
            "expected": len(annotated_rows),
            "status": "PASS",
        },
        {
            "check": "confident_signalp_subset_of_raw_signalp",
            "value": len(complete_confident),
            "expected": f"<= {len(complete_raw)}",
            "status": (
                "PASS"
                if len(complete_confident) <= len(complete_raw)
                else "FAIL"
            ),
        },
        {
            "check": "significant_secretome_subset_of_confident_secretome",
            "value": len(significant_confident),
            "expected": f"<= {len(complete_confident)}",
            "status": (
                "PASS"
                if len(significant_confident)
                <= len(complete_confident)
                else "FAIL"
            ),
        },
        {
            "check": "strong_shrunk_secretome_subset_of_significant",
            "value": len(strong_shrunk_confident),
            "expected": f"<= {len(significant_confident)}",
            "status": (
                "PASS"
                if len(strong_shrunk_confident)
                <= len(significant_confident)
                else "FAIL"
            ),
        },
    ]

    write_tsv(
        VALIDATION,
        ["check", "value", "expected", "status"],
        validation_rows,
    )

    failed = [
        row
        for row in validation_rows
        if row["status"] != "PASS"
    ]

    if failed:
        fail(
            "one or more validation checks failed; inspect "
            f"{VALIDATION}"
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)

    run_info_rows = [
        {
            "parameter": "script",
            "value": "workflow/scripts/23_run_fusarium_secretome_analysis.py",
        },
        {
            "parameter": "run_utc",
            "value": datetime.now(timezone.utc).isoformat(),
        },
        {
            "parameter": "input",
            "value": str(INPUT.relative_to(ROOT)),
        },
        {
            "parameter": "primary_secretome_definition",
            "value": "signalp_confident_positive == TRUE",
        },
        {
            "parameter": "secondary_secretome_definition",
            "value": "signalp_raw_positive == TRUE",
        },
        {
            "parameter": "significance_definition",
            "value": "raw_padj < 0.05",
        },
        {
            "parameter": "strong_primary_definition",
            "value": (
                "raw_padj < 0.05 and "
                "abs(shrunk_log2FoldChange) > 1"
            ),
        },
        {
            "parameter": "new_statistical_test",
            "value": "FALSE",
        },
    ]

    write_tsv(
        RUN_INFO,
        ["parameter", "value"],
        run_info_rows,
    )

    checksum_paths = [
        INPUT,
        summary_path,
        VALIDATION,
        RUN_INFO,
        *output_paths.values(),
    ]

    with CHECKSUMS.open("w") as handle:
        for path in checksum_paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(ROOT)}\n"
            )

    print("Fusarium secretome analysis completed successfully.")
    print(f"Total annotation genes:          {len(annotated_rows):,}")
    print(f"Raw SignalP-positive genes:      {len(complete_raw):,}")
    print(f"Confident secretome genes:       {len(complete_confident):,}")
    print(f"Secretome genes in DESeq2:       {len(tested_confident):,}")
    print(f"Secretome genes with padj:       {len(padj_confident):,}")
    print(
        "Significant secretome genes:    "
        f"{len(significant_confident):,} "
        f"(up={count_direction(significant_confident, 'up'):,}, "
        f"down={count_direction(significant_confident, 'down'):,})"
    )
    print(
        "Strong secretome, raw LFC:      "
        f"{len(strong_raw_confident):,} "
        f"(up={count_direction(strong_raw_confident, 'up'):,}, "
        f"down={count_direction(strong_raw_confident, 'down'):,})"
    )
    print(
        "Strong secretome, shrunk LFC:   "
        f"{len(strong_shrunk_confident):,} "
        f"(up={count_direction(strong_shrunk_confident, 'up'):,}, "
        f"down={count_direction(strong_shrunk_confident, 'down'):,})"
    )
    print(f"Ambiguous SignalP genes:         {len(ambiguous_genes):,}")
    print(f"Validation:                      {VALIDATION}")
    print(f"Summary:                         {summary_path}")


if __name__ == "__main__":
    main()
