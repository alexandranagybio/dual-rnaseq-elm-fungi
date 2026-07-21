#!/usr/bin/env python3
"""
14_secretome_analysis.py

Integrate validated gene-level DESeq2 results with SignalP 6.0 annotations
for Ophiostoma novo-ulmi.

Definitions
-----------
SignalP-positive candidate:
    signalp_is_sp == TRUE

Significant DE gene:
    padj < 0.05

Strong DE gene:
    padj < 0.05 and abs(log2FoldChange) > 1

The script preserves both DE definitions and does not perform a new
statistical test.

Run from the repository root:
    python workflow/scripts/14_secretome_analysis.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

ANNOTATION = (
    ROOT
    / "results/ophiostoma/functional_annotation/"
      "ophiostoma_functional_annotation.tsv"
)

DESEQ2_TABLES = ROOT / "results/ophiostoma/deseq2_results/tables"

OUTDIR = ROOT / "results/ophiostoma/secretome_analysis"
TABLES_OUT = OUTDIR / "tables"
DIAGNOSTICS_OUT = OUTDIR / "diagnostics"

CONTRASTS = (
    "interaction_vs_self",
    "interaction_vs_onu",
    "onu_vs_self",
)

ALPHA = 0.05
STRONG_ABS_LOG2FC = 1.0

EXPECTED_ANNOTATION_GENES = 8640
EXPECTED_TESTED_GENES = 8560

REQUIRED_DE_COLUMNS = {
    "gene_id",
    "contrast",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}

REQUIRED_ANNOTATION_COLUMNS = {
    "gene_id",
    "mrna_id",
    "signalp_prediction",
    "signalp_is_sp",
    "signalp_other_score",
    "signalp_sp_score",
    "signalp_cs_position",
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
        rows = list(reader)
        return list(reader.fieldnames), rows


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


def parse_float(value: str, field: str, gene_id: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return float("nan")
    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid {field} for {gene_id}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return value != value


def read_annotation() -> tuple[list[str], dict[str, dict[str, str]]]:
    fields, rows = read_tsv(ANNOTATION)
    missing = REQUIRED_ANNOTATION_COLUMNS - set(fields)
    if missing:
        fail(f"annotation table is missing columns: {sorted(missing)}")

    annotation: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = row["gene_id"]
        if gene_id in annotation:
            fail(f"duplicated annotation gene ID: {gene_id}")
        if row["signalp_is_sp"] not in {"TRUE", "FALSE"}:
            fail(
                f"invalid signalp_is_sp value for {gene_id}: "
                f"{row['signalp_is_sp']}"
            )
        annotation[gene_id] = row

    if len(annotation) != EXPECTED_ANNOTATION_GENES:
        fail(
            f"expected {EXPECTED_ANNOTATION_GENES} annotation genes, "
            f"observed {len(annotation)}"
        )

    return fields, annotation


def analyse_contrast(
    contrast: str,
    annotation_fields: list[str],
    annotation: dict[str, dict[str, str]],
) -> dict[str, object]:
    input_path = DESEQ2_TABLES / f"{contrast}_all_genes.tsv"
    de_fields, de_rows = read_tsv(input_path)

    missing = REQUIRED_DE_COLUMNS - set(de_fields)
    if missing:
        fail(f"{input_path} is missing columns: {sorted(missing)}")

    seen: set[str] = set()
    annotated_rows: list[dict[str, object]] = []

    for de_row in de_rows:
        gene_id = de_row["gene_id"]

        if gene_id in seen:
            fail(f"duplicated DESeq2 gene ID in {contrast}: {gene_id}")
        seen.add(gene_id)

        if gene_id not in annotation:
            fail(
                f"DESeq2 gene absent from functional annotation in "
                f"{contrast}: {gene_id}"
            )

        row_contrast = de_row["contrast"]
        if row_contrast != contrast:
            fail(
                f"unexpected contrast value for {gene_id}: "
                f"{row_contrast}; expected {contrast}"
            )

        padj = parse_float(de_row["padj"], "padj", gene_id)
        log2fc = parse_float(
            de_row["log2FoldChange"], "log2FoldChange", gene_id
        )

        significant = not is_nan(padj) and padj < ALPHA
        strong = (
            significant
            and not is_nan(log2fc)
            and abs(log2fc) > STRONG_ABS_LOG2FC
        )
        signalp_positive = (
            annotation[gene_id]["signalp_is_sp"] == "TRUE"
        )

        if is_nan(log2fc):
            direction = "NA"
        elif log2fc > 0:
            direction = "up"
        elif log2fc < 0:
            direction = "down"
        else:
            direction = "unchanged"

        merged: dict[str, object] = dict(de_row)
        for field in annotation_fields:
            if field != "gene_id":
                merged[field] = annotation[gene_id][field]

        merged["de_significant_padj_lt_0.05"] = str(significant).upper()
        merged["de_strong_padj_lt_0.05_abs_log2fc_gt_1"] = (
            str(strong).upper()
        )
        merged["de_direction"] = direction
        merged["signalp_positive_candidate"] = (
            str(signalp_positive).upper()
        )
        merged["significant_signalp_candidate"] = (
            str(significant and signalp_positive).upper()
        )
        merged["strong_signalp_candidate"] = (
            str(strong and signalp_positive).upper()
        )

        annotated_rows.append(merged)

    if len(seen) != EXPECTED_TESTED_GENES:
        fail(
            f"{contrast}: expected {EXPECTED_TESTED_GENES} tested genes, "
            f"observed {len(seen)}"
        )

    output_fields = (
        de_fields
        + [f for f in annotation_fields if f != "gene_id"]
        + [
            "de_significant_padj_lt_0.05",
            "de_strong_padj_lt_0.05_abs_log2fc_gt_1",
            "de_direction",
            "signalp_positive_candidate",
            "significant_signalp_candidate",
            "strong_signalp_candidate",
        ]
    )

    signalp_rows = [
        row for row in annotated_rows
        if row["signalp_positive_candidate"] == "TRUE"
    ]
    significant_signalp = [
        row for row in annotated_rows
        if row["significant_signalp_candidate"] == "TRUE"
    ]
    strong_signalp = [
        row for row in annotated_rows
        if row["strong_signalp_candidate"] == "TRUE"
    ]

    significant_up = [
        row for row in significant_signalp
        if row["de_direction"] == "up"
    ]
    significant_down = [
        row for row in significant_signalp
        if row["de_direction"] == "down"
    ]
    strong_up = [
        row for row in strong_signalp
        if row["de_direction"] == "up"
    ]
    strong_down = [
        row for row in strong_signalp
        if row["de_direction"] == "down"
    ]

    # Stable, biologically useful ordering.
    def sort_key(row: dict[str, object]) -> tuple[float, float, str]:
        padj = parse_float(str(row["padj"]), "padj", str(row["gene_id"]))
        lfc = parse_float(
            str(row["log2FoldChange"]),
            "log2FoldChange",
            str(row["gene_id"]),
        )
        padj_key = padj if not is_nan(padj) else float("inf")
        lfc_key = -abs(lfc) if not is_nan(lfc) else 0.0
        return (padj_key, lfc_key, str(row["gene_id"]))

    for subset in (
        annotated_rows,
        signalp_rows,
        significant_signalp,
        strong_signalp,
        significant_up,
        significant_down,
        strong_up,
        strong_down,
    ):
        subset.sort(key=sort_key)

    write_tsv(
        TABLES_OUT / f"{contrast}_all_genes_with_signalp.tsv",
        output_fields,
        annotated_rows,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_all.tsv",
        output_fields,
        signalp_rows,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_significant.tsv",
        output_fields,
        significant_signalp,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_strong.tsv",
        output_fields,
        strong_signalp,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_significant_up.tsv",
        output_fields,
        significant_up,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_significant_down.tsv",
        output_fields,
        significant_down,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_strong_up.tsv",
        output_fields,
        strong_up,
    )
    write_tsv(
        TABLES_OUT / f"{contrast}_signalp_positive_strong_down.tsv",
        output_fields,
        strong_down,
    )

    tested_signalp = len(signalp_rows)
    tested_non_signalp = len(annotated_rows) - tested_signalp
    significant_total = sum(
        row["de_significant_padj_lt_0.05"] == "TRUE"
        for row in annotated_rows
    )
    strong_total = sum(
        row["de_strong_padj_lt_0.05_abs_log2fc_gt_1"] == "TRUE"
        for row in annotated_rows
    )

    return {
        "contrast": contrast,
        "tested_genes": len(annotated_rows),
        "signalp_positive_tested_genes": tested_signalp,
        "signalp_negative_tested_genes": tested_non_signalp,
        "significant_de_genes_padj_lt_0.05": significant_total,
        "significant_signalp_positive_genes": len(significant_signalp),
        "significant_signalp_positive_up": len(significant_up),
        "significant_signalp_positive_down": len(significant_down),
        "strong_de_genes_padj_lt_0.05_abs_log2fc_gt_1": strong_total,
        "strong_signalp_positive_genes": len(strong_signalp),
        "strong_signalp_positive_up": len(strong_up),
        "strong_signalp_positive_down": len(strong_down),
        "signalp_positive_fraction_of_tested": (
            tested_signalp / len(annotated_rows)
        ),
        "signalp_positive_fraction_of_significant": (
            len(significant_signalp) / significant_total
            if significant_total else 0.0
        ),
        "signalp_positive_fraction_of_strong": (
            len(strong_signalp) / strong_total
            if strong_total else 0.0
        ),
        "input_file": str(input_path),
    }


def main() -> None:
    require_file(ANNOTATION)
    TABLES_OUT.mkdir(parents=True, exist_ok=True)
    DIAGNOSTICS_OUT.mkdir(parents=True, exist_ok=True)

    annotation_fields, annotation = read_annotation()
    complete_sp = sum(
        row["signalp_is_sp"] == "TRUE"
        for row in annotation.values()
    )
    if complete_sp != 591:
        fail(
            f"expected 591 SignalP-positive proteins in complete annotation, "
            f"observed {complete_sp}"
        )

    summaries = [
        analyse_contrast(
            contrast,
            annotation_fields,
            annotation,
        )
        for contrast in CONTRASTS
    ]

    summary_fields = list(summaries[0])
    write_tsv(
        TABLES_OUT / "secretome_summary.tsv",
        summary_fields,
        summaries,
    )

    tested_sets: dict[str, set[str]] = {}
    for contrast in CONTRASTS:
        _, rows = read_tsv(
            DESEQ2_TABLES / f"{contrast}_all_genes.tsv"
        )
        tested_sets[contrast] = {row["gene_id"] for row in rows}

    shared_tested = set.intersection(*tested_sets.values())
    union_tested = set.union(*tested_sets.values())
    excluded_from_testing = set(annotation) - union_tested

    diagnostics = [
        {
            "check": "annotation_gene_count",
            "value": len(annotation),
            "expected": EXPECTED_ANNOTATION_GENES,
            "status": "PASS",
        },
        {
            "check": "complete_annotation_signalp_positive",
            "value": complete_sp,
            "expected": 591,
            "status": "PASS",
        },
        {
            "check": "tested_gene_count_per_contrast",
            "value": ",".join(
                str(len(tested_sets[c])) for c in CONTRASTS
            ),
            "expected": EXPECTED_TESTED_GENES,
            "status": "PASS",
        },
        {
            "check": "tested_gene_sets_identical_across_contrasts",
            "value": len(shared_tested),
            "expected": len(union_tested),
            "status": (
                "PASS"
                if len(shared_tested) == len(union_tested)
                else "FAIL"
            ),
        },
        {
            "check": "genes_excluded_from_deseq2_testing",
            "value": len(excluded_from_testing),
            "expected": EXPECTED_ANNOTATION_GENES - EXPECTED_TESTED_GENES,
            "status": (
                "PASS"
                if len(excluded_from_testing)
                == EXPECTED_ANNOTATION_GENES - EXPECTED_TESTED_GENES
                else "FAIL"
            ),
        },
    ]

    if any(row["status"] == "FAIL" for row in diagnostics):
        fail("one or more cross-contrast validation checks failed")

    write_tsv(
        DIAGNOSTICS_OUT / "validation.tsv",
        ["check", "value", "expected", "status"],
        diagnostics,
    )

    excluded_rows = [
        {
            "gene_id": gene,
            "mrna_id": annotation[gene]["mrna_id"],
            "signalp_prediction": annotation[gene]["signalp_prediction"],
            "signalp_is_sp": annotation[gene]["signalp_is_sp"],
        }
        for gene in sorted(excluded_from_testing)
    ]
    write_tsv(
        DIAGNOSTICS_OUT / "genes_excluded_from_deseq2_testing.tsv",
        ["gene_id", "mrna_id", "signalp_prediction", "signalp_is_sp"],
        excluded_rows,
    )

    run_info = [
        {"field": "script", "value": "14_secretome_analysis.py"},
        {
            "field": "run_timestamp_utc",
            "value": datetime.now(timezone.utc).isoformat(),
        },
        {"field": "repository_root", "value": str(ROOT)},
        {"field": "annotation_file", "value": str(ANNOTATION)},
        {"field": "annotation_sha256", "value": sha256(ANNOTATION)},
        {
            "field": "signalp_definition",
            "value": "signalp_is_sp == TRUE",
        },
        {
            "field": "significant_de_definition",
            "value": "padj < 0.05",
        },
        {
            "field": "strong_de_definition",
            "value": "padj < 0.05 and abs(log2FoldChange) > 1",
        },
        {
            "field": "interpretation",
            "value": (
                "SignalP-positive candidates for classical Sec/SPI secretion; "
                "not a membrane-filtered complete secretome"
            ),
        },
    ]
    for contrast in CONTRASTS:
        path = DESEQ2_TABLES / f"{contrast}_all_genes.tsv"
        run_info.append(
            {
                "field": f"{contrast}_input_sha256",
                "value": sha256(path),
            }
        )

    write_tsv(
        OUTDIR / "run_info.tsv",
        ["field", "value"],
        run_info,
    )

    print("Ophiostoma SignalP-positive candidate analysis")
    print()
    print(f"Complete annotation genes: {len(annotation)}")
    print(f"Complete SignalP-positive proteins: {complete_sp}")
    print(f"Genes tested per contrast: {EXPECTED_TESTED_GENES}")
    print()
    for row in summaries:
        print(row["contrast"])
        print(
            "  Significant DE SignalP-positive: "
            f"{row['significant_signalp_positive_genes']} "
            f"(up {row['significant_signalp_positive_up']}, "
            f"down {row['significant_signalp_positive_down']})"
        )
        print(
            "  Strong DE SignalP-positive: "
            f"{row['strong_signalp_positive_genes']} "
            f"(up {row['strong_signalp_positive_up']}, "
            f"down {row['strong_signalp_positive_down']})"
        )
        print()
    print("Validation status: PASS")
    print(f"Output directory: {OUTDIR}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
