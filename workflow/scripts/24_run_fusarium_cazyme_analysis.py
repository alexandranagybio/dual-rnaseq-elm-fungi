#!/usr/bin/env python3
"""
24_run_fusarium_cazyme_analysis.py

Generate validated gene-level CAZyme subsets and summaries for the rebuilt
Fusarium RNA-seq analysis.

Primary CAZyme definition
-------------------------
High-confidence CAZyme gene:
    dbcan_high_confidence == TRUE

This corresponds to at least one predicted protein associated with the Trinity
gene being supported by at least two dbCAN methods.

Secondary definition
--------------------
Any-hit CAZyme candidate:
    dbcan_any_hit == TRUE

Differential-expression definitions
-----------------------------------
Significant gene:
    significant_padj_lt_0.05 == TRUE

Strong gene using the raw effect:
    significant_raw_abs_lfc_gt_1 == TRUE

Strong gene using the shrunken effect:
    significant_shrunk_abs_lfc_gt_1 == TRUE

The script performs no new differential-expression test and no enrichment
test. It operates directly on the validated publication annotation table.

Run from the repository root:

    python workflow/scripts/24_run_fusarium_cazyme_analysis.py
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

INPUT = (
    ROOT
    / "results/fusarium/publication_annotation/tables/"
      "fusarium_publication_annotation_full.tsv"
)

OUTDIR = ROOT / "results/fusarium/cazyme_analysis"
TABLES_OUT = OUTDIR / "tables"
DIAGNOSTICS_OUT = OUTDIR / "diagnostics"

EXPECTED_TOTAL_GENES = 15192
EXPECTED_DESEQ2_GENES = 10414
EXPECTED_PADJ_GENES = 9000
EXPECTED_SIGNIFICANT_GENES = 2973

EXPECTED_ANY_HIT_GENES = 2015
EXPECTED_HIGH_CONFIDENCE_GENES = 420
EXPECTED_HIGH_CONFIDENCE_DESEQ2 = 404
EXPECTED_HIGH_CONFIDENCE_PADJ = 378
EXPECTED_SIGNIFICANT_HIGH_CONFIDENCE = 142
EXPECTED_STRONG_RAW_HIGH_CONFIDENCE = 8
EXPECTED_STRONG_SHRUNK_HIGH_CONFIDENCE = 5

EXPECTED_SIGNIFICANT_UP = 80
EXPECTED_SIGNIFICANT_DOWN = 62
EXPECTED_STRONG_RAW_UP = 8
EXPECTED_STRONG_RAW_DOWN = 0
EXPECTED_STRONG_SHRUNK_UP = 5
EXPECTED_STRONG_SHRUNK_DOWN = 0

BOOLEAN_COLUMNS = {
    "in_deseq2_dataset",
    "padj_available",
    "significant_padj_lt_0.05",
    "significant_raw_abs_lfc_gt_1",
    "significant_shrunk_abs_lfc_gt_1",
    "dbcan_any_hit",
    "dbcan_high_confidence",
}

REQUIRED_COLUMNS = {
    "gene_id",
    "regulation",
    "raw_padj",
    "raw_log2FoldChange",
    "shrunk_log2FoldChange",
    "in_deseq2_dataset",
    "padj_available",
    "significant_padj_lt_0.05",
    "significant_raw_abs_lfc_gt_1",
    "significant_shrunk_abs_lfc_gt_1",
    "dbcan_any_hit_protein_count",
    "dbcan_any_hit",
    "dbcan_high_confidence_protein_count",
    "dbcan_high_confidence",
    "dbcan_max_n_tools",
    "dbcan_ec",
    "dbcan_hmm",
    "dbcan_sub",
    "dbcan_diamond",
    "dbcan_recommended",
    "dbcan_substrate",
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


def validate_boolean(
    value: str,
    field: str,
    gene_id: str,
) -> None:
    if value not in {"TRUE", "FALSE"}:
        fail(
            f"invalid Boolean value in {field} for {gene_id}: {value}"
        )


def parse_int(
    value: str,
    field: str,
    gene_id: str,
) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid integer in {field} for {gene_id}: {value}"
        ) from exc


def parse_float(
    value: str,
    field: str,
    gene_id: str,
) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return float("nan")

    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid numeric value in {field} "
            f"for {gene_id}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return value != value


def sort_key(
    row: dict[str, object],
) -> tuple[float, float, str]:
    gene_id = str(row["gene_id"])

    padj = parse_float(
        str(row["raw_padj"]),
        "raw_padj",
        gene_id,
    )

    lfc = parse_float(
        str(row["shrunk_log2FoldChange"]),
        "shrunk_log2FoldChange",
        gene_id,
    )

    padj_key = padj if not is_nan(padj) else float("inf")
    lfc_key = -abs(lfc) if not is_nan(lfc) else 0.0

    return padj_key, lfc_key, gene_id


def count_regulation(
    rows: list[dict[str, str]],
    label: str,
) -> int:
    return sum(row["regulation"] == label for row in rows)


def extract_cazyme_families(value: str) -> set[str]:
    """
    Extract unique CAZy family identifiers from dbcan_recommended.

    Examples retained:
        AA1_3
        AA3_3
        CBM18
        CE4
        GH5_31
        GT48

    Repeated family assignments within one Trinity gene are counted once.
    """
    if value in {"", "NA", "nan", "NaN"}:
        return set()

    return set(
        re.findall(
            r"\b(?:AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?\b",
            value,
        )
    )


def family_class(family: str) -> str:
    match = re.match(r"^(AA|CBM|CE|GH|GT|PL)", family)

    if match is None:
        return "OTHER"

    return match.group(1)


def make_family_summary(
    all_high_confidence: list[dict[str, str]],
    significant: list[dict[str, str]],
    strong_raw: list[dict[str, str]],
    strong_shrunk: list[dict[str, str]],
) -> list[dict[str, object]]:
    all_counts: Counter[str] = Counter()
    significant_counts: Counter[str] = Counter()
    significant_up_counts: Counter[str] = Counter()
    significant_down_counts: Counter[str] = Counter()
    strong_raw_counts: Counter[str] = Counter()
    strong_shrunk_counts: Counter[str] = Counter()

    for row in all_high_confidence:
        for family in extract_cazyme_families(
            row["dbcan_recommended"]
        ):
            all_counts[family] += 1

    for row in significant:
        families = extract_cazyme_families(
            row["dbcan_recommended"]
        )

        for family in families:
            significant_counts[family] += 1

            if row["regulation"] == "upregulated":
                significant_up_counts[family] += 1
            elif row["regulation"] == "downregulated":
                significant_down_counts[family] += 1

    for row in strong_raw:
        for family in extract_cazyme_families(
            row["dbcan_recommended"]
        ):
            strong_raw_counts[family] += 1

    for row in strong_shrunk:
        for family in extract_cazyme_families(
            row["dbcan_recommended"]
        ):
            strong_shrunk_counts[family] += 1

    families = sorted(
        all_counts,
        key=lambda family: (
            -significant_counts[family],
            -all_counts[family],
            family,
        ),
    )

    return [
        {
            "cazyme_family": family,
            "cazyme_class": family_class(family),
            "high_confidence_genes": all_counts[family],
            "significant_genes": significant_counts[family],
            "significant_up": significant_up_counts[family],
            "significant_down": significant_down_counts[family],
            "strong_raw_lfc_genes": strong_raw_counts[family],
            "strong_shrunk_lfc_genes": strong_shrunk_counts[family],
        }
        for family in families
    ]


def make_class_summary(
    family_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    class_totals: dict[str, Counter[str]] = {}

    for row in family_rows:
        cazyme_class = str(row["cazyme_class"])

        if cazyme_class not in class_totals:
            class_totals[cazyme_class] = Counter()

        totals = class_totals[cazyme_class]

        totals["families"] += 1
        totals["high_confidence_genes"] += int(
            row["high_confidence_genes"]
        )
        totals["significant_genes"] += int(
            row["significant_genes"]
        )
        totals["significant_up"] += int(
            row["significant_up"]
        )
        totals["significant_down"] += int(
            row["significant_down"]
        )
        totals["strong_raw_lfc_genes"] += int(
            row["strong_raw_lfc_genes"]
        )
        totals["strong_shrunk_lfc_genes"] += int(
            row["strong_shrunk_lfc_genes"]
        )

    return [
        {
            "cazyme_class": cazyme_class,
            **dict(class_totals[cazyme_class]),
        }
        for cazyme_class in sorted(class_totals)
    ]


def make_substrate_summary(
    rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    """
    Produce a descriptive substrate table.

    dbCAN substrate predictions may be empty for most genes. Values are split
    on semicolons and deduplicated within each Trinity gene.
    """
    total_counts: Counter[str] = Counter()
    significant_counts: Counter[str] = Counter()
    significant_up_counts: Counter[str] = Counter()
    significant_down_counts: Counter[str] = Counter()

    for row in rows:
        value = row["dbcan_substrate"].strip()

        if value in {"", "NA", "NaN", "nan"}:
            continue

        substrates = {
            item.strip()
            for item in value.split(";")
            if item.strip()
        }

        for substrate in substrates:
            total_counts[substrate] += 1

            if row["significant_padj_lt_0.05"] == "TRUE":
                significant_counts[substrate] += 1

                if row["regulation"] == "upregulated":
                    significant_up_counts[substrate] += 1
                elif row["regulation"] == "downregulated":
                    significant_down_counts[substrate] += 1

    return [
        {
            "substrate": substrate,
            "high_confidence_genes": total_counts[substrate],
            "significant_genes": significant_counts[substrate],
            "significant_up": significant_up_counts[substrate],
            "significant_down": significant_down_counts[substrate],
        }
        for substrate in sorted(
            total_counts,
            key=lambda item: (
                -significant_counts[item],
                -total_counts[item],
                item,
            ),
        )
    ]


def validation_row(
    check: str,
    value: object,
    expected: object,
    passed: bool,
) -> dict[str, object]:
    return {
        "check": check,
        "value": value,
        "expected": expected,
        "status": "PASS" if passed else "FAIL",
    }


def main() -> None:
    require_file(INPUT)

    TABLES_OUT.mkdir(parents=True, exist_ok=True)
    DIAGNOSTICS_OUT.mkdir(parents=True, exist_ok=True)

    fields, rows = read_tsv(INPUT)

    missing = REQUIRED_COLUMNS - set(fields)

    if missing:
        fail(
            "publication annotation table is missing columns: "
            f"{sorted(missing)}"
        )

    if len(rows) != EXPECTED_TOTAL_GENES:
        fail(
            f"expected {EXPECTED_TOTAL_GENES} gene rows, "
            f"observed {len(rows)}"
        )

    seen: set[str] = set()

    any_hit_count_consistency = 0
    high_confidence_count_consistency = 0
    high_confidence_subset_consistency = 0

    for row in rows:
        gene_id = row["gene_id"]

        if not gene_id:
            fail("encountered empty gene_id")

        if gene_id in seen:
            fail(f"duplicated gene_id: {gene_id}")

        seen.add(gene_id)

        for field in BOOLEAN_COLUMNS:
            validate_boolean(row[field], field, gene_id)

        any_count = parse_int(
            row["dbcan_any_hit_protein_count"],
            "dbcan_any_hit_protein_count",
            gene_id,
        )

        high_count = parse_int(
            row["dbcan_high_confidence_protein_count"],
            "dbcan_high_confidence_protein_count",
            gene_id,
        )

        max_tools = parse_int(
            row["dbcan_max_n_tools"],
            "dbcan_max_n_tools",
            gene_id,
        )

        any_flag = row["dbcan_any_hit"] == "TRUE"
        high_flag = row["dbcan_high_confidence"] == "TRUE"

        if any_flag == (any_count > 0):
            any_hit_count_consistency += 1

        if high_flag == (high_count > 0):
            high_confidence_count_consistency += 1

        if not high_flag or any_flag:
            high_confidence_subset_consistency += 1

        if any_flag and max_tools < 1:
            fail(
                f"dbCAN any-hit gene has dbcan_max_n_tools < 1: "
                f"{gene_id}"
            )

        if high_flag and max_tools < 2:
            fail(
                f"high-confidence CAZyme gene has "
                f"dbcan_max_n_tools < 2: {gene_id}"
            )

        if high_count > any_count:
            fail(
                f"high-confidence protein count exceeds any-hit "
                f"protein count for {gene_id}"
            )

    rows.sort(key=sort_key)

    deseq2_rows = [
        row
        for row in rows
        if row["in_deseq2_dataset"] == "TRUE"
    ]

    padj_rows = [
        row
        for row in rows
        if row["padj_available"] == "TRUE"
    ]

    all_significant = [
        row
        for row in rows
        if row["significant_padj_lt_0.05"] == "TRUE"
    ]

    any_hit = [
        row
        for row in rows
        if row["dbcan_any_hit"] == "TRUE"
    ]

    high_confidence = [
        row
        for row in rows
        if row["dbcan_high_confidence"] == "TRUE"
    ]

    high_confidence_deseq2 = [
        row
        for row in high_confidence
        if row["in_deseq2_dataset"] == "TRUE"
    ]

    high_confidence_padj = [
        row
        for row in high_confidence
        if row["padj_available"] == "TRUE"
    ]

    significant_high_confidence = [
        row
        for row in high_confidence
        if row["significant_padj_lt_0.05"] == "TRUE"
    ]

    strong_raw_high_confidence = [
        row
        for row in high_confidence
        if row["significant_raw_abs_lfc_gt_1"] == "TRUE"
    ]

    strong_shrunk_high_confidence = [
        row
        for row in high_confidence
        if row["significant_shrunk_abs_lfc_gt_1"] == "TRUE"
    ]

    significant_up = [
        row
        for row in significant_high_confidence
        if row["regulation"] == "upregulated"
    ]

    significant_down = [
        row
        for row in significant_high_confidence
        if row["regulation"] == "downregulated"
    ]

    strong_raw_up = [
        row
        for row in strong_raw_high_confidence
        if row["regulation"] == "upregulated"
    ]

    strong_raw_down = [
        row
        for row in strong_raw_high_confidence
        if row["regulation"] == "downregulated"
    ]

    strong_shrunk_up = [
        row
        for row in strong_shrunk_high_confidence
        if row["regulation"] == "upregulated"
    ]

    strong_shrunk_down = [
        row
        for row in strong_shrunk_high_confidence
        if row["regulation"] == "downregulated"
    ]

    outputs = {
        "fusarium_cazymes_any_hit_complete.tsv": any_hit,
        "fusarium_cazymes_high_confidence_complete.tsv": (
            high_confidence
        ),
        "fusarium_cazymes_high_confidence_deseq2_dataset.tsv": (
            high_confidence_deseq2
        ),
        "fusarium_cazymes_high_confidence_padj_available.tsv": (
            high_confidence_padj
        ),
        "fusarium_cazymes_high_confidence_significant.tsv": (
            significant_high_confidence
        ),
        "fusarium_cazymes_high_confidence_significant_up.tsv": (
            significant_up
        ),
        "fusarium_cazymes_high_confidence_significant_down.tsv": (
            significant_down
        ),
        "fusarium_cazymes_high_confidence_"
        "significant_raw_lfc1.tsv": (
            strong_raw_high_confidence
        ),
        "fusarium_cazymes_high_confidence_"
        "significant_raw_lfc1_up.tsv": strong_raw_up,
        "fusarium_cazymes_high_confidence_"
        "significant_raw_lfc1_down.tsv": strong_raw_down,
        "fusarium_cazymes_high_confidence_"
        "significant_shrunk_lfc1.tsv": (
            strong_shrunk_high_confidence
        ),
        "fusarium_cazymes_high_confidence_"
        "significant_shrunk_lfc1_up.tsv": strong_shrunk_up,
        "fusarium_cazymes_high_confidence_"
        "significant_shrunk_lfc1_down.tsv": strong_shrunk_down,
    }

    for filename, subset in outputs.items():
        write_tsv(
            TABLES_OUT / filename,
            fields,
            subset,
        )

    summary = [
        {
            "category": "all_annotation_genes",
            "total": len(rows),
            "significant_up": "",
            "significant_down": "",
        },
        {
            "category": "genes_in_deseq2_dataset",
            "total": len(deseq2_rows),
            "significant_up": "",
            "significant_down": "",
        },
        {
            "category": "genes_with_padj",
            "total": len(padj_rows),
            "significant_up": "",
            "significant_down": "",
        },
        {
            "category": "all_significant_genes",
            "total": len(all_significant),
            "significant_up": count_regulation(
                all_significant,
                "upregulated",
            ),
            "significant_down": count_regulation(
                all_significant,
                "downregulated",
            ),
        },
        {
            "category": "dbcan_any_hit_genes",
            "total": len(any_hit),
            "significant_up": "",
            "significant_down": "",
        },
        {
            "category": "high_confidence_cazyme_genes",
            "total": len(high_confidence),
            "significant_up": "",
            "significant_down": "",
        },
        {
            "category": "high_confidence_cazymes_in_deseq2_dataset",
            "total": len(high_confidence_deseq2),
            "significant_up": len(significant_up),
            "significant_down": len(significant_down),
        },
        {
            "category": "high_confidence_cazymes_with_padj",
            "total": len(high_confidence_padj),
            "significant_up": len(significant_up),
            "significant_down": len(significant_down),
        },
        {
            "category": "significant_high_confidence_cazymes",
            "total": len(significant_high_confidence),
            "significant_up": len(significant_up),
            "significant_down": len(significant_down),
        },
        {
            "category": "strong_raw_lfc_high_confidence_cazymes",
            "total": len(strong_raw_high_confidence),
            "significant_up": len(strong_raw_up),
            "significant_down": len(strong_raw_down),
        },
        {
            "category": "strong_shrunk_lfc_high_confidence_cazymes",
            "total": len(strong_shrunk_high_confidence),
            "significant_up": len(strong_shrunk_up),
            "significant_down": len(strong_shrunk_down),
        },
    ]

    write_tsv(
        TABLES_OUT / "fusarium_cazyme_summary.tsv",
        [
            "category",
            "total",
            "significant_up",
            "significant_down",
        ],
        summary,
    )

    family_rows = make_family_summary(
        high_confidence,
        significant_high_confidence,
        strong_raw_high_confidence,
        strong_shrunk_high_confidence,
    )

    write_tsv(
        TABLES_OUT / "fusarium_cazyme_family_summary.tsv",
        [
            "cazyme_family",
            "cazyme_class",
            "high_confidence_genes",
            "significant_genes",
            "significant_up",
            "significant_down",
            "strong_raw_lfc_genes",
            "strong_shrunk_lfc_genes",
        ],
        family_rows,
    )

    class_rows = make_class_summary(family_rows)

    write_tsv(
        TABLES_OUT / "fusarium_cazyme_class_summary.tsv",
        [
            "cazyme_class",
            "families",
            "high_confidence_genes",
            "significant_genes",
            "significant_up",
            "significant_down",
            "strong_raw_lfc_genes",
            "strong_shrunk_lfc_genes",
        ],
        class_rows,
    )

    substrate_rows = make_substrate_summary(high_confidence)

    write_tsv(
        TABLES_OUT / "fusarium_cazyme_substrate_summary.tsv",
        [
            "substrate",
            "high_confidence_genes",
            "significant_genes",
            "significant_up",
            "significant_down",
        ],
        substrate_rows,
    )

    validation = [
        validation_row(
            "total_gene_rows",
            len(rows),
            EXPECTED_TOTAL_GENES,
            len(rows) == EXPECTED_TOTAL_GENES,
        ),
        validation_row(
            "unique_gene_ids",
            len(seen),
            EXPECTED_TOTAL_GENES,
            len(seen) == EXPECTED_TOTAL_GENES,
        ),
        validation_row(
            "genes_in_deseq2_dataset",
            len(deseq2_rows),
            EXPECTED_DESEQ2_GENES,
            len(deseq2_rows) == EXPECTED_DESEQ2_GENES,
        ),
        validation_row(
            "genes_with_padj",
            len(padj_rows),
            EXPECTED_PADJ_GENES,
            len(padj_rows) == EXPECTED_PADJ_GENES,
        ),
        validation_row(
            "significant_genes",
            len(all_significant),
            EXPECTED_SIGNIFICANT_GENES,
            len(all_significant) == EXPECTED_SIGNIFICANT_GENES,
        ),
        validation_row(
            "dbcan_any_hit_gene_flag_matches_counts",
            any_hit_count_consistency,
            EXPECTED_TOTAL_GENES,
            any_hit_count_consistency == EXPECTED_TOTAL_GENES,
        ),
        validation_row(
            "dbcan_high_confidence_gene_flag_matches_counts",
            high_confidence_count_consistency,
            EXPECTED_TOTAL_GENES,
            high_confidence_count_consistency
            == EXPECTED_TOTAL_GENES,
        ),
        validation_row(
            "high_confidence_subset_of_any_hit",
            high_confidence_subset_consistency,
            EXPECTED_TOTAL_GENES,
            high_confidence_subset_consistency
            == EXPECTED_TOTAL_GENES,
        ),
        validation_row(
            "dbcan_any_hit_genes",
            len(any_hit),
            EXPECTED_ANY_HIT_GENES,
            len(any_hit) == EXPECTED_ANY_HIT_GENES,
        ),
        validation_row(
            "high_confidence_cazyme_genes",
            len(high_confidence),
            EXPECTED_HIGH_CONFIDENCE_GENES,
            len(high_confidence)
            == EXPECTED_HIGH_CONFIDENCE_GENES,
        ),
        validation_row(
            "high_confidence_cazymes_in_deseq2",
            len(high_confidence_deseq2),
            EXPECTED_HIGH_CONFIDENCE_DESEQ2,
            len(high_confidence_deseq2)
            == EXPECTED_HIGH_CONFIDENCE_DESEQ2,
        ),
        validation_row(
            "high_confidence_cazymes_with_padj",
            len(high_confidence_padj),
            EXPECTED_HIGH_CONFIDENCE_PADJ,
            len(high_confidence_padj)
            == EXPECTED_HIGH_CONFIDENCE_PADJ,
        ),
        validation_row(
            "significant_high_confidence_cazymes",
            len(significant_high_confidence),
            EXPECTED_SIGNIFICANT_HIGH_CONFIDENCE,
            len(significant_high_confidence)
            == EXPECTED_SIGNIFICANT_HIGH_CONFIDENCE,
        ),
        validation_row(
            "significant_high_confidence_direction_sum",
            len(significant_up) + len(significant_down),
            len(significant_high_confidence),
            len(significant_up) + len(significant_down)
            == len(significant_high_confidence),
        ),
        validation_row(
            "significant_high_confidence_up",
            len(significant_up),
            EXPECTED_SIGNIFICANT_UP,
            len(significant_up) == EXPECTED_SIGNIFICANT_UP,
        ),
        validation_row(
            "significant_high_confidence_down",
            len(significant_down),
            EXPECTED_SIGNIFICANT_DOWN,
            len(significant_down) == EXPECTED_SIGNIFICANT_DOWN,
        ),
        validation_row(
            "strong_raw_high_confidence_cazymes",
            len(strong_raw_high_confidence),
            EXPECTED_STRONG_RAW_HIGH_CONFIDENCE,
            len(strong_raw_high_confidence)
            == EXPECTED_STRONG_RAW_HIGH_CONFIDENCE,
        ),
        validation_row(
            "strong_raw_high_confidence_up",
            len(strong_raw_up),
            EXPECTED_STRONG_RAW_UP,
            len(strong_raw_up) == EXPECTED_STRONG_RAW_UP,
        ),
        validation_row(
            "strong_raw_high_confidence_down",
            len(strong_raw_down),
            EXPECTED_STRONG_RAW_DOWN,
            len(strong_raw_down) == EXPECTED_STRONG_RAW_DOWN,
        ),
        validation_row(
            "strong_shrunk_high_confidence_cazymes",
            len(strong_shrunk_high_confidence),
            EXPECTED_STRONG_SHRUNK_HIGH_CONFIDENCE,
            len(strong_shrunk_high_confidence)
            == EXPECTED_STRONG_SHRUNK_HIGH_CONFIDENCE,
        ),
        validation_row(
            "strong_shrunk_high_confidence_up",
            len(strong_shrunk_up),
            EXPECTED_STRONG_SHRUNK_UP,
            len(strong_shrunk_up)
            == EXPECTED_STRONG_SHRUNK_UP,
        ),
        validation_row(
            "strong_shrunk_high_confidence_down",
            len(strong_shrunk_down),
            EXPECTED_STRONG_SHRUNK_DOWN,
            len(strong_shrunk_down)
            == EXPECTED_STRONG_SHRUNK_DOWN,
        ),
        validation_row(
            "strong_shrunk_subset_of_significant",
            len(strong_shrunk_high_confidence),
            f"<= {len(significant_high_confidence)}",
            set(
                row["gene_id"]
                for row in strong_shrunk_high_confidence
            ).issubset(
                row["gene_id"]
                for row in significant_high_confidence
            ),
        ),
    ]

    if any(row["status"] == "FAIL" for row in validation):
        failed = [
            str(row["check"])
            for row in validation
            if row["status"] == "FAIL"
        ]

        fail(
            "one or more validation checks failed: "
            + ", ".join(failed)
        )

    write_tsv(
        DIAGNOSTICS_OUT / "validation.tsv",
        ["check", "value", "expected", "status"],
        validation,
    )

    not_tested_high_confidence = [
        row
        for row in high_confidence
        if row["in_deseq2_dataset"] == "FALSE"
    ]

    write_tsv(
        DIAGNOSTICS_OUT
        / "high_confidence_cazymes_excluded_from_deseq2.tsv",
        fields,
        not_tested_high_confidence,
    )

    run_info = [
        {
            "field": "script",
            "value": "24_run_fusarium_cazyme_analysis.py",
        },
        {
            "field": "run_timestamp_utc",
            "value": datetime.now(timezone.utc).isoformat(),
        },
        {
            "field": "repository_root",
            "value": str(ROOT),
        },
        {
            "field": "input_file",
            "value": str(INPUT),
        },
        {
            "field": "input_sha256",
            "value": sha256(INPUT),
        },
        {
            "field": "primary_cazyme_definition",
            "value": (
                "dbcan_high_confidence == TRUE; one or more "
                "gene-associated proteins supported by at least "
                "two dbCAN methods"
            ),
        },
        {
            "field": "secondary_cazyme_definition",
            "value": (
                "dbcan_any_hit == TRUE; one or more gene-associated "
                "proteins supported by at least one dbCAN method"
            ),
        },
        {
            "field": "significant_definition",
            "value": "raw_padj < 0.05",
        },
        {
            "field": "strong_raw_definition",
            "value": (
                "raw_padj < 0.05 and "
                "abs(raw_log2FoldChange) > 1"
            ),
        },
        {
            "field": "strong_shrunk_definition",
            "value": (
                "raw_padj < 0.05 and "
                "abs(shrunk_log2FoldChange) > 1"
            ),
        },
        {
            "field": "family_summary_source",
            "value": (
                "dbcan_recommended; semicolon-delimited assignments "
                "parsed and deduplicated within each Trinity gene"
            ),
        },
        {
            "field": "statistical_testing",
            "value": (
                "No new statistical test performed; tables are "
                "descriptive subsets of validated DESeq2 results"
            ),
        },
    ]

    write_tsv(
        OUTDIR / "run_info.tsv",
        ["field", "value"],
        run_info,
    )

    checksum_paths = [
        INPUT,
        DIAGNOSTICS_OUT / "validation.tsv",
        TABLES_OUT / "fusarium_cazyme_summary.tsv",
        TABLES_OUT / "fusarium_cazyme_family_summary.tsv",
        TABLES_OUT / "fusarium_cazyme_class_summary.tsv",
        TABLES_OUT / "fusarium_cazyme_substrate_summary.tsv",
    ]

    checksum_rows = [
        {
            "sha256": sha256(path),
            "file": str(path.relative_to(ROOT)),
        }
        for path in checksum_paths
    ]

    write_tsv(
        OUTDIR / "checksums.sha256",
        ["sha256", "file"],
        checksum_rows,
    )

    print("Fusarium CAZyme analysis completed successfully.")
    print(f"Total annotation genes:              {len(rows):,}")
    print(f"dbCAN any-hit genes:                 {len(any_hit):,}")
    print(
        "High-confidence CAZyme genes:        "
        f"{len(high_confidence):,}"
    )
    print(
        "High-confidence genes in DESeq2:     "
        f"{len(high_confidence_deseq2):,}"
    )
    print(
        "High-confidence genes with padj:     "
        f"{len(high_confidence_padj):,}"
    )
    print(
        "Significant high-confidence genes:   "
        f"{len(significant_high_confidence):,} "
        f"(up={len(significant_up):,}, "
        f"down={len(significant_down):,})"
    )
    print(
        "Strong high-confidence, raw LFC:     "
        f"{len(strong_raw_high_confidence):,} "
        f"(up={len(strong_raw_up):,}, "
        f"down={len(strong_raw_down):,})"
    )
    print(
        "Strong high-confidence, shrunk LFC:  "
        f"{len(strong_shrunk_high_confidence):,} "
        f"(up={len(strong_shrunk_up):,}, "
        f"down={len(strong_shrunk_down):,})"
    )
    print(
        "Distinct recommended CAZyme families: "
        f"{len(family_rows):,}"
    )
    print(
        "Validation:                          "
        f"{DIAGNOSTICS_OUT / 'validation.tsv'}"
    )
    print(
        "Summary:                             "
        f"{TABLES_OUT / 'fusarium_cazyme_summary.tsv'}"
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
