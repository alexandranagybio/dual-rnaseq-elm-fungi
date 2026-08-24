#!/usr/bin/env python3
"""
15_cazyme_analysis.py

Integrate validated gene-level DESeq2 results with dbCAN 5.2.8 annotations
for Ophiostoma novo-ulmi.

Primary CAZyme definition
-------------------------
High-confidence CAZyme:
    dbcan_high_confidence == TRUE
    equivalent to dbCAN #ofTools >= 2

Secondary annotation retained
-----------------------------
Any dbCAN hit:
    dbcan_any_hit == TRUE
    equivalent to dbCAN #ofTools >= 1

Differential-expression definitions
-----------------------------------
Significant:
    padj < 0.05

Raw-LFC sensitivity subset:
    padj < 0.05 and abs(log2FoldChange) > 1

Shrunken-LFC sensitivity subset:
    padj < 0.05 and abs(shrunk_log2FoldChange) > 1

The raw log2 fold-change sign determines direction. Shrunken fold changes are
used for stabilized effect-size filtering and ranking.

The script performs no new statistical test. It merges validated annotations
with the canonical DESeq2 tables containing apeglm-shrunken fold changes.

Run from the repository root:
    python workflow/scripts/15_cazyme_analysis.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

ANNOTATION = (
    ROOT
    / "results/ophiostoma/functional_annotation/"
      "ophiostoma_functional_annotation.tsv"
)

DESEQ2_TABLES = ROOT / "results/ophiostoma/deseq2_results/tables"

OUTDIR = ROOT / "results/ophiostoma/cazyme_analysis"
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
EXPECTED_ANY_DBCAN = 1440
EXPECTED_HIGH_CONFIDENCE_DBCAN = 315

REQUIRED_DE_COLUMNS = {
    "gene_id",
    "contrast",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "shrunk_log2FoldChange",
    "shrunk_lfcSE",
}

REQUIRED_ANNOTATION_COLUMNS = {
    "gene_id",
    "mrna_id",
    "dbcan_any_hit",
    "dbcan_high_confidence",
    "dbcan_n_tools",
    "dbcan_ec",
    "dbcan_hmm",
    "dbcan_sub",
    "dbcan_diamond",
    "dbcan_recommended",
    "dbcan_substrate",
}

DEPRECATED_OUTPUT_SUFFIXES = (
    "_high_confidence_cazyme_strong.tsv",
    "_high_confidence_cazyme_strong_up.tsv",
    "_high_confidence_cazyme_strong_down.tsv",
)


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
            f"ERROR: invalid {field} for {gene_id}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return value != value


def as_bool_string(value: bool) -> str:
    return str(value).upper()


def read_annotation() -> tuple[list[str], dict[str, dict[str, str]]]:
    fields, rows = read_tsv(ANNOTATION)

    missing = REQUIRED_ANNOTATION_COLUMNS - set(fields)

    if missing:
        fail(
            "annotation table is missing columns: "
            f"{sorted(missing)}"
        )

    annotation: dict[str, dict[str, str]] = {}

    for row in rows:
        gene_id = row["gene_id"]

        if gene_id in annotation:
            fail(f"duplicated annotation gene ID: {gene_id}")

        if row["dbcan_any_hit"] not in {"TRUE", "FALSE"}:
            fail(
                f"invalid dbcan_any_hit value for {gene_id}: "
                f"{row['dbcan_any_hit']}"
            )

        if row["dbcan_high_confidence"] not in {"TRUE", "FALSE"}:
            fail(
                f"invalid dbcan_high_confidence value for {gene_id}: "
                f"{row['dbcan_high_confidence']}"
            )

        try:
            n_tools = int(row["dbcan_n_tools"])
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: invalid dbcan_n_tools for {gene_id}: "
                f"{row['dbcan_n_tools']}"
            ) from exc

        expected_any = n_tools >= 1
        expected_high_confidence = n_tools >= 2

        if (row["dbcan_any_hit"] == "TRUE") != expected_any:
            fail(
                f"inconsistent dbCAN any-hit flag for {gene_id}: "
                f"n_tools={n_tools}, "
                f"dbcan_any_hit={row['dbcan_any_hit']}"
            )

        if (
            row["dbcan_high_confidence"] == "TRUE"
        ) != expected_high_confidence:
            fail(
                f"inconsistent high-confidence flag for {gene_id}: "
                f"n_tools={n_tools}, "
                "dbcan_high_confidence="
                f"{row['dbcan_high_confidence']}"
            )

        annotation[gene_id] = row

    if len(annotation) != EXPECTED_ANNOTATION_GENES:
        fail(
            f"expected {EXPECTED_ANNOTATION_GENES} annotation genes, "
            f"observed {len(annotation)}"
        )

    return fields, annotation


def sort_key(
    row: dict[str, object],
) -> tuple[float, float, float, str]:
    gene_id = str(row["gene_id"])

    padj = parse_float(
        str(row["padj"]),
        "padj",
        gene_id,
    )

    shrunk_lfc = parse_float(
        str(row["shrunk_log2FoldChange"]),
        "shrunk_log2FoldChange",
        gene_id,
    )

    raw_lfc = parse_float(
        str(row["log2FoldChange"]),
        "log2FoldChange",
        gene_id,
    )

    padj_key = padj if not is_nan(padj) else float("inf")
    shrunk_key = (
        -abs(shrunk_lfc)
        if not is_nan(shrunk_lfc)
        else 0.0
    )
    raw_key = (
        -abs(raw_lfc)
        if not is_nan(raw_lfc)
        else 0.0
    )

    return padj_key, shrunk_key, raw_key, gene_id


def remove_deprecated_outputs() -> None:
    for contrast in CONTRASTS:
        for suffix in DEPRECATED_OUTPUT_SUFFIXES:
            path = TABLES_OUT / f"{contrast}{suffix}"

            if path.exists():
                path.unlink()


def analyse_contrast(
    contrast: str,
    annotation_fields: list[str],
    annotation: dict[str, dict[str, str]],
) -> dict[str, object]:
    input_path = DESEQ2_TABLES / f"{contrast}_lfcshrunk.tsv"

    de_fields, de_rows = read_tsv(input_path)

    missing = REQUIRED_DE_COLUMNS - set(de_fields)

    if missing:
        fail(
            f"{input_path} is missing columns: "
            f"{sorted(missing)}"
        )

    seen: set[str] = set()
    merged_rows: list[dict[str, object]] = []

    for de_row in de_rows:
        gene_id = de_row["gene_id"]

        if gene_id in seen:
            fail(
                f"duplicated DESeq2 gene ID in {contrast}: "
                f"{gene_id}"
            )

        seen.add(gene_id)

        if gene_id not in annotation:
            fail(
                "DESeq2 gene absent from functional annotation in "
                f"{contrast}: {gene_id}"
            )

        if de_row["contrast"] != contrast:
            fail(
                f"unexpected contrast for {gene_id}: "
                f"{de_row['contrast']}; expected {contrast}"
            )

        padj = parse_float(
            de_row["padj"],
            "padj",
            gene_id,
        )

        raw_log2fc = parse_float(
            de_row["log2FoldChange"],
            "log2FoldChange",
            gene_id,
        )

        shrunk_log2fc = parse_float(
            de_row["shrunk_log2FoldChange"],
            "shrunk_log2FoldChange",
            gene_id,
        )

        if is_nan(shrunk_log2fc):
            fail(
                f"missing shrunken log2 fold-change for {gene_id} "
                f"in {contrast}"
            )

        significant = (
            not is_nan(padj)
            and padj < ALPHA
        )

        significant_raw_lfc1 = (
            significant
            and not is_nan(raw_log2fc)
            and abs(raw_log2fc) > STRONG_ABS_LOG2FC
        )

        significant_shrunk_lfc1 = (
            significant
            and abs(shrunk_log2fc) > STRONG_ABS_LOG2FC
        )

        if is_nan(raw_log2fc):
            direction = "NA"
        elif raw_log2fc > 0:
            direction = "up"
        elif raw_log2fc < 0:
            direction = "down"
        else:
            direction = "unchanged"

        any_dbcan = (
            annotation[gene_id]["dbcan_any_hit"] == "TRUE"
        )

        high_confidence = (
            annotation[gene_id]["dbcan_high_confidence"]
            == "TRUE"
        )

        merged: dict[str, object] = dict(de_row)

        for field in annotation_fields:
            if field != "gene_id":
                merged[field] = annotation[gene_id][field]

        merged["de_significant_padj_lt_0.05"] = (
            as_bool_string(significant)
        )

        merged[
            "de_significant_raw_lfc1"
        ] = as_bool_string(significant_raw_lfc1)

        merged[
            "de_significant_shrunk_lfc1"
        ] = as_bool_string(significant_shrunk_lfc1)

        merged["de_direction"] = direction

        merged["dbcan_any_hit_candidate"] = (
            as_bool_string(any_dbcan)
        )

        merged["high_confidence_cazyme"] = (
            as_bool_string(high_confidence)
        )

        merged[
            "significant_high_confidence_cazyme"
        ] = as_bool_string(
            significant and high_confidence
        )

        merged[
            "significant_raw_lfc1_high_confidence_cazyme"
        ] = as_bool_string(
            significant_raw_lfc1 and high_confidence
        )

        merged[
            "significant_shrunk_lfc1_high_confidence_cazyme"
        ] = as_bool_string(
            significant_shrunk_lfc1 and high_confidence
        )

        merged_rows.append(merged)

    if len(seen) != EXPECTED_TESTED_GENES:
        fail(
            f"{contrast}: expected "
            f"{EXPECTED_TESTED_GENES} tested genes, "
            f"observed {len(seen)}"
        )

    output_fields = (
        de_fields
        + [
            field
            for field in annotation_fields
            if field != "gene_id"
        ]
        + [
            "de_significant_padj_lt_0.05",
            "de_significant_raw_lfc1",
            "de_significant_shrunk_lfc1",
            "de_direction",
            "dbcan_any_hit_candidate",
            "high_confidence_cazyme",
            "significant_high_confidence_cazyme",
            "significant_raw_lfc1_high_confidence_cazyme",
            "significant_shrunk_lfc1_high_confidence_cazyme",
        ]
    )

    any_hit_rows = [
        row
        for row in merged_rows
        if row["dbcan_any_hit_candidate"] == "TRUE"
    ]

    high_confidence_rows = [
        row
        for row in merged_rows
        if row["high_confidence_cazyme"] == "TRUE"
    ]

    significant_rows = [
        row
        for row in high_confidence_rows
        if row[
            "significant_high_confidence_cazyme"
        ] == "TRUE"
    ]

    significant_up = [
        row
        for row in significant_rows
        if row["de_direction"] == "up"
    ]

    significant_down = [
        row
        for row in significant_rows
        if row["de_direction"] == "down"
    ]

    raw_lfc1_rows = [
        row
        for row in high_confidence_rows
        if row[
            "significant_raw_lfc1_high_confidence_cazyme"
        ] == "TRUE"
    ]

    raw_lfc1_up = [
        row
        for row in raw_lfc1_rows
        if row["de_direction"] == "up"
    ]

    raw_lfc1_down = [
        row
        for row in raw_lfc1_rows
        if row["de_direction"] == "down"
    ]

    shrunk_lfc1_rows = [
        row
        for row in high_confidence_rows
        if row[
            "significant_shrunk_lfc1_high_confidence_cazyme"
        ] == "TRUE"
    ]

    shrunk_lfc1_up = [
        row
        for row in shrunk_lfc1_rows
        if row["de_direction"] == "up"
    ]

    shrunk_lfc1_down = [
        row
        for row in shrunk_lfc1_rows
        if row["de_direction"] == "down"
    ]

    subsets = (
        merged_rows,
        any_hit_rows,
        high_confidence_rows,
        significant_rows,
        significant_up,
        significant_down,
        raw_lfc1_rows,
        raw_lfc1_up,
        raw_lfc1_down,
        shrunk_lfc1_rows,
        shrunk_lfc1_up,
        shrunk_lfc1_down,
    )

    for subset in subsets:
        subset.sort(key=sort_key)

    outputs = {
        f"{contrast}_all_genes_with_dbcan.tsv":
            merged_rows,

        f"{contrast}_dbcan_any_hit_all.tsv":
            any_hit_rows,

        f"{contrast}_high_confidence_cazyme_all.tsv":
            high_confidence_rows,

        f"{contrast}_high_confidence_cazyme_significant.tsv":
            significant_rows,

        f"{contrast}_high_confidence_cazyme_significant_up.tsv":
            significant_up,

        f"{contrast}_high_confidence_cazyme_significant_down.tsv":
            significant_down,

        f"{contrast}_high_confidence_cazyme_"
        "significant_raw_lfc1.tsv":
            raw_lfc1_rows,

        f"{contrast}_high_confidence_cazyme_"
        "significant_raw_lfc1_up.tsv":
            raw_lfc1_up,

        f"{contrast}_high_confidence_cazyme_"
        "significant_raw_lfc1_down.tsv":
            raw_lfc1_down,

        f"{contrast}_high_confidence_cazyme_"
        "significant_shrunk_lfc1.tsv":
            shrunk_lfc1_rows,

        f"{contrast}_high_confidence_cazyme_"
        "significant_shrunk_lfc1_up.tsv":
            shrunk_lfc1_up,

        f"{contrast}_high_confidence_cazyme_"
        "significant_shrunk_lfc1_down.tsv":
            shrunk_lfc1_down,
    }

    for filename, rows in outputs.items():
        write_tsv(
            TABLES_OUT / filename,
            output_fields,
            rows,
        )

    significant_total = sum(
        row["de_significant_padj_lt_0.05"] == "TRUE"
        for row in merged_rows
    )

    raw_lfc1_total = sum(
        row["de_significant_raw_lfc1"] == "TRUE"
        for row in merged_rows
    )

    shrunk_lfc1_total = sum(
        row["de_significant_shrunk_lfc1"] == "TRUE"
        for row in merged_rows
    )

    tested_high_confidence = len(high_confidence_rows)
    tested_any_hit = len(any_hit_rows)

    return {
        "contrast": contrast,
        "tested_genes": len(merged_rows),
        "dbcan_any_hit_tested_genes": tested_any_hit,
        "high_confidence_cazyme_tested_genes":
            tested_high_confidence,
        "non_high_confidence_tested_genes": (
            len(merged_rows) - tested_high_confidence
        ),
        "significant_de_genes_padj_lt_0.05":
            significant_total,
        "significant_high_confidence_cazymes":
            len(significant_rows),
        "significant_high_confidence_cazymes_up":
            len(significant_up),
        "significant_high_confidence_cazymes_down":
            len(significant_down),
        "raw_lfc1_de_genes":
            raw_lfc1_total,
        "raw_lfc1_high_confidence_cazymes":
            len(raw_lfc1_rows),
        "raw_lfc1_high_confidence_cazymes_up":
            len(raw_lfc1_up),
        "raw_lfc1_high_confidence_cazymes_down":
            len(raw_lfc1_down),
        "shrunk_lfc1_de_genes":
            shrunk_lfc1_total,
        "shrunk_lfc1_high_confidence_cazymes":
            len(shrunk_lfc1_rows),
        "shrunk_lfc1_high_confidence_cazymes_up":
            len(shrunk_lfc1_up),
        "shrunk_lfc1_high_confidence_cazymes_down":
            len(shrunk_lfc1_down),
        "high_confidence_cazyme_fraction_of_tested": (
            tested_high_confidence / len(merged_rows)
        ),
        "high_confidence_cazyme_fraction_of_significant": (
            len(significant_rows) / significant_total
            if significant_total
            else 0.0
        ),
        "high_confidence_cazyme_fraction_of_raw_lfc1": (
            len(raw_lfc1_rows) / raw_lfc1_total
            if raw_lfc1_total
            else 0.0
        ),
        "high_confidence_cazyme_fraction_of_shrunk_lfc1": (
            len(shrunk_lfc1_rows) / shrunk_lfc1_total
            if shrunk_lfc1_total
            else 0.0
        ),
        "input_file": str(input_path),
    }


def main() -> None:
    require_file(ANNOTATION)

    TABLES_OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    DIAGNOSTICS_OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    remove_deprecated_outputs()

    annotation_fields, annotation = read_annotation()

    complete_any = sum(
        row["dbcan_any_hit"] == "TRUE"
        for row in annotation.values()
    )

    complete_high_confidence = sum(
        row["dbcan_high_confidence"] == "TRUE"
        for row in annotation.values()
    )

    if complete_any != EXPECTED_ANY_DBCAN:
        fail(
            f"expected {EXPECTED_ANY_DBCAN} genes with any "
            f"dbCAN hit, observed {complete_any}"
        )

    if (
        complete_high_confidence
        != EXPECTED_HIGH_CONFIDENCE_DBCAN
    ):
        fail(
            f"expected {EXPECTED_HIGH_CONFIDENCE_DBCAN} "
            "high-confidence CAZymes, observed "
            f"{complete_high_confidence}"
        )

    summaries = [
        analyse_contrast(
            contrast,
            annotation_fields,
            annotation,
        )
        for contrast in CONTRASTS
    ]

    write_tsv(
        TABLES_OUT / "cazyme_summary.tsv",
        list(summaries[0]),
        summaries,
    )

    tested_sets: dict[str, set[str]] = {}

    for contrast in CONTRASTS:
        _, rows = read_tsv(
            DESEQ2_TABLES
            / f"{contrast}_lfcshrunk.tsv"
        )

        tested_sets[contrast] = {
            row["gene_id"]
            for row in rows
        }

    shared_tested = set.intersection(
        *tested_sets.values()
    )

    union_tested = set.union(
        *tested_sets.values()
    )

    excluded_from_testing = (
        set(annotation) - union_tested
    )

    excluded_high_confidence = sum(
        annotation[gene]["dbcan_high_confidence"]
        == "TRUE"
        for gene in excluded_from_testing
    )

    excluded_any_hit = sum(
        annotation[gene]["dbcan_any_hit"] == "TRUE"
        for gene in excluded_from_testing
    )

    diagnostics = [
        {
            "check": "annotation_gene_count",
            "value": len(annotation),
            "expected": EXPECTED_ANNOTATION_GENES,
            "status": "PASS",
        },
        {
            "check": "complete_annotation_dbcan_any_hit",
            "value": complete_any,
            "expected": EXPECTED_ANY_DBCAN,
            "status": "PASS",
        },
        {
            "check":
                "complete_annotation_high_confidence_cazymes",
            "value": complete_high_confidence,
            "expected": EXPECTED_HIGH_CONFIDENCE_DBCAN,
            "status": "PASS",
        },
        {
            "check": "tested_gene_count_per_contrast",
            "value": ",".join(
                str(len(tested_sets[contrast]))
                for contrast in CONTRASTS
            ),
            "expected": EXPECTED_TESTED_GENES,
            "status": "PASS",
        },
        {
            "check":
                "tested_gene_sets_identical_across_contrasts",
            "value": len(shared_tested),
            "expected": len(union_tested),
            "status": (
                "PASS"
                if len(shared_tested) == len(union_tested)
                else "FAIL"
            ),
        },
        {
            "check":
                "genes_excluded_from_deseq2_testing",
            "value": len(excluded_from_testing),
            "expected": (
                EXPECTED_ANNOTATION_GENES
                - EXPECTED_TESTED_GENES
            ),
            "status": (
                "PASS"
                if len(excluded_from_testing)
                == (
                    EXPECTED_ANNOTATION_GENES
                    - EXPECTED_TESTED_GENES
                )
                else "FAIL"
            ),
        },
        {
            "check":
                "excluded_genes_with_any_dbcan_hit",
            "value": excluded_any_hit,
            "expected": "descriptive",
            "status": "PASS",
        },
        {
            "check":
                "excluded_high_confidence_cazymes",
            "value": excluded_high_confidence,
            "expected": "descriptive",
            "status": "PASS",
        },
    ]

    if any(
        row["status"] == "FAIL"
        for row in diagnostics
    ):
        fail(
            "one or more cross-contrast validation "
            "checks failed"
        )

    write_tsv(
        DIAGNOSTICS_OUT / "validation.tsv",
        ["check", "value", "expected", "status"],
        diagnostics,
    )

    excluded_rows = [
        {
            "gene_id": gene,
            "mrna_id": annotation[gene]["mrna_id"],
            "dbcan_any_hit":
                annotation[gene]["dbcan_any_hit"],
            "dbcan_high_confidence":
                annotation[gene][
                    "dbcan_high_confidence"
                ],
            "dbcan_n_tools":
                annotation[gene]["dbcan_n_tools"],
            "dbcan_recommended":
                annotation[gene]["dbcan_recommended"],
        }
        for gene in sorted(excluded_from_testing)
    ]

    write_tsv(
        DIAGNOSTICS_OUT
        / "genes_excluded_from_deseq2_testing.tsv",
        [
            "gene_id",
            "mrna_id",
            "dbcan_any_hit",
            "dbcan_high_confidence",
            "dbcan_n_tools",
            "dbcan_recommended",
        ],
        excluded_rows,
    )

    run_info = [
        {
            "field": "script",
            "value": "15_cazyme_analysis.py",
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
            "field": "annotation_file",
            "value": str(ANNOTATION),
        },
        {
            "field": "annotation_sha256",
            "value": sha256(ANNOTATION),
        },
        {
            "field": "primary_cazyme_definition",
            "value": (
                "dbcan_high_confidence == TRUE; "
                "equivalent to dbCAN #ofTools >= 2"
            ),
        },
        {
            "field": "secondary_dbcan_definition",
            "value": (
                "dbcan_any_hit == TRUE; "
                "#ofTools >= 1"
            ),
        },
        {
            "field": "significant_de_definition",
            "value": "padj < 0.05",
        },
        {
            "field": "raw_lfc1_definition",
            "value": (
                "padj < 0.05 and "
                "abs(log2FoldChange) > 1"
            ),
        },
        {
            "field": "shrunk_lfc1_definition",
            "value": (
                "padj < 0.05 and "
                "abs(shrunk_log2FoldChange) > 1"
            ),
        },
        {
            "field": "direction_definition",
            "value": (
                "sign of raw log2FoldChange"
            ),
        },
        {
            "field": "effect_size_interpretation",
            "value": (
                "shrunk_log2FoldChange is used for "
                "stabilized effect-size filtering and ranking"
            ),
        },
        {
            "field": "interpretation",
            "value": (
                "High-confidence CAZyme candidates "
                "supported by at least two dbCAN methods; "
                "one-method hits retained in secondary tables"
            ),
        },
    ]

    for contrast in CONTRASTS:
        path = (
            DESEQ2_TABLES
            / f"{contrast}_lfcshrunk.tsv"
        )

        run_info.append(
            {
                "field":
                    f"{contrast}_input_sha256",
                "value": sha256(path),
            }
        )

    write_tsv(
        OUTDIR / "run_info.tsv",
        ["field", "value"],
        run_info,
    )

    print(
        "Ophiostoma high-confidence "
        "CAZyme analysis"
    )
    print()
    print(
        f"Complete annotation genes: "
        f"{len(annotation)}"
    )
    print(
        f"Genes with any dbCAN hit: "
        f"{complete_any}"
    )
    print(
        "Complete high-confidence CAZymes "
        "(#ofTools >= 2): "
        f"{complete_high_confidence}"
    )
    print(
        f"Genes tested per contrast: "
        f"{EXPECTED_TESTED_GENES}"
    )
    print()

    for row in summaries:
        print(row["contrast"])

        print(
            "  Significant high-confidence CAZymes: "
            f"{row['significant_high_confidence_cazymes']} "
            "(up "
            f"{row['significant_high_confidence_cazymes_up']}, "
            "down "
            f"{row['significant_high_confidence_cazymes_down']})"
        )

        print(
            "  Raw |LFC| > 1 high-confidence CAZymes: "
            f"{row['raw_lfc1_high_confidence_cazymes']} "
            "(up "
            f"{row['raw_lfc1_high_confidence_cazymes_up']}, "
            "down "
            f"{row['raw_lfc1_high_confidence_cazymes_down']})"
        )

        print(
            "  Shrunken |LFC| > 1 "
            "high-confidence CAZymes: "
            f"{row['shrunk_lfc1_high_confidence_cazymes']} "
            "(up "
            f"{row['shrunk_lfc1_high_confidence_cazymes_up']}, "
            "down "
            f"{row['shrunk_lfc1_high_confidence_cazymes_down']})"
        )

        print()

    print("Validation status: PASS")
    print(f"Output directory: {OUTDIR}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
