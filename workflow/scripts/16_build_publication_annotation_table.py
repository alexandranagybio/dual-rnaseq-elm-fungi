#!/usr/bin/env python3
"""
16_build_publication_annotation_table.py

Build validated Ophiostoma novo-ulmi publication-ready annotation tables by
combining:

1. Structural/JGI gene annotation, including aggregated GO and KEGG fields
2. SignalP 6.0 and dbCAN 5.2.8 functional annotations
3. Validated gene-level DESeq2 results for all three contrasts

Outputs
-------
results/ophiostoma/publication_annotation/
├── run_info.tsv
├── diagnostics/
│   ├── validation.tsv
│   └── genes_excluded_from_deseq2_testing.tsv
└── tables/
    ├── ophiostoma_publication_annotation_long.tsv
    ├── ophiostoma_publication_annotation_wide.tsv
    ├── interaction_vs_self_publication_annotation.tsv
    ├── interaction_vs_onu_publication_annotation.tsv
    └── onu_vs_self_publication_annotation.tsv

Definitions inherited from validated upstream scripts
------------------------------------------------------
Significant DE:
    padj < 0.05

Strong DE:
    padj < 0.05 and abs(log2FoldChange) > 1

SignalP-positive:
    signalp_is_sp == TRUE

High-confidence CAZyme:
    dbcan_high_confidence == TRUE
    equivalent to dbCAN #ofTools >= 2

The script performs no new statistical tests and does not reinterpret the
upstream annotations.

Run from repository root:
    python workflow/scripts/16_build_publication_annotation_table.py
"""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

STRUCTURAL_ANNOTATION = (
    ROOT
    / "results/ophiostoma/annotation/"
      "ophiostoma_gene_annotation_map.tsv"
)

FUNCTIONAL_ANNOTATION = (
    ROOT
    / "results/ophiostoma/functional_annotation/"
      "ophiostoma_functional_annotation.tsv"
)

DESEQ2_TABLES = ROOT / "results/ophiostoma/deseq2_results/tables"

OUTDIR = ROOT / "results/ophiostoma/publication_annotation"
TABLES_OUT = OUTDIR / "tables"
DIAGNOSTICS_OUT = OUTDIR / "diagnostics"

CONTRASTS = (
    "interaction_vs_self",
    "interaction_vs_onu",
    "onu_vs_self",
)

ALPHA = 0.05
STRONG_ABS_LOG2FC = 1.0

EXPECTED_STRUCTURAL_GENES = 8640
EXPECTED_FUNCTIONAL_GENES = 8640
EXPECTED_TESTED_GENES = 8560
EXPECTED_SIGNALP_POSITIVE = 591
EXPECTED_DBCAN_ANY_HIT = 1440
EXPECTED_HIGH_CONFIDENCE_CAZYMES = 315

EXPECTED_EGGNOG_ANNOTATED = 8231
EXPECTED_EGGNOG_UNANNOTATED = 409
EXPECTED_TESTED_EGGNOG_ANNOTATED = 8185
REQUIRED_STRUCTURAL_COLUMNS = {
    "gene_id",
    "protein_id",
    "transcript_id",
    "seqid",
    "start",
    "end",
    "strand",
    "go_count",
    "go_ids",
    "go_names",
    "go_types",
    "kegg_count",
    "kegg_ecnum",
    "kegg_definition",
    "kegg_pathway",
    "kegg_pathway_class",
    "kegg_pathway_type",
    "has_go",
    "has_kegg",
}

REQUIRED_FUNCTIONAL_COLUMNS = {
    "gene_id",
    "eggnog_annotated",
    "mrna_id",
    "signalp_prediction",
    "signalp_is_sp",
    "signalp_other_score",
    "signalp_sp_score",
    "signalp_cs_position",
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
    return math.isnan(value)


def read_unique_by_gene(
    path: Path,
    required_columns: set[str],
    label: str,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    fields, rows = read_tsv(path)
    missing = required_columns - set(fields)
    if missing:
        fail(f"{label} is missing columns: {sorted(missing)}")

    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        gene_id = row["gene_id"]
        if not gene_id:
            fail(f"blank gene_id in {label}")
        if gene_id in indexed:
            fail(f"duplicated gene_id in {label}: {gene_id}")
        indexed[gene_id] = row

    return fields, indexed


def de_flags(row: dict[str, str]) -> dict[str, str]:
    gene_id = row["gene_id"]
    padj = parse_float(row["padj"], "padj", gene_id)
    log2fc = parse_float(
        row["log2FoldChange"],
        "log2FoldChange",
        gene_id,
    )

    significant = not is_nan(padj) and padj < ALPHA
    strong = (
        significant
        and not is_nan(log2fc)
        and abs(log2fc) > STRONG_ABS_LOG2FC
    )

    if is_nan(log2fc):
        direction = "NA"
    elif log2fc > 0:
        direction = "up"
    elif log2fc < 0:
        direction = "down"
    else:
        direction = "unchanged"

    return {
        "de_significant_padj_lt_0.05": str(significant).upper(),
        "de_strong_padj_lt_0.05_abs_log2fc_gt_1": str(strong).upper(),
        "de_direction": direction,
    }


def main() -> None:
    for path in (STRUCTURAL_ANNOTATION, FUNCTIONAL_ANNOTATION):
        require_file(path)

    TABLES_OUT.mkdir(parents=True, exist_ok=True)
    DIAGNOSTICS_OUT.mkdir(parents=True, exist_ok=True)

    structural_fields, structural = read_unique_by_gene(
        STRUCTURAL_ANNOTATION,
        REQUIRED_STRUCTURAL_COLUMNS,
        "structural annotation",
    )

    functional_fields, functional = read_unique_by_gene(
        FUNCTIONAL_ANNOTATION,
        REQUIRED_FUNCTIONAL_COLUMNS,
        "functional annotation",
    )

    if len(structural) != EXPECTED_STRUCTURAL_GENES:
        fail(
            f"expected {EXPECTED_STRUCTURAL_GENES} structural genes, "
            f"observed {len(structural)}"
        )

    if len(functional) != EXPECTED_FUNCTIONAL_GENES:
        fail(
            f"expected {EXPECTED_FUNCTIONAL_GENES} functional genes, "
            f"observed {len(functional)}"
        )

    structural_ids = set(structural)
    functional_ids = set(functional)

    if structural_ids != functional_ids:
        only_structural = sorted(structural_ids - functional_ids)
        only_functional = sorted(functional_ids - structural_ids)
        fail(
            "structural and functional gene sets differ; "
            f"only structural={len(only_structural)}, "
            f"only functional={len(only_functional)}"
        )

    signalp_count = sum(
        row["signalp_is_sp"] == "TRUE"
        for row in functional.values()
    )
    dbcan_any_count = sum(
        row["dbcan_any_hit"] == "TRUE"
        for row in functional.values()
    )
    dbcan_high_count = sum(
        row["dbcan_high_confidence"] == "TRUE"
        for row in functional.values()
    )


    eggnog_annotated_count = sum(
        row["eggnog_annotated"] == "TRUE"
        for row in functional.values()
    )

    eggnog_unannotated_count = (
        len(functional) - eggnog_annotated_count
    )

    if eggnog_annotated_count != EXPECTED_EGGNOG_ANNOTATED:
        fail(
            f"expected {EXPECTED_EGGNOG_ANNOTATED} "
            "eggNOG-annotated genes, observed "
            f"{eggnog_annotated_count}"
        )

    if (
        eggnog_unannotated_count
        != EXPECTED_EGGNOG_UNANNOTATED
    ):
        fail(
            f"expected {EXPECTED_EGGNOG_UNANNOTATED} genes "
            "without final eggNOG annotation, observed "
            f"{eggnog_unannotated_count}"
        )

    if signalp_count != EXPECTED_SIGNALP_POSITIVE:
        fail(
            f"expected {EXPECTED_SIGNALP_POSITIVE} SignalP-positive genes, "
            f"observed {signalp_count}"
        )
    if dbcan_any_count != EXPECTED_DBCAN_ANY_HIT:
        fail(
            f"expected {EXPECTED_DBCAN_ANY_HIT} genes with any dbCAN hit, "
            f"observed {dbcan_any_count}"
        )
    if dbcan_high_count != EXPECTED_HIGH_CONFIDENCE_CAZYMES:
        fail(
            f"expected {EXPECTED_HIGH_CONFIDENCE_CAZYMES} "
            f"high-confidence CAZymes, observed {dbcan_high_count}"
        )

    de_by_contrast: dict[str, dict[str, dict[str, str]]] = {}
    de_fields_reference: list[str] | None = None

    for contrast in CONTRASTS:
        path = DESEQ2_TABLES / f"{contrast}_all_genes.tsv"
        de_fields, de_index = read_unique_by_gene(
            path,
            REQUIRED_DE_COLUMNS,
            f"DESeq2 table {contrast}",
        )

        if len(de_index) != EXPECTED_TESTED_GENES:
            fail(
                f"{contrast}: expected {EXPECTED_TESTED_GENES} tested genes, "
                f"observed {len(de_index)}"
            )

        if not set(de_index).issubset(structural_ids):
            missing = sorted(set(de_index) - structural_ids)
            fail(
                f"{contrast}: {len(missing)} DESeq2 genes are absent from "
                "the annotation universe"
            )

        for gene_id, row in de_index.items():
            if row["contrast"] != contrast:
                fail(
                    f"{contrast}: unexpected contrast value for "
                    f"{gene_id}: {row['contrast']}"
                )

        if de_fields_reference is None:
            de_fields_reference = de_fields
        elif de_fields != de_fields_reference:
            fail("DESeq2 column order differs across contrasts")

        de_by_contrast[contrast] = de_index

    tested_sets = {
        contrast: set(index)
        for contrast, index in de_by_contrast.items()
    }
    shared_tested = set.intersection(*tested_sets.values())
    union_tested = set.union(*tested_sets.values())

    if shared_tested != union_tested:
        fail("tested gene sets are not identical across contrasts")

    excluded_genes = structural_ids - union_tested

    if len(excluded_genes) != (
        EXPECTED_STRUCTURAL_GENES - EXPECTED_TESTED_GENES
    ):
        fail(
            f"expected {EXPECTED_STRUCTURAL_GENES - EXPECTED_TESTED_GENES} "
            f"genes excluded from DESeq2, observed {len(excluded_genes)}"
        )
    tested_eggnog_annotated_by_contrast = {
        contrast: sum(
            functional[gene_id]["eggnog_annotated"] == "TRUE"
            for gene_id in de_by_contrast[contrast]
        )
        for contrast in CONTRASTS
    }

    tested_eggnog_counts = set(
        tested_eggnog_annotated_by_contrast.values()
    )

    if len(tested_eggnog_counts) != 1:
        fail(
            "eggNOG-annotated tested-gene counts differ among "
            f"contrasts: {tested_eggnog_annotated_by_contrast}"
        )

    tested_eggnog_annotated_count = next(
        iter(tested_eggnog_counts)
    )

    if (
        tested_eggnog_annotated_count
        != EXPECTED_TESTED_EGGNOG_ANNOTATED
    ):
        fail(
            f"expected {EXPECTED_TESTED_EGGNOG_ANNOTATED} "
            "eggNOG-annotated tested genes per contrast, "
            f"observed {tested_eggnog_annotated_count}"
        )



    # Explicitly rename the functional mRNA identifier to avoid confusion
    # with the numeric JGI transcript_id in the structural map.
    functional_output_fields = [
        "signalp_mrna_id" if field == "mrna_id" else field
        for field in functional_fields
        if field != "gene_id"
    ]

    de_output_fields = [
        field
        for field in de_fields_reference
        if field not in {"gene_id"}
    ] + [
        "de_significant_padj_lt_0.05",
        "de_strong_padj_lt_0.05_abs_log2fc_gt_1",
        "de_direction",
    ]

    long_fields = (
        structural_fields
        + functional_output_fields
        + de_output_fields
    )

    long_rows: list[dict[str, object]] = []
    contrast_rows: dict[str, list[dict[str, object]]] = {
        contrast: [] for contrast in CONTRASTS
    }

    for contrast in CONTRASTS:
        for gene_id in sorted(
            de_by_contrast[contrast],
            key=lambda value: int(value.split("_", 1)[1]),
        ):
            row: dict[str, object] = dict(structural[gene_id])

            for field in functional_fields:
                if field == "gene_id":
                    continue
                output_field = (
                    "signalp_mrna_id" if field == "mrna_id" else field
                )
                row[output_field] = functional[gene_id][field]

            de_row = de_by_contrast[contrast][gene_id]
            for field in de_fields_reference:
                if field != "gene_id":
                    row[field] = de_row[field]

            row.update(de_flags(de_row))

            long_rows.append(row)
            contrast_rows[contrast].append(row)

    write_tsv(
        TABLES_OUT / "ophiostoma_publication_annotation_long.tsv",
        long_fields,
        long_rows,
    )

    for contrast in CONTRASTS:
        write_tsv(
            TABLES_OUT / f"{contrast}_publication_annotation.tsv",
            long_fields,
            contrast_rows[contrast],
        )

    # Wide table: all 8640 structural genes, including the 80 genes not tested.
    wide_de_fields: list[str] = []
    base_de_fields = [
        field
        for field in de_fields_reference
        if field not in {"gene_id", "contrast"}
    ] + [
        "de_significant_padj_lt_0.05",
        "de_strong_padj_lt_0.05_abs_log2fc_gt_1",
        "de_direction",
    ]

    for contrast in CONTRASTS:
        wide_de_fields.extend(
            f"{contrast}__{field}"
            for field in base_de_fields
        )

    wide_fields = (
        structural_fields
        + functional_output_fields
        + ["tested_by_deseq2"]
        + wide_de_fields
    )

    wide_rows: list[dict[str, object]] = []

    for gene_id in sorted(
        structural,
        key=lambda value: int(value.split("_", 1)[1]),
    ):
        row: dict[str, object] = dict(structural[gene_id])

        for field in functional_fields:
            if field == "gene_id":
                continue
            output_field = (
                "signalp_mrna_id" if field == "mrna_id" else field
            )
            row[output_field] = functional[gene_id][field]

        row["tested_by_deseq2"] = str(
            gene_id in union_tested
        ).upper()

        for contrast in CONTRASTS:
            prefix = f"{contrast}__"
            if gene_id in de_by_contrast[contrast]:
                de_row = de_by_contrast[contrast][gene_id]
                flags = de_flags(de_row)

                for field in base_de_fields:
                    if field in flags:
                        value = flags[field]
                    else:
                        value = de_row[field]
                    row[prefix + field] = value
            else:
                for field in base_de_fields:
                    row[prefix + field] = ""

        wide_rows.append(row)

    write_tsv(
        TABLES_OUT / "ophiostoma_publication_annotation_wide.tsv",
        wide_fields,
        wide_rows,
    )

    excluded_rows = []
    for gene_id in sorted(
        excluded_genes,
        key=lambda value: int(value.split("_", 1)[1]),
    ):
        excluded_rows.append(
            {
                "gene_id": gene_id,
                "protein_id": structural[gene_id]["protein_id"],
                "transcript_id": structural[gene_id]["transcript_id"],
                "signalp_is_sp": functional[gene_id]["signalp_is_sp"],
                "dbcan_any_hit": functional[gene_id]["dbcan_any_hit"],
                "dbcan_high_confidence": (
                    functional[gene_id]["dbcan_high_confidence"]
                ),
                "has_go": structural[gene_id]["has_go"],
                "has_kegg": structural[gene_id]["has_kegg"],
            }
        )

    write_tsv(
        DIAGNOSTICS_OUT / "genes_excluded_from_deseq2_testing.tsv",
        [
            "gene_id",
            "protein_id",
            "transcript_id",
            "signalp_is_sp",
            "dbcan_any_hit",
            "dbcan_high_confidence",
            "has_go",
            "has_kegg",
        ],
        excluded_rows,
    )

    go_annotated = sum(
        row["has_go"].lower() == "yes"
        for row in structural.values()
    )
    kegg_annotated = sum(
        row["has_kegg"].lower() == "yes"
        for row in structural.values()
    )

    validation = [
        {
            "check": "complete_eggnog_annotated",
            "value": eggnog_annotated_count,
            "expected": EXPECTED_EGGNOG_ANNOTATED,
            "status": "PASS",
        },
        {
            "check": "complete_eggnog_unannotated",
            "value": eggnog_unannotated_count,
            "expected": EXPECTED_EGGNOG_UNANNOTATED,
            "status": "PASS",
        },
        {
            "check": "tested_eggnog_annotated_per_contrast",
            "value": ",".join(
                str(
                    tested_eggnog_annotated_by_contrast[
                        contrast
                    ]
                )
                for contrast in CONTRASTS
            ),
            "expected": EXPECTED_TESTED_EGGNOG_ANNOTATED,
            "status": "PASS",
        },
        {
            "check": "structural_annotation_gene_count",
            "value": len(structural),
            "expected": EXPECTED_STRUCTURAL_GENES,
            "status": "PASS",
        },
        {
            "check": "functional_annotation_gene_count",
            "value": len(functional),
            "expected": EXPECTED_FUNCTIONAL_GENES,
            "status": "PASS",
        },
        {
            "check": "structural_and_functional_gene_sets_identical",
            "value": len(structural_ids & functional_ids),
            "expected": EXPECTED_STRUCTURAL_GENES,
            "status": "PASS",
        },
        {
            "check": "signalp_positive_complete",
            "value": signalp_count,
            "expected": EXPECTED_SIGNALP_POSITIVE,
            "status": "PASS",
        },
        {
            "check": "dbcan_any_hit_complete",
            "value": dbcan_any_count,
            "expected": EXPECTED_DBCAN_ANY_HIT,
            "status": "PASS",
        },
        {
            "check": "high_confidence_cazymes_complete",
            "value": dbcan_high_count,
            "expected": EXPECTED_HIGH_CONFIDENCE_CAZYMES,
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
            "check": "tested_gene_sets_identical",
            "value": len(shared_tested),
            "expected": EXPECTED_TESTED_GENES,
            "status": "PASS",
        },
        {
            "check": "genes_excluded_from_deseq2_testing",
            "value": len(excluded_genes),
            "expected": (
                EXPECTED_STRUCTURAL_GENES - EXPECTED_TESTED_GENES
            ),
            "status": "PASS",
        },
        {
            "check": "long_table_rows",
            "value": len(long_rows),
            "expected": EXPECTED_TESTED_GENES * len(CONTRASTS),
            "status": "PASS",
        },
        {
            "check": "wide_table_rows",
            "value": len(wide_rows),
            "expected": EXPECTED_STRUCTURAL_GENES,
            "status": "PASS",
        },
        {
            "check": "complete_genes_with_GO",
            "value": go_annotated,
            "expected": "descriptive",
            "status": "PASS",
        },
        {
            "check": "complete_genes_with_KEGG",
            "value": kegg_annotated,
            "expected": "descriptive",
            "status": "PASS",
        },
    ]

    write_tsv(
        DIAGNOSTICS_OUT / "validation.tsv",
        ["check", "value", "expected", "status"],
        validation,
    )

    run_info = [
        {
            "field": "script",
            "value": "16_build_publication_annotation_table.py",
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
            "field": "structural_annotation_file",
            "value": str(STRUCTURAL_ANNOTATION),
        },
        {
            "field": "structural_annotation_sha256",
            "value": sha256(STRUCTURAL_ANNOTATION),
        },
        {
            "field": "functional_annotation_file",
            "value": str(FUNCTIONAL_ANNOTATION),
        },
        {
            "field": "functional_annotation_sha256",
            "value": sha256(FUNCTIONAL_ANNOTATION),
        },
        {
            "field": "significant_de_definition",
            "value": "padj < 0.05",
        },
        {
            "field": "strong_de_definition",
            "value": (
                "padj < 0.05 and abs(log2FoldChange) > 1"
            ),
        },
        {
            "field": "signalp_positive_definition",
            "value": "signalp_is_sp == TRUE",
        },
        {
            "field": "high_confidence_cazyme_definition",
            "value": (
                "dbcan_high_confidence == TRUE; dbCAN #ofTools >= 2"
            ),
        },
        {
            "field": "long_table_scope",
            "value": (
                "8560 DESeq2-tested genes x 3 contrasts"
            ),
        },
        {
            "field": "wide_table_scope",
            "value": (
                "all 8640 structural genes; blank DE fields for "
                "80 genes excluded before DESeq2 inference"
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

    print("Ophiostoma publication annotation table build")
    print()
    print(f"Structural genes: {len(structural)}")
    print(f"Functional annotation genes: {len(functional)}")
    print(f"DESeq2-tested genes per contrast: {len(shared_tested)}")
    print(f"Genes excluded before DESeq2: {len(excluded_genes)}")
    print(f"Genes with GO annotation: {go_annotated}")
    print(f"Genes with KEGG annotation: {kegg_annotated}")
    print()
    print(f"Long table rows: {len(long_rows)}")
    print(f"Wide table rows: {len(wide_rows)}")
    print()
    print("Validation status: PASS")
    print(f"Output directory: {OUTDIR}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
