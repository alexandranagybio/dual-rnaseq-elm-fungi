#!/usr/bin/env python3
"""
Build Fusarium publication annotation tables by joining the validated
gene-level functional annotation table to raw and shrunken DESeq2 results.

The full table retains the complete Trinity gene universe. DESeq2 status is
represented explicitly so that genes absent from DESeq2, genes with NA padj,
and nonsignificant genes are never conflated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


EXPECTED = {
    "annotation_gene_rows": 15192,
    "raw_deseq2_rows": 10414,
    "shrunk_deseq2_rows": 10414,
    "genes_with_padj": 9000,
    "genes_with_na_padj": 1414,
    "significant_padj_lt_0.05": 2973,
    "significant_raw_abs_lfc_gt_1": 149,
    "significant_shrunk_abs_lfc_gt_1": 104,
    "significant_with_predicted_protein": 2834,
    "deseq2_ids_not_in_annotation": 0,
    "annotation_ids_not_in_deseq2": 4778,
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def parse_float(text: str) -> float | None:
    value = text.strip()
    if value == "" or value.upper() == "NA":
        return None
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"cannot parse numeric value: {text!r}") from exc


def read_tsv(path: Path) -> tuple[List[str], List[Dict[str, str]]]:
    if not path.is_file():
        fail(f"missing input file: {path}")

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"missing header: {path}")
        rows = list(reader)

    return list(reader.fieldnames), rows


def index_unique(
    rows: Sequence[Mapping[str, str]],
    key: str,
    label: str,
) -> Dict[str, Mapping[str, str]]:
    result: Dict[str, Mapping[str, str]] = {}
    duplicates: List[str] = []

    for row in rows:
        identifier = row.get(key, "").strip()
        if not identifier:
            fail(f"blank {key} in {label}")
        if identifier in result:
            duplicates.append(identifier)
        result[identifier] = row

    if duplicates:
        fail(
            f"duplicate {key} values in {label}: "
            + ", ".join(sorted(set(duplicates))[:10])
        )

    return result


def write_tsv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def natural_gene_key(gene_id: str) -> tuple:
    import re

    parts = re.split(r"(\d+)", gene_id)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Fusarium publication annotation tables."
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=Path(
            "results/fusarium/functional_annotation/tables/"
            "fusarium_gene_functional_annotation.tsv"
        ),
    )
    parser.add_argument(
        "--raw",
        type=Path,
        default=Path(
            "results/fusarium/deseq2_results/"
            "fusarium_interaction_vs_self_raw.tsv"
        ),
    )
    parser.add_argument(
        "--shrunk",
        type=Path,
        default=Path(
            "results/fusarium/deseq2_results/"
            "fusarium_interaction_vs_self_lfcshrunk.tsv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/fusarium/publication_annotation"),
    )
    args = parser.parse_args()

    annotation_path = args.annotation.resolve()
    raw_path = args.raw.resolve()
    shrunk_path = args.shrunk.resolve()
    output_dir = args.output_dir.resolve()

    tables_dir = output_dir / "tables"
    diagnostics_dir = output_dir / "diagnostics"
    tables_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    print("Reading Fusarium gene functional annotation...")
    annotation_fields, annotation_rows = read_tsv(annotation_path)

    print("Reading raw DESeq2 results...")
    raw_fields, raw_rows = read_tsv(raw_path)

    print("Reading shrunken DESeq2 results...")
    shrunk_fields, shrunk_rows = read_tsv(shrunk_path)

    required_annotation = {
        "gene_id",
        "has_predicted_protein",
        "signalp_confident_positive",
        "dbcan_any_hit",
        "dbcan_high_confidence",
    }
    required_raw = {
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    }
    required_shrunk = {
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "pvalue",
        "padj",
    }

    missing = required_annotation - set(annotation_fields)
    if missing:
        fail(f"annotation table missing columns: {sorted(missing)}")

    missing = required_raw - set(raw_fields)
    if missing:
        fail(f"raw DESeq2 table missing columns: {sorted(missing)}")

    missing = required_shrunk - set(shrunk_fields)
    if missing:
        fail(f"shrunken DESeq2 table missing columns: {sorted(missing)}")

    annotation_by_id = index_unique(
        annotation_rows, "gene_id", "annotation table"
    )
    raw_by_id = index_unique(raw_rows, "gene_id", "raw DESeq2 table")
    shrunk_by_id = index_unique(
        shrunk_rows, "gene_id", "shrunken DESeq2 table"
    )

    raw_ids = list(raw_by_id)
    shrunk_ids = list(shrunk_by_id)
    annotation_ids = set(annotation_by_id)

    if set(raw_ids) != set(shrunk_ids):
        fail("raw and shrunken DESeq2 gene sets differ")
    if raw_ids != shrunk_ids:
        fail("raw and shrunken DESeq2 gene order differs")

    deseq2_not_annotation = sorted(
        set(raw_ids) - annotation_ids, key=natural_gene_key
    )
    annotation_not_deseq2 = sorted(
        annotation_ids - set(raw_ids), key=natural_gene_key
    )

    comparison_columns = [
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
    ]

    output_fields = (
        ["gene_id"]
        + comparison_columns
        + [field for field in annotation_fields if field != "gene_id"]
    )

    full_rows: List[Dict[str, object]] = []
    tested_rows: List[Dict[str, object]] = []
    significant_rows: List[Dict[str, object]] = []
    significant_raw_lfc1_rows: List[Dict[str, object]] = []
    significant_shrunk_lfc1_rows: List[Dict[str, object]] = []

    genes_with_padj = 0
    genes_with_na_padj = 0
    significant = 0
    significant_raw_lfc1 = 0
    significant_shrunk_lfc1 = 0
    significant_with_protein = 0

    for gene_id in sorted(annotation_by_id, key=natural_gene_key):
        annotation = dict(annotation_by_id[gene_id])
        raw = raw_by_id.get(gene_id)
        shrunk = shrunk_by_id.get(gene_id)

        if raw is None:
            joined: Dict[str, object] = {
                "gene_id": gene_id,
                "in_deseq2_dataset": "FALSE",
                "padj_available": "FALSE",
                "deseq2_status": "not_in_deseq2_dataset",
                "regulation": "not_tested",
                "significant_padj_lt_0.05": "FALSE",
                "significant_raw_abs_lfc_gt_1": "FALSE",
                "significant_shrunk_abs_lfc_gt_1": "FALSE",
                "raw_baseMean": "",
                "raw_log2FoldChange": "",
                "raw_lfcSE": "",
                "raw_stat": "",
                "raw_pvalue": "",
                "raw_padj": "",
                "shrunk_log2FoldChange": "",
                "shrunk_lfcSE": "",
            }
        else:
            if shrunk is None:
                fail(f"missing shrunken row for {gene_id}")

            raw_pvalue = parse_float(raw["pvalue"])
            raw_padj = parse_float(raw["padj"])
            shrunk_pvalue = parse_float(shrunk["pvalue"])
            shrunk_padj = parse_float(shrunk["padj"])
            raw_lfc = parse_float(raw["log2FoldChange"])
            shrunk_lfc = parse_float(shrunk["log2FoldChange"])

            if raw_pvalue != shrunk_pvalue or raw_padj != shrunk_padj:
                fail(f"raw and shrunken p-values differ for {gene_id}")
            if raw_lfc is None or shrunk_lfc is None:
                fail(f"missing fold change for {gene_id}")

            padj_available = raw_padj is not None
            sig = padj_available and raw_padj < 0.05
            sig_raw_lfc1 = sig and abs(raw_lfc) > 1
            sig_shrunk_lfc1 = sig and abs(shrunk_lfc) > 1

            if padj_available:
                genes_with_padj += 1
                status = "significant" if sig else "not_significant"
            else:
                genes_with_na_padj += 1
                status = "independent_filtered_or_undefined_padj"

            if sig:
                significant += 1
                regulation = (
                    "upregulated" if raw_lfc > 0 else "downregulated"
                )
                if (
                    annotation["has_predicted_protein"].strip().upper()
                    == "TRUE"
                ):
                    significant_with_protein += 1
            elif padj_available:
                regulation = "not_significant"
            else:
                regulation = "not_tested"

            if sig_raw_lfc1:
                significant_raw_lfc1 += 1
            if sig_shrunk_lfc1:
                significant_shrunk_lfc1 += 1

            joined = {
                "gene_id": gene_id,
                "in_deseq2_dataset": "TRUE",
                "padj_available": bool_text(padj_available),
                "deseq2_status": status,
                "regulation": regulation,
                "significant_padj_lt_0.05": bool_text(sig),
                "significant_raw_abs_lfc_gt_1": bool_text(sig_raw_lfc1),
                "significant_shrunk_abs_lfc_gt_1": bool_text(
                    sig_shrunk_lfc1
                ),
                "raw_baseMean": raw["baseMean"],
                "raw_log2FoldChange": raw["log2FoldChange"],
                "raw_lfcSE": raw["lfcSE"],
                "raw_stat": raw["stat"],
                "raw_pvalue": raw["pvalue"],
                "raw_padj": raw["padj"],
                "shrunk_log2FoldChange": shrunk["log2FoldChange"],
                "shrunk_lfcSE": shrunk["lfcSE"],
            }

        joined.update(
            {
                field: annotation[field]
                for field in annotation_fields
                if field != "gene_id"
            }
        )
        full_rows.append(joined)

        if joined["in_deseq2_dataset"] == "TRUE":
            tested_rows.append(joined)
        if joined["significant_padj_lt_0.05"] == "TRUE":
            significant_rows.append(joined)
        if joined["significant_raw_abs_lfc_gt_1"] == "TRUE":
            significant_raw_lfc1_rows.append(joined)
        if joined["significant_shrunk_abs_lfc_gt_1"] == "TRUE":
            significant_shrunk_lfc1_rows.append(joined)

    output_paths = {
        "full": tables_dir / "fusarium_publication_annotation_full.tsv",
        "deseq2_dataset": (
            tables_dir / "fusarium_publication_annotation_deseq2_dataset.tsv"
        ),
        "significant": (
            tables_dir / "fusarium_publication_annotation_significant.tsv"
        ),
        "significant_raw_lfc1": (
            tables_dir
            / "fusarium_publication_annotation_significant_raw_lfc1.tsv"
        ),
        "significant_shrunk_lfc1": (
            tables_dir
            / "fusarium_publication_annotation_significant_shrunk_lfc1.tsv"
        ),
    }

    write_tsv(output_paths["full"], output_fields, full_rows)
    write_tsv(
        output_paths["deseq2_dataset"], output_fields, tested_rows
    )
    write_tsv(
        output_paths["significant"], output_fields, significant_rows
    )
    write_tsv(
        output_paths["significant_raw_lfc1"],
        output_fields,
        significant_raw_lfc1_rows,
    )
    write_tsv(
        output_paths["significant_shrunk_lfc1"],
        output_fields,
        significant_shrunk_lfc1_rows,
    )

    write_tsv(
        diagnostics_dir / "deseq2_ids_not_in_annotation.tsv",
        ["gene_id"],
        ({"gene_id": gene_id} for gene_id in deseq2_not_annotation),
    )
    write_tsv(
        diagnostics_dir / "annotation_ids_not_in_deseq2.tsv",
        ["gene_id"],
        ({"gene_id": gene_id} for gene_id in annotation_not_deseq2),
    )

    observed = {
        "annotation_gene_rows": len(annotation_rows),
        "raw_deseq2_rows": len(raw_rows),
        "shrunk_deseq2_rows": len(shrunk_rows),
        "genes_with_padj": genes_with_padj,
        "genes_with_na_padj": genes_with_na_padj,
        "significant_padj_lt_0.05": significant,
        "significant_raw_abs_lfc_gt_1": significant_raw_lfc1,
        "significant_shrunk_abs_lfc_gt_1": significant_shrunk_lfc1,
        "significant_with_predicted_protein": significant_with_protein,
        "deseq2_ids_not_in_annotation": len(deseq2_not_annotation),
        "annotation_ids_not_in_deseq2": len(annotation_not_deseq2),
    }

    validation_rows = []
    all_pass = True
    for check, expected in EXPECTED.items():
        value = observed[check]
        status = "PASS" if value == expected else "FAIL"
        if status == "FAIL":
            all_pass = False
        validation_rows.append(
            {
                "check": check,
                "value": value,
                "expected": expected,
                "status": status,
            }
        )

    validation_path = diagnostics_dir / "validation.tsv"
    write_tsv(
        validation_path,
        ["check", "value", "expected", "status"],
        validation_rows,
    )

    run_info_path = output_dir / "run_info.tsv"
    run_info_rows = [
        {"key": "generated_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"key": "python_version", "value": platform.python_version()},
        {"key": "script", "value": str(Path(__file__).resolve())},
        {"key": "annotation_input", "value": str(annotation_path)},
        {"key": "raw_deseq2_input", "value": str(raw_path)},
        {"key": "shrunk_deseq2_input", "value": str(shrunk_path)},
        {"key": "significance_rule", "value": "raw padj < 0.05"},
        {"key": "raw_lfc1_rule", "value": "raw padj < 0.05 and abs(raw LFC) > 1"},
        {
            "key": "shrunk_lfc1_rule",
            "value": "raw padj < 0.05 and abs(shrunken LFC) > 1",
        },
    ]
    write_tsv(run_info_path, ["key", "value"], run_info_rows)

    checksum_targets = [
        diagnostics_dir / "annotation_ids_not_in_deseq2.tsv",
        diagnostics_dir / "deseq2_ids_not_in_annotation.tsv",
        validation_path,
        run_info_path,
        *output_paths.values(),
    ]
    checksum_path = output_dir / "checksums.sha256"
    with checksum_path.open("w") as handle:
        for path in sorted(checksum_targets):
            relative = path.relative_to(output_dir)
            handle.write(f"{sha256(path)}  {relative}\n")

    print()
    print("Fusarium publication annotation build")
    print()
    print(f"Annotation genes:                  {len(annotation_rows):,}")
    print(f"Genes in DESeq2 dataset:           {len(raw_rows):,}")
    print(f"Genes with non-NA padj:            {genes_with_padj:,}")
    print(f"Genes with NA padj:                {genes_with_na_padj:,}")
    print(f"Significant genes:                 {significant:,}")
    print(f"Significant with protein:          {significant_with_protein:,}")
    print(f"Significant + |raw LFC| > 1:       {significant_raw_lfc1:,}")
    print(f"Significant + |shrunk LFC| > 1:    {significant_shrunk_lfc1:,}")
    print(f"Annotation genes absent from DESeq2: {len(annotation_not_deseq2):,}")
    print()
    print("Validation status:", "PASS" if all_pass else "FAIL")
    print(f"Output directory: {output_dir}")

    if not all_pass:
        fail(f"validation failed; inspect {validation_path}")


if __name__ == "__main__":
    main()
