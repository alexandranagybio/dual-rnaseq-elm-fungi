#!/usr/bin/env python3
"""
Build and validate the Ophiostoma novo-ulmi gene-to-annotation map.

This script:

1. Parses only `gene` features from the original GFF3.
2. Extracts gene_id, protein_id, transcript_id, JGI name, and coordinates.
3. Validates that exactly 8,640 unique genes are present.
4. Validates that every gene maps to exactly one protein ID.
5. Joins the structural map to JGI GO and KEGG annotation tables.
6. Annotates the three DESeq2 contrast tables.
7. Reports unmatched identifiers and duplicated mappings.

The script does not modify counts or rerun DESeq2.

Default paths assume the repository layout used for the rebuilt workflow and the
historical annotation files under ~/rnaseq.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, TextIO


EXPECTED_GENE_COUNT = 8640


DEFAULT_GFF3 = Path(
    "data/annotation/ophiostoma/annotation.gff3"
)

DEFAULT_GO = Path(
    "data/annotation/ophiostoma/GO/"
    "Ophnu1_GeneCatalog_proteins_20170425_GO.tab"
)

DEFAULT_KEGG = Path(
    "data/annotation/ophiostoma/KEGG/"
    "Ophnu1_GeneCatalog_proteins_20170425_KEGG.tab"
)

DEFAULT_DE_DIR = Path(
    "results/ophiostoma/deseq2_results/tables"
)

DEFAULT_DE_DIR = Path("results/ophiostoma/deseq2_results")
DEFAULT_OUT_DIR = Path("results/ophiostoma/annotation")



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and validate the Ophiostoma annotation map."
    )
    parser.add_argument("--gff3", type=Path, default=DEFAULT_GFF3)
    parser.add_argument("--go", type=Path, default=DEFAULT_GO)
    parser.add_argument("--kegg", type=Path, default=DEFAULT_KEGG)
    parser.add_argument("--de-dir", type=Path, default=DEFAULT_DE_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--expected-genes",
        type=int,
        default=EXPECTED_GENE_COUNT,
        help="Expected number of unique gene features in the GFF3.",
    )
    parser.add_argument(
        "--de-files",
        nargs="*",
        type=Path,
        default=None,
        help=(
            "Optional explicit DESeq2 result tables. When omitted, the script "
            "discovers the three contrast tables in --de-dir."
        ),
    )
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def parse_gff3_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def natural_gene_key(gene_id: str) -> tuple:
    match = re.fullmatch(r"(.+?)(\d+)", gene_id)
    if match:
        return (match.group(1), int(match.group(2)))
    return (gene_id, -1)


def parse_gene_features(
    gff3_path: Path,
) -> tuple[list[dict[str, str]], dict[str, list[str]], dict[str, list[str]]]:
    """
    Return:
      - one canonical record per unique gene ID;
      - duplicate gene feature records;
      - protein IDs associated with multiple genes.
    """
    records_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)

    with open_text(gff3_path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed GFF3 line {line_no}: expected 9 fields, "
                    f"observed {len(fields)}."
                )

            seqid, source, feature_type, start, end, score, strand, phase, raw_attrs = fields

            if feature_type != "gene":
                continue

            attrs = parse_gff3_attributes(raw_attrs)
            gene_id = attrs.get("ID", "").strip()
            protein_id = attrs.get("proteinId", "").strip()
            transcript_id = attrs.get("transcriptId", "").strip()
            jgi_name = attrs.get("Name", "").strip()
            portal_id = attrs.get("portal_id", "").strip()

            if not gene_id:
                raise ValueError(f"Gene feature at line {line_no} has no ID attribute.")

            record = {
                "gene_id": gene_id,
                "protein_id": protein_id,
                "transcript_id": transcript_id,
                "jgi_name": jgi_name,
                "portal_id": portal_id,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "source": source,
                "gff3_line": str(line_no),
            }
            records_by_gene[gene_id].append(record)

    if not records_by_gene:
        raise ValueError(f"No gene features found in {gff3_path}")

    duplicate_gene_records: dict[str, list[str]] = {}
    canonical_records: list[dict[str, str]] = []

    for gene_id, records in records_by_gene.items():
        if len(records) > 1:
            duplicate_gene_records[gene_id] = [r["gff3_line"] for r in records]

        unique_proteins = {r["protein_id"] for r in records if r["protein_id"]}
        if len(unique_proteins) > 1:
            raise ValueError(
                f"Gene {gene_id} maps to multiple protein IDs: "
                f"{sorted(unique_proteins)}"
            )

        canonical_records.append(records[0])

    genes_by_protein: dict[str, list[str]] = defaultdict(list)
    for record in canonical_records:
        if record["protein_id"]:
            genes_by_protein[record["protein_id"]].append(record["gene_id"])

    reused_protein_ids = {
        protein_id: sorted(gene_ids, key=natural_gene_key)
        for protein_id, gene_ids in genes_by_protein.items()
        if len(set(gene_ids)) > 1
    }

    canonical_records.sort(key=lambda row: natural_gene_key(row["gene_id"]))
    return canonical_records, duplicate_gene_records, reused_protein_ids


def sniff_delimiter(path: Path) -> str:
    with open_text(path) as handle:
        sample = "".join(handle.readline() for _ in range(10))
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        return dialect.delimiter
    except csv.Error:
        return "\t"


def read_annotation_table(
    path: Path,
    annotation_name: str,
) -> tuple[
    dict[str, list[dict[str, str]]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[str],
]:
    """
    Read a GO or KEGG table indexed by proteinId.

    Returns:
      mappings_by_protein,
      rows lacking a protein ID,
      exact duplicate rows,
      original field names.
    """
    delimiter = sniff_delimiter(path)

    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"Empty annotation table: {path}")

        first = first.lstrip("\ufeff")
        if first.startswith("#"):
            first = first[1:]

        remaining = handle.read()

    lines = [first] + remaining.splitlines(keepends=True)
    reader = csv.DictReader(lines, delimiter=delimiter)

    if not reader.fieldnames:
        raise ValueError(f"No header found in {path}")

    fieldnames = [str(name).strip().lstrip("#") for name in reader.fieldnames]
    reader.fieldnames = fieldnames

    protein_candidates = [
        name
        for name in fieldnames
        if name.lower().replace("_", "") in {"proteinid", "protein"}
    ]
    if not protein_candidates:
        raise ValueError(
            f"{annotation_name} table has no recognizable proteinId column. "
            f"Columns: {fieldnames}"
        )
    protein_col = protein_candidates[0]

    mappings: dict[str, list[dict[str, str]]] = defaultdict(list)
    missing_protein_rows: list[dict[str, str]] = []
    duplicate_rows: list[dict[str, str]] = []
    seen_rows: set[tuple[tuple[str, str], ...]] = set()

    for row_number, raw_row in enumerate(reader, start=2):
        row = {
            str(key).strip().lstrip("#"): (value or "").strip()
            for key, value in raw_row.items()
            if key is not None
        }
        row["_source_row"] = str(row_number)

        protein_id = row.get(protein_col, "").strip()
        if not protein_id:
            missing_protein_rows.append(row)
            continue

        row_signature = tuple(
            sorted(
                (key, value)
                for key, value in row.items()
                if key != "_source_row"
            )
        )
        if row_signature in seen_rows:
            duplicate_rows.append(row)
            continue

        seen_rows.add(row_signature)
        mappings[protein_id].append(row)

    return mappings, missing_protein_rows, duplicate_rows, fieldnames


def unique_join(values: Iterable[str], separator: str = "; ") -> str:
    observed: set[str] = set()
    ordered: list[str] = []
    for raw_value in values:
        value = str(raw_value).strip()
        if not value or value in observed:
            continue
        observed.add(value)
        ordered.append(value)

    def natural_key(text: str) -> tuple:
        return tuple(
            (0, int(part)) if part.isdigit() else (1, part.casefold())
            for part in re.split(r"(\d+)", text)
        )

    return separator.join(sorted(ordered, key=natural_key))


def find_column(fieldnames: list[str], candidates: list[str]) -> str | None:
    normalized = {
        name.lower().replace("_", "").replace("-", ""): name
        for name in fieldnames
    }
    for candidate in candidates:
        key = candidate.lower().replace("_", "").replace("-", "")
        if key in normalized:
            return normalized[key]
    return None


def summarize_go(rows: list[dict[str, str]], fieldnames: list[str]) -> dict[str, str]:
    go_acc_col = find_column(fieldnames, ["goAcc", "go_id", "go"])
    go_name_col = find_column(fieldnames, ["goName", "term", "description"])
    go_type_col = find_column(fieldnames, ["gotermType", "ontology", "namespace"])

    return {
        "go_count": str(len(rows)),
        "go_ids": unique_join(
            row.get(go_acc_col, "") for row in rows
        ) if go_acc_col else "",
        "go_names": unique_join(
            row.get(go_name_col, "") for row in rows
        ) if go_name_col else "",
        "go_types": unique_join(
            row.get(go_type_col, "") for row in rows
        ) if go_type_col else "",
    }


def summarize_kegg(rows: list[dict[str, str]], fieldnames: list[str]) -> dict[str, str]:
    non_protein_fields = [
        field
        for field in fieldnames
        if field.lower().replace("_", "") not in {"proteinid", "protein"}
    ]

    summary: dict[str, str] = {"kegg_count": str(len(rows))}

    # Preserve all KEGG columns as collapsed annotation fields so that the
    # script remains compatible with the exact JGI KEGG table schema.
    for field in non_protein_fields:
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", field).strip("_").lower()
        summary[f"kegg_{safe_name}"] = unique_join(
            row.get(field, "") for row in rows
        )

    return summary


def build_annotation_map(
    gene_records: list[dict[str, str]],
    go_map: dict[str, list[dict[str, str]]],
    go_fields: list[str],
    kegg_map: dict[str, list[dict[str, str]]],
    kegg_fields: list[str],
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    output_rows: list[dict[str, str]] = []
    annotation_by_gene: dict[str, dict[str, str]] = {}

    for gene in gene_records:
        protein_id = gene["protein_id"]
        go_summary = summarize_go(go_map.get(protein_id, []), go_fields)
        kegg_summary = summarize_kegg(kegg_map.get(protein_id, []), kegg_fields)

        combined = {
            **gene,
            **go_summary,
            **kegg_summary,
            "has_go": "yes" if go_map.get(protein_id) else "no",
            "has_kegg": "yes" if kegg_map.get(protein_id) else "no",
        }
        output_rows.append(combined)
        annotation_by_gene[gene["gene_id"]] = combined

    return output_rows, annotation_by_gene


def discover_de_files(de_dir: Path) -> list[Path]:
    if not de_dir.is_dir():
        raise FileNotFoundError(f"DESeq2 result directory not found: {de_dir}")

    tables_dir = de_dir / "tables"
    search_dir = tables_dir if tables_dir.is_dir() else de_dir

    candidates = sorted(
        path
        for path in search_dir.iterdir()
        if path.is_file()
        and path.name.lower().endswith("_all_genes.tsv")
    )

    scored: list[tuple[int, Path]] = []
    for path in candidates:
        name = path.name.lower()
        score = 0
        if "interaction" in name:
            score += 2
        if "self" in name:
            score += 2
        if "onu" in name:
            score += 2
        if "deseq" in name or "result" in name:
            score += 1
        scored.append((score, path))

    likely = [path for score, path in scored if score >= 3]

    if len(likely) == 3:
        contrast_order = {
            "interaction_vs_self_all_genes.tsv": 0,
            "onu_vs_self_all_genes.tsv": 1,
            "interaction_vs_onu_all_genes.tsv": 2,
        }
        return sorted(
            likely,
            key=lambda path: (
                contrast_order.get(path.name.lower(), 99),
                path.name.lower(),
            ),
        )

    # Fall back to inspecting table headers for a gene identifier and DESeq2
    # result columns.
    inspected: list[Path] = []
    for path in candidates:
        delimiter = sniff_delimiter(path)
        with open_text(path) as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            try:
                header = next(reader)
            except StopIteration:
                continue

        normalized = {
            col.strip().lower().replace("_", "").replace(".", "")
            for col in header
        }
        has_gene = any(
            col in normalized for col in {"geneid", "gene", "id"}
        )
        has_lfc = "log2foldchange" in normalized
        has_padj = "padj" in normalized

        if has_gene and has_lfc and has_padj:
            inspected.append(path)

    if len(inspected) != 3:
        found = "\n  ".join(str(path) for path in candidates) or "(none)"
        raise ValueError(
            "Could not identify exactly three DESeq2 contrast tables.\n"
            f"Candidate files:\n  {found}\n"
            "Provide them explicitly with --de-files."
        )

    return sorted(inspected)


def identify_gene_column(fieldnames: list[str]) -> str:
    normalized = {
        name.lower().replace("_", "").replace(".", ""): name
        for name in fieldnames
    }
    for candidate in ("geneid", "gene", "id"):
        if candidate in normalized:
            return normalized[candidate]
    raise ValueError(
        f"Could not identify a gene ID column. Columns: {fieldnames}"
    )


def read_table(path: Path) -> tuple[list[dict[str, str]], list[str], str]:
    delimiter = sniff_delimiter(path)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        fieldnames = [str(name).strip().lstrip("\ufeff") for name in reader.fieldnames]
        reader.fieldnames = fieldnames
        rows = [
            {
                str(key).strip(): (value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            for row in reader
        ]
    return rows, fieldnames, delimiter


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def annotate_de_table(
    de_path: Path,
    out_dir: Path,
    annotation_by_gene: dict[str, dict[str, str]],
) -> dict[str, object]:
    rows, de_fields, _ = read_table(de_path)
    gene_col = identify_gene_column(de_fields)

    annotated_rows: list[dict[str, str]] = []
    unmatched_gene_ids: list[str] = []

    annotation_fields: list[str] = []
    if annotation_by_gene:
        first_annotation = next(iter(annotation_by_gene.values()))
        annotation_fields = list(first_annotation.keys())

    for row in rows:
        gene_id = row.get(gene_col, "").strip()
        annotation = annotation_by_gene.get(gene_id)

        if annotation is None:
            unmatched_gene_ids.append(gene_id)
            annotation = {field: "" for field in annotation_fields}

        # Avoid writing gene_id twice when the DE table already contains it.
        annotation_without_duplicate_gene = {
            key: value
            for key, value in annotation.items()
            if key != "gene_id" or gene_col != "gene_id"
        }

        annotated_rows.append({**row, **annotation_without_duplicate_gene})

    output_fields = list(de_fields)
    for field in annotation_fields:
        if field == "gene_id" and gene_col == "gene_id":
            continue
        if field not in output_fields:
            output_fields.append(field)

    output_path = out_dir / f"{de_path.stem}.annotated.tsv"
    write_tsv(output_path, annotated_rows, output_fields)

    unmatched_path = out_dir / f"{de_path.stem}.unmatched_gene_ids.tsv"
    unmatched_rows = [{"gene_id": gene_id} for gene_id in sorted(set(unmatched_gene_ids), key=natural_gene_key)]
    write_tsv(unmatched_path, unmatched_rows, ["gene_id"])

    return {
        "input_file": str(de_path),
        "output_file": str(output_path),
        "rows": len(rows),
        "matched_rows": len(rows) - len(unmatched_gene_ids),
        "unmatched_rows": len(unmatched_gene_ids),
        "unique_unmatched_gene_ids": len(set(unmatched_gene_ids)),
        "gene_id_column": gene_col,
    }


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()

    args.gff3 = args.gff3.expanduser().resolve()
    args.go = args.go.expanduser().resolve()
    args.kegg = args.kegg.expanduser().resolve()
    args.de_dir = args.de_dir.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()

    require_file(args.gff3, "GFF3")
    require_file(args.go, "GO table")
    require_file(args.kegg, "KEGG table")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Parsing gene features from GFF3...")
    gene_records, duplicate_gene_records, reused_protein_ids = parse_gene_features(args.gff3)

    unique_gene_count = len(gene_records)
    genes_missing_protein_id = [
        row["gene_id"] for row in gene_records if not row["protein_id"]
    ]

    if unique_gene_count != args.expected_genes:
        raise ValueError(
            f"Expected {args.expected_genes:,} unique genes, "
            f"but found {unique_gene_count:,}."
        )

    if duplicate_gene_records:
        raise ValueError(
            f"Found {len(duplicate_gene_records):,} duplicated gene IDs in the GFF3. "
            "See the validation report after resolving the structural annotation."
        )

    if genes_missing_protein_id:
        raise ValueError(
            f"{len(genes_missing_protein_id):,} genes have no proteinId attribute. "
            f"First examples: {genes_missing_protein_id[:10]}"
        )

    if reused_protein_ids:
        examples = list(reused_protein_ids.items())[:10]
        raise ValueError(
            f"{len(reused_protein_ids):,} protein IDs map to more than one gene. "
            f"Examples: {examples}"
        )

    print(f"Validated {unique_gene_count:,} unique genes.")
    print("Validated exactly one protein ID per gene.")

    print("Reading GO annotations...")
    go_map, go_missing_protein, go_duplicate_rows, go_fields = read_annotation_table(
        args.go, "GO"
    )

    print("Reading KEGG annotations...")
    kegg_map, kegg_missing_protein, kegg_duplicate_rows, kegg_fields = read_annotation_table(
        args.kegg, "KEGG"
    )

    annotation_rows, annotation_by_gene = build_annotation_map(
        gene_records,
        go_map,
        go_fields,
        kegg_map,
        kegg_fields,
    )

    annotation_fields = list(annotation_rows[0].keys())
    map_path = args.out_dir / "ophiostoma_gene_annotation_map.tsv"
    write_tsv(map_path, annotation_rows, annotation_fields)

    structural_proteins = {row["protein_id"] for row in gene_records}
    go_proteins = set(go_map)
    kegg_proteins = set(kegg_map)

    unmatched_go_proteins = sorted(go_proteins - structural_proteins)
    unmatched_kegg_proteins = sorted(kegg_proteins - structural_proteins)
    genes_without_go = sorted(
        (
            row["gene_id"]
            for row in gene_records
            if row["protein_id"] not in go_map
        ),
        key=natural_gene_key,
    )
    genes_without_kegg = sorted(
        (
            row["gene_id"]
            for row in gene_records
            if row["protein_id"] not in kegg_map
        ),
        key=natural_gene_key,
    )

    write_tsv(
        args.out_dir / "genes_without_GO.tsv",
        [
            {
                "gene_id": gene_id,
                "protein_id": annotation_by_gene[gene_id]["protein_id"],
            }
            for gene_id in genes_without_go
        ],
        ["gene_id", "protein_id"],
    )
    write_tsv(
        args.out_dir / "genes_without_KEGG.tsv",
        [
            {
                "gene_id": gene_id,
                "protein_id": annotation_by_gene[gene_id]["protein_id"],
            }
            for gene_id in genes_without_kegg
        ],
        ["gene_id", "protein_id"],
    )
    write_tsv(
        args.out_dir / "GO_protein_ids_not_in_GFF3.tsv",
        [{"protein_id": protein_id} for protein_id in unmatched_go_proteins],
        ["protein_id"],
    )
    write_tsv(
        args.out_dir / "KEGG_protein_ids_not_in_GFF3.tsv",
        [{"protein_id": protein_id} for protein_id in unmatched_kegg_proteins],
        ["protein_id"],
    )

    duplicate_report_rows: list[dict[str, str]] = []

    for gene_id, line_numbers in duplicate_gene_records.items():
        duplicate_report_rows.append(
            {
                "mapping_type": "duplicate_gene_feature",
                "identifier": gene_id,
                "details": f"GFF3 lines: {', '.join(line_numbers)}",
            }
        )

    for protein_id, gene_ids in reused_protein_ids.items():
        duplicate_report_rows.append(
            {
                "mapping_type": "protein_id_maps_to_multiple_genes",
                "identifier": protein_id,
                "details": "; ".join(gene_ids),
            }
        )

    for row in go_duplicate_rows:
        duplicate_report_rows.append(
            {
                "mapping_type": "exact_duplicate_GO_row",
                "identifier": row.get("proteinId", row.get("protein_id", "")),
                "details": json.dumps(
                    {k: v for k, v in row.items() if k != "_source_row"},
                    sort_keys=True,
                ),
            }
        )

    for row in kegg_duplicate_rows:
        duplicate_report_rows.append(
            {
                "mapping_type": "exact_duplicate_KEGG_row",
                "identifier": row.get("proteinId", row.get("protein_id", "")),
                "details": json.dumps(
                    {k: v for k, v in row.items() if k != "_source_row"},
                    sort_keys=True,
                ),
            }
        )

    write_tsv(
        args.out_dir / "duplicated_mappings.tsv",
        duplicate_report_rows,
        ["mapping_type", "identifier", "details"],
    )

    if args.de_files:
        de_files = [path.expanduser().resolve() for path in args.de_files]
        for path in de_files:
            require_file(path, "DESeq2 result table")
        if len(de_files) != 3:
            raise ValueError(
                f"--de-files must contain exactly three tables; received {len(de_files)}."
            )
    else:
        de_files = discover_de_files(args.de_dir)

    print("Annotating DESeq2 contrast tables:")
    de_reports: list[dict[str, object]] = []
    for de_file in de_files:
        print(f"  - {de_file}")
        de_reports.append(
            annotate_de_table(de_file, args.out_dir, annotation_by_gene)
        )

    validation = {
        "run_timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "inputs": {
            "gff3": str(args.gff3),
            "go": str(args.go),
            "kegg": str(args.kegg),
            "de_files": [str(path) for path in de_files],
        },
        "input_sha256": {
            "gff3": sha256sum(args.gff3),
            "go": sha256sum(args.go),
            "kegg": sha256sum(args.kegg),
            **{
                f"de_file_{index + 1}": sha256sum(path)
                for index, path in enumerate(de_files)
            },
        },
        "structural_annotation": {
            "expected_unique_genes": args.expected_genes,
            "observed_unique_genes": unique_gene_count,
            "genes_missing_protein_id": len(genes_missing_protein_id),
            "duplicated_gene_ids": len(duplicate_gene_records),
            "protein_ids_mapping_to_multiple_genes": len(reused_protein_ids),
            "validation_passed": (
                unique_gene_count == args.expected_genes
                and not genes_missing_protein_id
                and not duplicate_gene_records
                and not reused_protein_ids
            ),
        },
        "functional_annotation": {
            "genes_with_GO": unique_gene_count - len(genes_without_go),
            "genes_without_GO": len(genes_without_go),
            "genes_with_KEGG": unique_gene_count - len(genes_without_kegg),
            "genes_without_KEGG": len(genes_without_kegg),
            "GO_protein_ids_not_in_GFF3": len(unmatched_go_proteins),
            "KEGG_protein_ids_not_in_GFF3": len(unmatched_kegg_proteins),
            "GO_rows_missing_protein_id": len(go_missing_protein),
            "KEGG_rows_missing_protein_id": len(kegg_missing_protein),
            "exact_duplicate_GO_rows": len(go_duplicate_rows),
            "exact_duplicate_KEGG_rows": len(kegg_duplicate_rows),
        },
        "deseq2_tables": de_reports,
        "outputs": {
            "annotation_map": str(map_path),
            "output_directory": str(args.out_dir),
        },
    }

    validation_path = args.out_dir / "09_build_ophiostoma_annotation_map.validation.json"
    validation_path.write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    checksum_path = args.out_dir / "09_build_ophiostoma_annotation_map.sha256"
    output_files = sorted(
        path for path in args.out_dir.iterdir()
        if path.is_file() and path != checksum_path
    )
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in output_files:
            handle.write(f"{sha256sum(path)}  {path.name}\n")

    print()
    print("Annotation map completed successfully.")
    print(f"Unique genes:             {unique_gene_count:,}")
    print(f"Genes with GO:            {unique_gene_count - len(genes_without_go):,}")
    print(f"Genes without GO:         {len(genes_without_go):,}")
    print(f"Genes with KEGG:          {unique_gene_count - len(genes_without_kegg):,}")
    print(f"Genes without KEGG:       {len(genes_without_kegg):,}")
    print(f"Duplicate GO rows:        {len(go_duplicate_rows):,}")
    print(f"Duplicate KEGG rows:      {len(kegg_duplicate_rows):,}")
    print(f"Results written to:       {args.out_dir}")

    for report in de_reports:
        print(
            f"  {Path(str(report['input_file'])).name}: "
            f"{report['matched_rows']:,}/{report['rows']:,} rows matched"
        )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
