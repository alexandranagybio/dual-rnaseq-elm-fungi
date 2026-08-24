#!/usr/bin/env python3
"""
Build the complete Ophiostoma novo-ulmi functional annotation table.

This script integrates functional annotations for the complete authoritative
Ophiostoma protein set:

1. eggNOG-mapper orthology and functional annotation
2. SignalP 6.0 signal-peptide predictions
3. dbCAN carbohydrate-active enzyme annotation

The authoritative protein identifiers use the form:

    mRNA_<number>

and map deterministically to DESeq2/GFF3 gene identifiers:

    gene_<number>

All 8,640 authoritative proteins are retained. Proteins without a final
eggNOG annotation row receive:

    eggnog_annotated = FALSE

with blank eggNOG annotation fields.

Run from the repository root:

    python workflow/scripts/13_build_ophiostoma_functional_annotation.py
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path.cwd()

FASTA = (
    ROOT
    / "data/reference/ophiostoma/"
      "ophiostoma.proteins.clean.fa"
)

EGGNOG = (
    ROOT
    / "results/ophiostoma/eggnog_complete/"
      "ophiostoma_complete.emapper.annotations"
)

SIGNALP = (
    ROOT
    / "results/ophiostoma/signalp_raw/"
      "prediction_results.txt"
)

DBCAN = (
    ROOT
    / "results/ophiostoma/dbcan_raw/"
      "overview.tsv"
)

OUTDIR = (
    ROOT
    / "results/ophiostoma/functional_annotation"
)

OUT_TSV = (
    OUTDIR
    / "ophiostoma_functional_annotation.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "validation_summary.txt"
)

OUT_VALIDATION = (
    OUTDIR
    / "validation.tsv"
)

OUT_RUN_INFO = (
    OUTDIR
    / "run_info.tsv"
)

OUT_CHECKSUMS = (
    OUTDIR
    / "checksums.sha256"
)

EXPECTED_PROTEINS = 8640
EXPECTED_EGGNOG_ANNOTATED = 8231
EXPECTED_SIGNALP_POSITIVE = 591
EXPECTED_DBCAN_ANY_HIT = 1440
EXPECTED_DBCAN_HIGH_CONFIDENCE = 315

MRNA_RE = re.compile(r"^mRNA_(\d+)$")


EGGNOG_SOURCE_COLUMNS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]


EGGNOG_OUTPUT_MAP = {
    "seed_ortholog": "eggnog_seed_ortholog",
    "evalue": "eggnog_evalue",
    "score": "eggnog_score",
    "eggNOG_OGs": "eggnog_eggnog_ogs",
    "max_annot_lvl": "eggnog_max_annot_lvl",
    "COG_category": "eggnog_cog_category",
    "Description": "eggnog_description",
    "Preferred_name": "eggnog_preferred_name",
    "GOs": "eggnog_gos",
    "EC": "eggnog_ec",
    "KEGG_ko": "eggnog_kegg_ko",
    "KEGG_Pathway": "eggnog_kegg_pathway",
    "KEGG_Module": "eggnog_kegg_module",
    "KEGG_Reaction": "eggnog_kegg_reaction",
    "KEGG_rclass": "eggnog_kegg_rclass",
    "BRITE": "eggnog_brite",
    "KEGG_TC": "eggnog_kegg_tc",
    "CAZy": "eggnog_cazy",
    "BiGG_Reaction": "eggnog_bigg_reaction",
    "PFAMs": "eggnog_pfams",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty {label}: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


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


def natural_mrna_key(mrna_id: str) -> int:
    match = MRNA_RE.fullmatch(mrna_id)

    if not match:
        fail(f"unexpected protein identifier: {mrna_id}")

    return int(match.group(1))


def gene_id(mrna_id: str) -> str:
    return f"gene_{natural_mrna_key(mrna_id)}"


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []

    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split()[0])

    if not identifiers:
        fail(f"no FASTA identifiers found in {path}")

    if len(identifiers) != len(set(identifiers)):
        fail("duplicated FASTA identifiers")

    for identifier in identifiers:
        natural_mrna_key(identifier)

    return identifiers


def normalize_annotation_value(value: str | None) -> str:
    if value is None:
        return ""

    cleaned = value.strip()

    if cleaned in {"", "-", "NA", "N/A", "None", "none"}:
        return ""

    return cleaned


def read_eggnog(path: Path) -> dict[str, dict[str, str]]:
    """
    Read the eggNOG-mapper .emapper.annotations output.

    The true header starts with '#query'. Other comment and metadata lines
    beginning with '#' are ignored.
    """
    header: list[str] | None = None
    annotation_rows: list[list[str]] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#query\t"):
                header = line[1:].split("\t")
                continue

            if line.startswith("#"):
                continue

            annotation_rows.append(line.split("\t"))

    if header is None:
        fail(f"eggNOG header beginning with #query not found: {path}")

    missing_columns = set(EGGNOG_SOURCE_COLUMNS) - set(header)

    if missing_columns:
        fail(
            "eggNOG annotation table is missing columns: "
            f"{sorted(missing_columns)}"
        )

    query_index = header.index("query")
    output: dict[str, dict[str, str]] = {}

    for values in annotation_rows:
        if len(values) != len(header):
            fail(
                "eggNOG row has unexpected column count: "
                f"expected {len(header)}, observed {len(values)}"
            )

        source_row = dict(zip(header, values))
        query = source_row["query"].strip()

        if not query:
            fail("blank query identifier in eggNOG annotation table")

        natural_mrna_key(query)

        if query in output:
            fail(f"duplicated eggNOG query identifier: {query}")

        annotation = {
            "eggnog_annotated": "TRUE",
        }

        for source_field, output_field in EGGNOG_OUTPUT_MAP.items():
            annotation[output_field] = normalize_annotation_value(
                source_row.get(source_field, "")
            )

        output[query] = annotation

    return output


def blank_eggnog() -> dict[str, str]:
    row = {
        "eggnog_annotated": "FALSE",
    }

    for output_field in EGGNOG_OUTPUT_MAP.values():
        row[output_field] = ""

    return row


def read_signalp(path: Path) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}

    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue

            if len(row) < 4:
                fail(
                    "SignalP row contains fewer than four columns: "
                    f"{row}"
                )

            protein_id, prediction, other_score, sp_score = row[:4]
            cleavage_site = row[4] if len(row) > 4 else ""

            natural_mrna_key(protein_id)

            if protein_id in output:
                fail(f"duplicated SignalP identifier: {protein_id}")

            if prediction not in {"OTHER", "SP"}:
                fail(
                    "unexpected SignalP prediction for "
                    f"{protein_id}: {prediction}"
                )

            output[protein_id] = {
                "signalp_prediction": prediction,
                "signalp_is_sp": str(prediction == "SP").upper(),
                "signalp_other_score": other_score,
                "signalp_sp_score": sp_score,
                "signalp_cs_position": cleavage_site,
            }

    return output


def read_dbcan(path: Path) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            fail(f"no dbCAN header found in {path}")

        required = {
            "Gene ID",
            "#ofTools",
            "EC#",
            "dbCAN_hmm",
            "dbCAN_sub",
            "DIAMOND",
            "Recommend Results",
            "Substrate",
        }

        missing = required - set(reader.fieldnames)

        if missing:
            fail(f"dbCAN table is missing columns: {sorted(missing)}")

        for row in reader:
            protein_id = row["Gene ID"].strip()
            natural_mrna_key(protein_id)

            if protein_id in output:
                fail(f"duplicated dbCAN identifier: {protein_id}")

            try:
                number_of_tools = int(row["#ofTools"])
            except ValueError as exc:
                fail(
                    f"invalid dbCAN #ofTools for {protein_id}: "
                    f"{row['#ofTools']}"
                )

            output[protein_id] = {
                "dbcan_any_hit": "TRUE",
                "dbcan_high_confidence": str(
                    number_of_tools >= 2
                ).upper(),
                "dbcan_n_tools": str(number_of_tools),
                "dbcan_ec": row["EC#"],
                "dbcan_hmm": row["dbCAN_hmm"],
                "dbcan_sub": row["dbCAN_sub"],
                "dbcan_diamond": row["DIAMOND"],
                "dbcan_recommended": row["Recommend Results"],
                "dbcan_substrate": row["Substrate"],
            }

    return output


def blank_dbcan() -> dict[str, str]:
    return {
        "dbcan_any_hit": "FALSE",
        "dbcan_high_confidence": "FALSE",
        "dbcan_n_tools": "0",
        "dbcan_ec": "-",
        "dbcan_hmm": "-",
        "dbcan_sub": "-",
        "dbcan_diamond": "-",
        "dbcan_recommended": "-",
        "dbcan_substrate": "-",
    }


def main() -> None:
    inputs = {
        "authoritative protein FASTA": FASTA,
        "eggNOG annotation": EGGNOG,
        "SignalP annotation": SIGNALP,
        "dbCAN annotation": DBCAN,
    }

    for label, path in inputs.items():
        require_file(path, label)

    ids = fasta_ids(FASTA)
    fasta_set = set(ids)

    eggnog = read_eggnog(EGGNOG)
    signalp = read_signalp(SIGNALP)
    dbcan = read_dbcan(DBCAN)

    if len(ids) != EXPECTED_PROTEINS:
        fail(
            f"expected {EXPECTED_PROTEINS} proteins, "
            f"observed {len(ids)}"
        )

    if len(eggnog) != EXPECTED_EGGNOG_ANNOTATED:
        fail(
            f"expected {EXPECTED_EGGNOG_ANNOTATED} final eggNOG "
            f"annotation rows, observed {len(eggnog)}"
        )

    if not set(eggnog).issubset(fasta_set):
        unexpected = sorted(
            set(eggnog) - fasta_set,
            key=natural_mrna_key,
        )

        fail(
            f"eggNOG contains {len(unexpected)} identifiers absent "
            f"from FASTA; examples: {unexpected[:10]}"
        )

    if set(signalp) != fasta_set:
        only_fasta = sorted(
            fasta_set - set(signalp),
            key=natural_mrna_key,
        )

        only_signalp = sorted(
            set(signalp) - fasta_set,
            key=natural_mrna_key,
        )

        fail(
            "SignalP identifiers do not exactly match FASTA IDs; "
            f"only FASTA={len(only_fasta)}, "
            f"only SignalP={len(only_signalp)}"
        )

    if not set(dbcan).issubset(fasta_set):
        unexpected = sorted(
            set(dbcan) - fasta_set,
            key=natural_mrna_key,
        )

        fail(
            f"dbCAN contains {len(unexpected)} identifiers absent "
            f"from FASTA; examples: {unexpected[:10]}"
        )

    rows: list[dict[str, str]] = []

    for protein_id in sorted(ids, key=natural_mrna_key):
        row = {
            "gene_id": gene_id(protein_id),
            "mrna_id": protein_id,
        }

        row.update(
            eggnog.get(
                protein_id,
                blank_eggnog(),
            )
        )

        row.update(signalp[protein_id])

        row.update(
            dbcan.get(
                protein_id,
                blank_dbcan(),
            )
        )

        rows.append(row)

    if len(rows) != EXPECTED_PROTEINS:
        fail(
            f"expected {EXPECTED_PROTEINS} output rows, "
            f"observed {len(rows)}"
        )

    output_gene_ids = [row["gene_id"] for row in rows]

    if len(output_gene_ids) != len(set(output_gene_ids)):
        fail("duplicated output gene identifiers")

    eggnog_count = sum(
        row["eggnog_annotated"] == "TRUE"
        for row in rows
    )

    signalp_positive_count = sum(
        row["signalp_is_sp"] == "TRUE"
        for row in rows
    )

    dbcan_any_count = sum(
        row["dbcan_any_hit"] == "TRUE"
        for row in rows
    )

    dbcan_high_count = sum(
        row["dbcan_high_confidence"] == "TRUE"
        for row in rows
    )

    if eggnog_count != EXPECTED_EGGNOG_ANNOTATED:
        fail(
            f"expected {EXPECTED_EGGNOG_ANNOTATED} eggNOG-annotated "
            f"proteins in output, observed {eggnog_count}"
        )

    if signalp_positive_count != EXPECTED_SIGNALP_POSITIVE:
        fail(
            f"expected {EXPECTED_SIGNALP_POSITIVE} SignalP-positive "
            f"proteins, observed {signalp_positive_count}"
        )

    if dbcan_any_count != EXPECTED_DBCAN_ANY_HIT:
        fail(
            f"expected {EXPECTED_DBCAN_ANY_HIT} proteins with any "
            f"dbCAN hit, observed {dbcan_any_count}"
        )

    if dbcan_high_count != EXPECTED_DBCAN_HIGH_CONFIDENCE:
        fail(
            f"expected {EXPECTED_DBCAN_HIGH_CONFIDENCE} "
            f"high-confidence dbCAN proteins, observed "
            f"{dbcan_high_count}"
        )

    support = Counter(
        row["dbcan_n_tools"]
        for row in rows
    )

    OUTDIR.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys())

    write_tsv(
        OUT_TSV,
        fieldnames,
        rows,
    )

    validation_rows = [
        {
            "check": "authoritative_fasta_protein_count",
            "value": len(ids),
            "expected": EXPECTED_PROTEINS,
            "status": "PASS",
        },
        {
            "check": "output_gene_count",
            "value": len(rows),
            "expected": EXPECTED_PROTEINS,
            "status": "PASS",
        },
        {
            "check": "unique_output_gene_count",
            "value": len(set(output_gene_ids)),
            "expected": EXPECTED_PROTEINS,
            "status": "PASS",
        },
        {
            "check": "eggnog_annotation_count",
            "value": eggnog_count,
            "expected": EXPECTED_EGGNOG_ANNOTATED,
            "status": "PASS",
        },
        {
            "check": "eggnog_unannotated_count",
            "value": len(rows) - eggnog_count,
            "expected": EXPECTED_PROTEINS - EXPECTED_EGGNOG_ANNOTATED,
            "status": "PASS",
        },
        {
            "check": "signalp_prediction_count",
            "value": len(signalp),
            "expected": EXPECTED_PROTEINS,
            "status": "PASS",
        },
        {
            "check": "signalp_positive_count",
            "value": signalp_positive_count,
            "expected": EXPECTED_SIGNALP_POSITIVE,
            "status": "PASS",
        },
        {
            "check": "dbcan_any_hit_count",
            "value": dbcan_any_count,
            "expected": EXPECTED_DBCAN_ANY_HIT,
            "status": "PASS",
        },
        {
            "check": "dbcan_high_confidence_count",
            "value": dbcan_high_count,
            "expected": EXPECTED_DBCAN_HIGH_CONFIDENCE,
            "status": "PASS",
        },
    ]

    write_tsv(
        OUT_VALIDATION,
        ["check", "value", "expected", "status"],
        validation_rows,
    )

    summary = [
        "Ophiostoma functional annotation validation",
        "",
        f"Authoritative FASTA proteins: {len(ids)}",
        f"Complete output genes: {len(rows)}",
        "",
        f"eggNOG final annotation rows: {eggnog_count}",
        f"Proteins without final eggNOG annotation: "
        f"{len(rows) - eggnog_count}",
        "",
        f"SignalP predictions: {len(signalp)}",
        f"SignalP SP predictions: {signalp_positive_count}",
        f"SignalP OTHER predictions: "
        f"{len(rows) - signalp_positive_count}",
        "",
        f"dbCAN proteins with >=1 method hit: {dbcan_any_count}",
        f"dbCAN high-confidence proteins (#ofTools >=2): "
        f"{dbcan_high_count}",
        f"Proteins with no dbCAN hit: "
        f"{len(rows) - dbcan_any_count}",
        "",
        "dbCAN support distribution across complete proteome:",
    ]

    for number_of_tools in sorted(support, key=int):
        summary.append(
            f"  {number_of_tools} tool(s): "
            f"{support[number_of_tools]}"
        )

    summary.extend(
        [
            "",
            "Validation status: PASS",
            f"Output: {OUT_TSV}",
            f"Validation table: {OUT_VALIDATION}",
        ]
    )

    OUT_SUMMARY.write_text(
        "\n".join(summary) + "\n"
    )

    timestamp = datetime.now(timezone.utc).isoformat()

    run_info_rows = [
        {
            "field": "script",
            "value": (
                "13_build_ophiostoma_functional_annotation.py"
            ),
        },
        {
            "field": "run_timestamp_utc",
            "value": timestamp,
        },
        {
            "field": "authoritative_fasta",
            "value": str(FASTA),
        },
        {
            "field": "eggnog_annotation",
            "value": str(EGGNOG),
        },
        {
            "field": "signalp_annotation",
            "value": str(SIGNALP),
        },
        {
            "field": "dbcan_annotation",
            "value": str(DBCAN),
        },
        {
            "field": "output_table",
            "value": str(OUT_TSV),
        },
        {
            "field": "authoritative_fasta_sha256",
            "value": sha256(FASTA),
        },
        {
            "field": "eggnog_annotation_sha256",
            "value": sha256(EGGNOG),
        },
        {
            "field": "signalp_annotation_sha256",
            "value": sha256(SIGNALP),
        },
        {
            "field": "dbcan_annotation_sha256",
            "value": sha256(DBCAN),
        },
        {
            "field": "output_table_sha256",
            "value": sha256(OUT_TSV),
        },
    ]

    write_tsv(
        OUT_RUN_INFO,
        ["field", "value"],
        run_info_rows,
    )

    checksum_paths = [
        FASTA,
        EGGNOG,
        SIGNALP,
        DBCAN,
        OUT_TSV,
        OUT_VALIDATION,
        OUT_SUMMARY,
        OUT_RUN_INFO,
    ]

    checksum_lines = [
        f"{sha256(path)}  {path}"
        for path in checksum_paths
    ]

    OUT_CHECKSUMS.write_text(
        "\n".join(checksum_lines) + "\n"
    )

    print("\n".join(summary))


if __name__ == "__main__":
    main()
