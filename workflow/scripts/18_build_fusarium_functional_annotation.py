#!/usr/bin/env python3
"""
18_build_fusarium_functional_annotation.py

Build and validate protein-, transcript-, and gene-level functional annotation
tables for the de novo Fusarium cf. salinense Trinity assembly.

Inputs
------
1. Trinity gene-transcript map
2. TransDecoder predicted protein FASTA
3. eggNOG-mapper annotation table
4. SignalP 6 prediction table
5. dbCAN overview table

Outputs
-------
results/fusarium/functional_annotation/
├── run_info.tsv
├── checksums.sha256
├── diagnostics/
│   ├── validation.tsv
│   ├── mapping_multiplicity.tsv
│   ├── signalp_ambiguous_no_cleavage_site.tsv
│   ├── eggnog_ids_not_in_proteome.tsv
│   └── dbcan_ids_not_in_proteome.tsv
└── tables/
    ├── fusarium_protein_functional_annotation.tsv
    ├── fusarium_transcript_functional_annotation.tsv
    └── fusarium_gene_functional_annotation.tsv

Definitions
-----------
SignalP raw positive:
    Prediction == SP

SignalP confident positive:
    Prediction == SP and a cleavage-site position is reported

dbCAN any hit:
    Protein is present in the dbCAN overview table

dbCAN high confidence:
    dbCAN #ofTools >= 2

eggNOG annotated:
    Protein has a row in the eggNOG-mapper annotations file

Important interpretation
------------------------
TransDecoder can predict multiple ORFs per Trinity transcript, and Trinity genes
can contain multiple transcripts. Protein-level annotations are therefore
collapsed deterministically to transcript and gene level. Gene-level Boolean
fields indicate whether at least one associated predicted protein satisfies the
criterion. Protein-hit counts retain the underlying multiplicity.

Run from repository root:
    python workflow/scripts/18_build_fusarium_functional_annotation.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO


EXPECTED_PROTEINS = 35_327
EXPECTED_TRANSCRIPTS = 27_979
EXPECTED_GENES = 15_192
EXPECTED_TRANSCRIPTS_WITH_PROTEIN = 19_807
EXPECTED_GENES_WITH_PROTEIN = 9_197
EXPECTED_EGGNOG_ANNOTATED_PROTEINS = 32_248
EXPECTED_DBCAN_PROTEINS = 4_935
EXPECTED_SIGNALP_ROWS = 35_327
EXPECTED_SIGNALP_RAW_SP = 1_706
EXPECTED_SIGNALP_CONFIDENT_SP = 1_703
EXPECTED_SIGNALP_AMBIGUOUS_SP = 3

DEFAULT_PROTEINS = Path(
    "data/external/fusarium_annotation/"
    "Fusarium_pure.Trinity.fasta.transdecoder.pep"
)
DEFAULT_GENE_TRANS_MAP = Path(
    "data/external/fusarium_assembly/"
    "Fusarium_pure.Trinity.fasta.gene_trans_map"
)
DEFAULT_EGGNOG = Path(
    "data/external/fusarium_annotation/"
    "fusarium_complete.emapper.annotations"
)
DEFAULT_SIGNALP = Path(
    "results/fusarium/annotation/signalp_complete/"
    "prediction_results.txt"
)
DEFAULT_DBCAN = Path(
    "data/external/fusarium_annotation/dbcan_overview.tsv"
)
DEFAULT_OUTDIR = Path(
    "results/fusarium/functional_annotation"
)

PROTEIN_ID_RE = re.compile(
    r"^(?P<transcript>TRINITY_.+_i\d+)\|m\.(?P<orf_number>\d+)$"
)
HEADER_META_RE = re.compile(
    r"type:(?P<orf_type>\S+)\s+"
    r"len:(?P<reported_length>\d+)\s+"
    r"\((?P<strand>[+-])\)\s+"
    r"(?P<source_transcript>TRINITY_\S+?):"
    r"(?P<orf_start>\d+)-(?P<orf_end>\d+)"
    r"\((?P<coord_strand>[+-])\)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build validated Fusarium protein-, transcript-, and gene-level "
            "functional annotation tables."
        )
    )
    parser.add_argument(
        "--proteins",
        type=Path,
        default=DEFAULT_PROTEINS,
        help="TransDecoder predicted protein FASTA.",
    )
    parser.add_argument(
        "--gene-trans-map",
        type=Path,
        default=DEFAULT_GENE_TRANS_MAP,
        help="Trinity gene-transcript map.",
    )
    parser.add_argument(
        "--eggnog",
        type=Path,
        default=DEFAULT_EGGNOG,
        help="eggNOG-mapper annotations file.",
    )
    parser.add_argument(
        "--signalp",
        type=Path,
        default=DEFAULT_SIGNALP,
        help="SignalP prediction_results.txt.",
    )
    parser.add_argument(
        "--dbcan",
        type=Path,
        default=DEFAULT_DBCAN,
        help="dbCAN overview.tsv.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory.",
    )
    return parser.parse_args()


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
    rows: Iterable[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", value)
    )


def bool_text(value: bool) -> str:
    return str(value).upper()


def clean_annotation(value: str) -> str:
    value = value.strip()
    return "" if value in {"", "-", "NA", "N/A", "nan", "NaN"} else value


def unique_join(values: Iterable[str], separator: str = "; ") -> str:
    observed: set[str] = set()
    retained: list[str] = []

    for raw in values:
        value = clean_annotation(str(raw))
        if not value:
            continue

        # eggNOG often stores multiple comma-separated values in one cell.
        # Preserve each complete cell here; downstream columns retain the
        # original annotation semantics without attempting risky re-parsing.
        if value not in observed:
            observed.add(value)
            retained.append(value)

    return separator.join(sorted(retained, key=natural_key))


def parse_fasta(
    path: Path,
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, list[str]],
]:
    proteins: dict[str, dict[str, str]] = {}
    proteins_by_transcript: dict[str, list[str]] = defaultdict(list)

    current_id: str | None = None
    current_header = ""
    sequence_parts: list[str] = []

    def store_current() -> None:
        nonlocal current_id, current_header, sequence_parts
        if current_id is None:
            return

        sequence = "".join(sequence_parts).replace(" ", "")
        if not sequence:
            fail(f"empty protein sequence: {current_id}")

        match = PROTEIN_ID_RE.fullmatch(current_id)
        if not match:
            fail(f"unexpected TransDecoder protein ID: {current_id}")

        transcript_id = match.group("transcript")
        orf_number = match.group("orf_number")
        meta_match = HEADER_META_RE.search(current_header)

        orf_type = ""
        reported_length = ""
        strand = ""
        source_transcript = ""
        orf_start = ""
        orf_end = ""

        if meta_match:
            orf_type = meta_match.group("orf_type")
            reported_length = meta_match.group("reported_length")
            strand = meta_match.group("strand")
            source_transcript = meta_match.group("source_transcript")
            orf_start = meta_match.group("orf_start")
            orf_end = meta_match.group("orf_end")

            if source_transcript != transcript_id:
                fail(
                    f"FASTA header transcript mismatch for {current_id}: "
                    f"{source_transcript} != {transcript_id}"
                )
            if meta_match.group("coord_strand") != strand:
                fail(f"FASTA header strand mismatch for {current_id}")

        sequence_without_terminal_stop = (
            sequence[:-1] if sequence.endswith("*") else sequence
        )
        internal_stop_count = sequence_without_terminal_stop.count("*")
        if internal_stop_count:
            fail(
                f"protein contains internal stop codon(s): "
                f"{current_id} ({internal_stop_count})"
            )

        protein_length_aa = len(sequence_without_terminal_stop)
        fasta_sequence_length = len(sequence)

        transdecoder_reported_length = (
            int(reported_length) if reported_length else None
        )

        # TransDecoder uses two systematic header-length conventions:
        #
        # complete / 5prime_partial:
        #     reported len == raw FASTA length
        #
        # internal / 3prime_partial:
        #     reported len == raw FASTA length + 1
        #
        # Preserve both values and validate the convention appropriate
        # for the reported ORF type.
        expected_reported_length = None

        if orf_type in {"complete", "5prime_partial"}:
            expected_reported_length = fasta_sequence_length
        elif orf_type in {"internal", "3prime_partial"}:
            expected_reported_length = fasta_sequence_length + 1
        elif orf_type:
            fail(
                f"unexpected TransDecoder ORF type for "
                f"{current_id}: {orf_type}"
            )

        if (
            transdecoder_reported_length is not None
            and expected_reported_length is not None
            and transdecoder_reported_length
            != expected_reported_length
        ):
            fail(
                f"unexpected TransDecoder length convention for "
                f"{current_id}: type={orf_type}, "
                f"reported={transdecoder_reported_length}, "
                f"raw_fasta={fasta_sequence_length}, "
                f"expected={expected_reported_length}"
            )

        proteins[current_id] = {
            "gene_id": "",
            "transcript_id": transcript_id,
            "protein_id": current_id,
            "transdecoder_orf_number": orf_number,
            "transdecoder_orf_type": orf_type,
            "protein_length_aa": str(protein_length_aa),
            "fasta_sequence_length": str(fasta_sequence_length),
            "transdecoder_reported_length": (
                str(transdecoder_reported_length)
                if transdecoder_reported_length is not None
                else ""
            ),
            "transdecoder_length_difference": (
                str(
                    transdecoder_reported_length
                    - fasta_sequence_length
                )
                if transdecoder_reported_length is not None
                else ""
            ),
            "terminal_stop_present": bool_text(sequence.endswith("*")),
            "transdecoder_strand": strand,
            "transdecoder_orf_start": orf_start,
            "transdecoder_orf_end": orf_end,
        }
        proteins_by_transcript[transcript_id].append(current_id)

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                store_current()
                header = line[1:]
                protein_id = header.split()[0]

                if protein_id in proteins:
                    fail(f"duplicated FASTA protein ID: {protein_id}")

                current_id = protein_id
                current_header = header
                sequence_parts = []
            else:
                if current_id is None:
                    fail(
                        f"sequence before first FASTA header at line "
                        f"{line_number}"
                    )
                sequence_parts.append(line.strip())

    store_current()

    for transcript_id in proteins_by_transcript:
        proteins_by_transcript[transcript_id].sort(key=natural_key)

    return proteins, proteins_by_transcript


def read_gene_transcript_map(
    path: Path,
) -> tuple[
    dict[str, str],
    dict[str, list[str]],
]:
    transcript_to_gene: dict[str, str] = {}
    transcripts_by_gene: dict[str, list[str]] = defaultdict(list)

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                fail(
                    f"malformed gene-transcript map line {line_number}: "
                    f"{line}"
                )

            gene_id, transcript_id = fields
            if not gene_id or not transcript_id:
                fail(
                    f"blank identifier in gene-transcript map line "
                    f"{line_number}"
                )
            if transcript_id in transcript_to_gene:
                fail(
                    f"transcript appears more than once in gene-transcript "
                    f"map: {transcript_id}"
                )

            transcript_to_gene[transcript_id] = gene_id
            transcripts_by_gene[gene_id].append(transcript_id)

    for gene_id in transcripts_by_gene:
        transcripts_by_gene[gene_id].sort(key=natural_key)

    return transcript_to_gene, transcripts_by_gene


def read_eggnog(
    path: Path,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    header: list[str] | None = None
    annotations: dict[str, dict[str, str]] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            if raw_line.startswith("##"):
                continue
            if raw_line.startswith("#"):
                header = raw_line[1:].rstrip("\n").split("\t")
                break

        if header is None:
            fail(f"eggNOG header not found: {path}")
        if len(header) != 21:
            fail(
                f"expected 21 eggNOG columns, observed {len(header)}"
            )
        if header[0] != "query":
            fail(
                f"unexpected first eggNOG column: {header[0]}"
            )

        reader = csv.DictReader(
            handle,
            fieldnames=header,
            delimiter="\t",
        )
        for line_number, row in enumerate(reader, start=2):
            if not row or not row.get("query"):
                continue
            if row["query"].startswith("#"):
                continue

            protein_id = row["query"].strip()
            if protein_id in annotations:
                fail(f"duplicated eggNOG query ID: {protein_id}")

            cleaned = {
                f"eggnog_{key.lower()}": clean_annotation(value or "")
                for key, value in row.items()
                if key != "query"
            }
            cleaned["eggnog_annotated"] = "TRUE"
            annotations[protein_id] = cleaned

    eggnog_fields = [
        "eggnog_annotated",
        *[
            f"eggnog_{field.lower()}"
            for field in header
            if field != "query"
        ],
    ]
    return eggnog_fields, annotations


def blank_eggnog(eggnog_fields: list[str]) -> dict[str, str]:
    return {
        field: ("FALSE" if field == "eggnog_annotated" else "")
        for field in eggnog_fields
    }


def read_signalp(path: Path) -> dict[str, dict[str, str]]:
    predictions: dict[str, dict[str, str]] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                fail(f"malformed SignalP row for ID: {row[0]}")

            protein_id = row[0].strip()
            prediction = row[1].strip()
            other_score = row[2].strip()
            sp_score = row[3].strip()
            cs_position = row[4].strip() if len(row) >= 5 else ""

            if protein_id in predictions:
                fail(f"duplicated SignalP ID: {protein_id}")
            if prediction not in {"OTHER", "SP"}:
                fail(
                    f"unexpected SignalP prediction for {protein_id}: "
                    f"{prediction}"
                )

            raw_sp = prediction == "SP"
            has_cleavage_site = bool(cs_position)
            confident_sp = raw_sp and has_cleavage_site

            if prediction == "OTHER":
                status = "other"
            elif confident_sp:
                status = "sp_with_cleavage_site"
            else:
                status = "ambiguous_sp_without_cleavage_site"

            predictions[protein_id] = {
                "signalp_prediction": prediction,
                "signalp_is_sp_raw": bool_text(raw_sp),
                "signalp_has_cleavage_site": bool_text(
                    has_cleavage_site
                ),
                "signalp_is_sp_confident": bool_text(confident_sp),
                "signalp_status": status,
                "signalp_other_score": other_score,
                "signalp_sp_score": sp_score,
                "signalp_cs_position": cs_position,
            }

    return predictions


def read_dbcan(path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {
            "Gene ID",
            "EC#",
            "dbCAN_hmm",
            "dbCAN_sub",
            "DIAMOND",
            "#ofTools",
            "Recommend Results",
            "Substrate",
        }

        if reader.fieldnames is None:
            fail(f"dbCAN header not found: {path}")
        missing = expected - set(reader.fieldnames)
        if missing:
            fail(f"dbCAN table missing columns: {sorted(missing)}")

        for row in reader:
            protein_id = row["Gene ID"].strip()
            if not protein_id:
                fail("blank Gene ID in dbCAN overview")
            if protein_id in annotations:
                fail(f"duplicated dbCAN Gene ID: {protein_id}")

            try:
                n_tools = int(row["#ofTools"])
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: invalid dbCAN #ofTools for "
                    f"{protein_id}: {row['#ofTools']}"
                ) from exc

            annotations[protein_id] = {
                "dbcan_any_hit": "TRUE",
                "dbcan_high_confidence": bool_text(n_tools >= 2),
                "dbcan_n_tools": str(n_tools),
                "dbcan_ec": clean_annotation(row["EC#"]),
                "dbcan_hmm": clean_annotation(row["dbCAN_hmm"]),
                "dbcan_sub": clean_annotation(row["dbCAN_sub"]),
                "dbcan_diamond": clean_annotation(row["DIAMOND"]),
                "dbcan_recommended": clean_annotation(
                    row["Recommend Results"]
                ),
                "dbcan_substrate": clean_annotation(row["Substrate"]),
            }

    return annotations


def blank_dbcan() -> dict[str, str]:
    return {
        "dbcan_any_hit": "FALSE",
        "dbcan_high_confidence": "FALSE",
        "dbcan_n_tools": "0",
        "dbcan_ec": "",
        "dbcan_hmm": "",
        "dbcan_sub": "",
        "dbcan_diamond": "",
        "dbcan_recommended": "",
        "dbcan_substrate": "",
    }


def count_true(rows: Iterable[dict[str, str]], field: str) -> int:
    return sum(row[field] == "TRUE" for row in rows)


def aggregate_annotation_fields(
    protein_rows: list[dict[str, str]],
    fields: list[str],
) -> dict[str, str]:
    return {
        field: unique_join(row.get(field, "") for row in protein_rows)
        for field in fields
    }


def main() -> int:
    args = parse_args()

    paths = {
        "proteins": args.proteins.expanduser().resolve(),
        "gene_transcript_map": args.gene_trans_map.expanduser().resolve(),
        "eggnog": args.eggnog.expanduser().resolve(),
        "signalp": args.signalp.expanduser().resolve(),
        "dbcan": args.dbcan.expanduser().resolve(),
    }
    outdir = args.outdir.expanduser().resolve()
    tables_dir = outdir / "tables"
    diagnostics_dir = outdir / "diagnostics"

    for label, path in paths.items():
        require_file(path, label)

    tables_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    print("Reading Trinity gene-transcript map...")
    transcript_to_gene, transcripts_by_gene = (
        read_gene_transcript_map(paths["gene_transcript_map"])
    )

    print("Reading TransDecoder proteins...")
    proteins, proteins_by_transcript = parse_fasta(paths["proteins"])

    protein_transcripts = set(proteins_by_transcript)
    mapped_transcripts = set(transcript_to_gene)
    mapped_genes = set(transcripts_by_gene)

    missing_protein_transcripts = sorted(
        protein_transcripts - mapped_transcripts,
        key=natural_key,
    )
    if missing_protein_transcripts:
        fail(
            f"{len(missing_protein_transcripts)} protein-bearing transcripts "
            "are absent from the Trinity gene-transcript map"
        )

    for protein_id, row in proteins.items():
        row["gene_id"] = transcript_to_gene[row["transcript_id"]]

    print("Reading eggNOG annotations...")
    eggnog_fields, eggnog = read_eggnog(paths["eggnog"])

    print("Reading SignalP predictions...")
    signalp = read_signalp(paths["signalp"])

    print("Reading dbCAN annotations...")
    dbcan = read_dbcan(paths["dbcan"])

    protein_ids = set(proteins)
    eggnog_ids = set(eggnog)
    signalp_ids = set(signalp)
    dbcan_ids = set(dbcan)

    if signalp_ids != protein_ids:
        fail(
            "SignalP IDs do not exactly match the TransDecoder proteome; "
            f"missing={len(protein_ids - signalp_ids)}, "
            f"unexpected={len(signalp_ids - protein_ids)}"
        )

    unexpected_eggnog_ids = sorted(
        eggnog_ids - protein_ids,
        key=natural_key,
    )
    unexpected_dbcan_ids = sorted(
        dbcan_ids - protein_ids,
        key=natural_key,
    )

    if unexpected_eggnog_ids:
        fail(
            f"eggNOG contains {len(unexpected_eggnog_ids)} IDs absent "
            "from the TransDecoder proteome"
        )
    if unexpected_dbcan_ids:
        fail(
            f"dbCAN contains {len(unexpected_dbcan_ids)} IDs absent "
            "from the TransDecoder proteome"
        )

    protein_rows: list[dict[str, str]] = []
    blank_eggnog_row = blank_eggnog(eggnog_fields)

    for protein_id in sorted(proteins, key=natural_key):
        row = dict(proteins[protein_id])
        row.update(eggnog.get(protein_id, blank_eggnog_row))
        row.update(signalp[protein_id])
        row.update(dbcan.get(protein_id, blank_dbcan()))
        protein_rows.append(row)

    structural_protein_fields = [
        "gene_id",
        "transcript_id",
        "protein_id",
        "transdecoder_orf_number",
        "transdecoder_orf_type",
        "protein_length_aa",
        "fasta_sequence_length",
        "transdecoder_reported_length",
        "transdecoder_length_difference",
        "terminal_stop_present",
        "transdecoder_strand",
        "transdecoder_orf_start",
        "transdecoder_orf_end",
    ]
    signalp_fields = [
        "signalp_prediction",
        "signalp_is_sp_raw",
        "signalp_has_cleavage_site",
        "signalp_is_sp_confident",
        "signalp_status",
        "signalp_other_score",
        "signalp_sp_score",
        "signalp_cs_position",
    ]
    dbcan_fields = list(blank_dbcan())
    protein_fields = (
        structural_protein_fields
        + eggnog_fields
        + signalp_fields
        + dbcan_fields
    )

    protein_table = (
        tables_dir / "fusarium_protein_functional_annotation.tsv"
    )
    write_tsv(protein_table, protein_fields, protein_rows)

    protein_rows_by_transcript: dict[str, list[dict[str, str]]] = (
        defaultdict(list)
    )
    protein_rows_by_gene: dict[str, list[dict[str, str]]] = (
        defaultdict(list)
    )
    for row in protein_rows:
        protein_rows_by_transcript[row["transcript_id"]].append(row)
        protein_rows_by_gene[row["gene_id"]].append(row)

    eggnog_text_fields = [
        field
        for field in eggnog_fields
        if field != "eggnog_annotated"
    ]
    dbcan_text_fields = [
        "dbcan_ec",
        "dbcan_hmm",
        "dbcan_sub",
        "dbcan_diamond",
        "dbcan_recommended",
        "dbcan_substrate",
    ]

    transcript_rows: list[dict[str, str]] = []
    for transcript_id in sorted(mapped_transcripts, key=natural_key):
        gene_id = transcript_to_gene[transcript_id]
        associated = protein_rows_by_transcript.get(transcript_id, [])

        row = {
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "has_predicted_protein": bool_text(bool(associated)),
            "predicted_protein_count": str(len(associated)),
            "protein_ids": unique_join(
                protein["protein_id"] for protein in associated
            ),
            "eggnog_annotated_protein_count": str(
                count_true(associated, "eggnog_annotated")
            ),
            "eggnog_annotated": bool_text(
                any(
                    protein["eggnog_annotated"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_raw_positive_protein_count": str(
                count_true(associated, "signalp_is_sp_raw")
            ),
            "signalp_raw_positive": bool_text(
                any(
                    protein["signalp_is_sp_raw"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_confident_positive_protein_count": str(
                count_true(associated, "signalp_is_sp_confident")
            ),
            "signalp_confident_positive": bool_text(
                any(
                    protein["signalp_is_sp_confident"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_ambiguous_protein_count": str(
                sum(
                    protein["signalp_status"]
                    == "ambiguous_sp_without_cleavage_site"
                    for protein in associated
                )
            ),
            "dbcan_any_hit_protein_count": str(
                count_true(associated, "dbcan_any_hit")
            ),
            "dbcan_any_hit": bool_text(
                any(
                    protein["dbcan_any_hit"] == "TRUE"
                    for protein in associated
                )
            ),
            "dbcan_high_confidence_protein_count": str(
                count_true(associated, "dbcan_high_confidence")
            ),
            "dbcan_high_confidence": bool_text(
                any(
                    protein["dbcan_high_confidence"] == "TRUE"
                    for protein in associated
                )
            ),
            "dbcan_max_n_tools": str(
                max(
                    (
                        int(protein["dbcan_n_tools"])
                        for protein in associated
                    ),
                    default=0,
                )
            ),
        }
        row.update(
            aggregate_annotation_fields(
                associated,
                eggnog_text_fields + dbcan_text_fields,
            )
        )
        transcript_rows.append(row)

    aggregate_base_fields = [
        "gene_id",
        "transcript_id",
        "has_predicted_protein",
        "predicted_protein_count",
        "protein_ids",
        "eggnog_annotated_protein_count",
        "eggnog_annotated",
        "signalp_raw_positive_protein_count",
        "signalp_raw_positive",
        "signalp_confident_positive_protein_count",
        "signalp_confident_positive",
        "signalp_ambiguous_protein_count",
        "dbcan_any_hit_protein_count",
        "dbcan_any_hit",
        "dbcan_high_confidence_protein_count",
        "dbcan_high_confidence",
        "dbcan_max_n_tools",
    ]
    aggregate_annotation_fields_order = (
        eggnog_text_fields + dbcan_text_fields
    )
    transcript_fields = (
        aggregate_base_fields + aggregate_annotation_fields_order
    )

    transcript_table = (
        tables_dir / "fusarium_transcript_functional_annotation.tsv"
    )
    write_tsv(transcript_table, transcript_fields, transcript_rows)

    gene_rows: list[dict[str, str]] = []
    for gene_id in sorted(mapped_genes, key=natural_key):
        all_transcripts = transcripts_by_gene[gene_id]
        associated = protein_rows_by_gene.get(gene_id, [])
        protein_coding_transcripts = sorted(
            {
                protein["transcript_id"]
                for protein in associated
            },
            key=natural_key,
        )

        row = {
            "gene_id": gene_id,
            "mapped_transcript_count": str(len(all_transcripts)),
            "transcript_ids": unique_join(all_transcripts),
            "protein_coding_transcript_count": str(
                len(protein_coding_transcripts)
            ),
            "protein_coding_transcript_ids": unique_join(
                protein_coding_transcripts
            ),
            "has_predicted_protein": bool_text(bool(associated)),
            "predicted_protein_count": str(len(associated)),
            "protein_ids": unique_join(
                protein["protein_id"] for protein in associated
            ),
            "eggnog_annotated_protein_count": str(
                count_true(associated, "eggnog_annotated")
            ),
            "eggnog_annotated": bool_text(
                any(
                    protein["eggnog_annotated"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_raw_positive_protein_count": str(
                count_true(associated, "signalp_is_sp_raw")
            ),
            "signalp_raw_positive": bool_text(
                any(
                    protein["signalp_is_sp_raw"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_confident_positive_protein_count": str(
                count_true(associated, "signalp_is_sp_confident")
            ),
            "signalp_confident_positive": bool_text(
                any(
                    protein["signalp_is_sp_confident"] == "TRUE"
                    for protein in associated
                )
            ),
            "signalp_ambiguous_protein_count": str(
                sum(
                    protein["signalp_status"]
                    == "ambiguous_sp_without_cleavage_site"
                    for protein in associated
                )
            ),
            "dbcan_any_hit_protein_count": str(
                count_true(associated, "dbcan_any_hit")
            ),
            "dbcan_any_hit": bool_text(
                any(
                    protein["dbcan_any_hit"] == "TRUE"
                    for protein in associated
                )
            ),
            "dbcan_high_confidence_protein_count": str(
                count_true(associated, "dbcan_high_confidence")
            ),
            "dbcan_high_confidence": bool_text(
                any(
                    protein["dbcan_high_confidence"] == "TRUE"
                    for protein in associated
                )
            ),
            "dbcan_max_n_tools": str(
                max(
                    (
                        int(protein["dbcan_n_tools"])
                        for protein in associated
                    ),
                    default=0,
                )
            ),
        }
        row.update(
            aggregate_annotation_fields(
                associated,
                eggnog_text_fields + dbcan_text_fields,
            )
        )
        gene_rows.append(row)

    gene_base_fields = [
        "gene_id",
        "mapped_transcript_count",
        "transcript_ids",
        "protein_coding_transcript_count",
        "protein_coding_transcript_ids",
        "has_predicted_protein",
        "predicted_protein_count",
        "protein_ids",
        "eggnog_annotated_protein_count",
        "eggnog_annotated",
        "signalp_raw_positive_protein_count",
        "signalp_raw_positive",
        "signalp_confident_positive_protein_count",
        "signalp_confident_positive",
        "signalp_ambiguous_protein_count",
        "dbcan_any_hit_protein_count",
        "dbcan_any_hit",
        "dbcan_high_confidence_protein_count",
        "dbcan_high_confidence",
        "dbcan_max_n_tools",
    ]
    gene_fields = (
        gene_base_fields + aggregate_annotation_fields_order
    )

    gene_table = (
        tables_dir / "fusarium_gene_functional_annotation.tsv"
    )
    write_tsv(gene_table, gene_fields, gene_rows)

    ambiguous_rows = [
        {
            "gene_id": row["gene_id"],
            "transcript_id": row["transcript_id"],
            "protein_id": row["protein_id"],
            "signalp_prediction": row["signalp_prediction"],
            "signalp_other_score": row["signalp_other_score"],
            "signalp_sp_score": row["signalp_sp_score"],
            "signalp_cs_position": row["signalp_cs_position"],
            "recommended_interpretation": (
                "raw SignalP-positive; exclude from conservative "
                "cleavage-site-supported secretome count"
            ),
        }
        for row in protein_rows
        if row["signalp_status"]
        == "ambiguous_sp_without_cleavage_site"
    ]
    write_tsv(
        diagnostics_dir
        / "signalp_ambiguous_no_cleavage_site.tsv",
        [
            "gene_id",
            "transcript_id",
            "protein_id",
            "signalp_prediction",
            "signalp_other_score",
            "signalp_sp_score",
            "signalp_cs_position",
            "recommended_interpretation",
        ],
        ambiguous_rows,
    )

    write_tsv(
        diagnostics_dir / "eggnog_ids_not_in_proteome.tsv",
        ["protein_id"],
        [
            {"protein_id": protein_id}
            for protein_id in unexpected_eggnog_ids
        ],
    )
    write_tsv(
        diagnostics_dir / "dbcan_ids_not_in_proteome.tsv",
        ["protein_id"],
        [
            {"protein_id": protein_id}
            for protein_id in unexpected_dbcan_ids
        ],
    )

    proteins_per_transcript = Counter(
        len(proteins_by_transcript.get(transcript_id, []))
        for transcript_id in mapped_transcripts
    )
    proteins_per_gene = Counter(
        len(protein_rows_by_gene.get(gene_id, []))
        for gene_id in mapped_genes
    )
    transcripts_per_gene = Counter(
        len(transcripts_by_gene[gene_id])
        for gene_id in mapped_genes
    )

    multiplicity_rows: list[dict[str, object]] = []
    for level, distribution in (
        ("proteins_per_transcript", proteins_per_transcript),
        ("proteins_per_gene", proteins_per_gene),
        ("transcripts_per_gene", transcripts_per_gene),
    ):
        for multiplicity in sorted(distribution):
            multiplicity_rows.append(
                {
                    "level": level,
                    "multiplicity": multiplicity,
                    "entity_count": distribution[multiplicity],
                }
            )

    write_tsv(
        diagnostics_dir / "mapping_multiplicity.tsv",
        ["level", "multiplicity", "entity_count"],
        multiplicity_rows,
    )

    protein_count = len(protein_rows)
    transcript_count = len(transcript_rows)
    gene_count = len(gene_rows)
    transcripts_with_protein = sum(
        row["has_predicted_protein"] == "TRUE"
        for row in transcript_rows
    )
    genes_with_protein = sum(
        row["has_predicted_protein"] == "TRUE"
        for row in gene_rows
    )
    eggnog_count = count_true(protein_rows, "eggnog_annotated")
    dbcan_count = count_true(protein_rows, "dbcan_any_hit")
    signalp_raw_count = count_true(
        protein_rows,
        "signalp_is_sp_raw",
    )
    signalp_confident_count = count_true(
        protein_rows,
        "signalp_is_sp_confident",
    )
    signalp_ambiguous_count = len(ambiguous_rows)

    validations = [
        (
            "protein_rows",
            protein_count,
            EXPECTED_PROTEINS,
        ),
        (
            "unique_protein_ids",
            len(protein_ids),
            EXPECTED_PROTEINS,
        ),
        (
            "transcript_rows",
            transcript_count,
            EXPECTED_TRANSCRIPTS,
        ),
        (
            "unique_transcript_ids",
            len(mapped_transcripts),
            EXPECTED_TRANSCRIPTS,
        ),
        (
            "gene_rows",
            gene_count,
            EXPECTED_GENES,
        ),
        (
            "unique_gene_ids",
            len(mapped_genes),
            EXPECTED_GENES,
        ),
        (
            "transcripts_with_predicted_protein",
            transcripts_with_protein,
            EXPECTED_TRANSCRIPTS_WITH_PROTEIN,
        ),
        (
            "genes_with_predicted_protein",
            genes_with_protein,
            EXPECTED_GENES_WITH_PROTEIN,
        ),
        (
            "protein_transcripts_absent_from_map",
            len(missing_protein_transcripts),
            0,
        ),
        (
            "signalp_rows",
            len(signalp),
            EXPECTED_SIGNALP_ROWS,
        ),
        (
            "signalp_raw_sp",
            signalp_raw_count,
            EXPECTED_SIGNALP_RAW_SP,
        ),
        (
            "signalp_confident_sp",
            signalp_confident_count,
            EXPECTED_SIGNALP_CONFIDENT_SP,
        ),
        (
            "signalp_ambiguous_sp_without_cleavage_site",
            signalp_ambiguous_count,
            EXPECTED_SIGNALP_AMBIGUOUS_SP,
        ),
        (
            "eggnog_annotated_proteins",
            eggnog_count,
            EXPECTED_EGGNOG_ANNOTATED_PROTEINS,
        ),
        (
            "dbcan_overview_proteins",
            dbcan_count,
            EXPECTED_DBCAN_PROTEINS,
        ),
        (
            "eggnog_ids_not_in_proteome",
            len(unexpected_eggnog_ids),
            0,
        ),
        (
            "dbcan_ids_not_in_proteome",
            len(unexpected_dbcan_ids),
            0,
        ),
    ]

    validation_rows = [
        {
            "check": check,
            "value": value,
            "expected": expected,
            "status": "PASS" if value == expected else "FAIL",
        }
        for check, value, expected in validations
    ]
    write_tsv(
        diagnostics_dir / "validation.tsv",
        ["check", "value", "expected", "status"],
        validation_rows,
    )

    failures = [
        row for row in validation_rows if row["status"] != "PASS"
    ]
    if failures:
        failed_checks = ", ".join(
            str(row["check"]) for row in failures
        )
        fail(f"validation failed: {failed_checks}")

    run_info_rows: list[dict[str, object]] = [
        {
            "field": "script",
            "value": "18_build_fusarium_functional_annotation.py",
        },
        {
            "field": "run_timestamp_utc",
            "value": datetime.now(timezone.utc).isoformat(),
        },
        {
            "field": "repository_root",
            "value": str(Path.cwd().resolve()),
        },
        {
            "field": "protein_annotation_scope",
            "value": (
                "all 35327 TransDecoder predicted proteins"
            ),
        },
        {
            "field": "transcript_annotation_scope",
            "value": (
                "all 27979 Trinity transcripts, including 8172 "
                "without a predicted protein"
            ),
        },
        {
            "field": "gene_annotation_scope",
            "value": (
                "all 15192 Trinity genes, including 5995 "
                "without a predicted protein"
            ),
        },
        {
            "field": "signalp_raw_positive_definition",
            "value": "Prediction == SP",
        },
        {
            "field": "signalp_confident_positive_definition",
            "value": (
                "Prediction == SP and cleavage-site position present"
            ),
        },
        {
            "field": "dbcan_any_hit_definition",
            "value": "protein present in dbCAN overview.tsv",
        },
        {
            "field": "dbcan_high_confidence_definition",
            "value": "dbCAN #ofTools >= 2",
        },
        {
            "field": "gene_level_collapse_definition",
            "value": (
                "Boolean presence if at least one associated predicted "
                "protein satisfies the criterion; multiplicity retained "
                "in protein-count fields"
            ),
        },
    ]

    for label, path in paths.items():
        run_info_rows.extend(
            [
                {
                    "field": f"{label}_file",
                    "value": str(path),
                },
                {
                    "field": f"{label}_sha256",
                    "value": sha256(path),
                },
            ]
        )

    run_info_rows.extend(
        [
            {
                "field": "protein_output",
                "value": str(protein_table),
            },
            {
                "field": "transcript_output",
                "value": str(transcript_table),
            },
            {
                "field": "gene_output",
                "value": str(gene_table),
            },
        ]
    )

    write_tsv(
        outdir / "run_info.tsv",
        ["field", "value"],
        run_info_rows,
    )

    checksum_path = outdir / "checksums.sha256"
    output_files = sorted(
        (
            path
            for path in outdir.rglob("*")
            if path.is_file() and path != checksum_path
        ),
        key=lambda path: str(path.relative_to(outdir)),
    )
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in output_files:
            relative = path.relative_to(outdir)
            handle.write(f"{sha256(path)}  {relative}\n")

    print()
    print("Fusarium functional annotation build")
    print()
    print(f"Proteins:                         {protein_count:,}")
    print(f"Trinity transcripts:             {transcript_count:,}")
    print(
        f"Transcripts with protein:         "
        f"{transcripts_with_protein:,}"
    )
    print(f"Trinity genes:                   {gene_count:,}")
    print(f"Genes with protein:               {genes_with_protein:,}")
    print(f"eggNOG-annotated proteins:        {eggnog_count:,}")
    print(f"dbCAN overview proteins:          {dbcan_count:,}")
    print(f"SignalP raw positives:            {signalp_raw_count:,}")
    print(
        f"SignalP confident positives:      "
        f"{signalp_confident_count:,}"
    )
    print(
        f"SignalP ambiguous positives:      "
        f"{signalp_ambiguous_count:,}"
    )
    print()
    print("Validation status: PASS")
    print(f"Output directory: {outdir}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        sys.exit("Interrupted")
