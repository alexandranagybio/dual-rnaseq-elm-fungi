#!/usr/bin/env python3
"""
05_prepare_deseq2_matrix.py

Convert the validated featureCounts output into:

1. A clean 8,640-gene integer count matrix.
2. A canonical sample metadata table.
3. A machine-readable validation report.

The script intentionally stops on:
- missing samples,
- duplicated gene IDs,
- non-integer or negative counts,
- unexpected gene-row count,
- unexpected sample columns,
- malformed featureCounts output.

Usage:
    python3 workflow/scripts/05_prepare_deseq2_matrix.py

Optional overrides:
    python3 workflow/scripts/05_prepare_deseq2_matrix.py \
        --input results/ophiostoma/gene_counts/ophiostoma_gene_featurecounts.txt \
        --outdir results/ophiostoma/gene_counts
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


EXPECTED_GENES = 8640

SAMPLE_METADATA: Sequence[Tuple[str, str, str]] = (
    ("149", "interaction", "1"),
    ("150", "interaction", "2"),
    ("151", "interaction", "3"),
    ("152", "self", "1"),
    ("153", "self", "2"),
    ("154", "self", "3"),
    ("155", "onu", "1"),
    ("156", "onu", "2"),
    ("157", "onu", "3"),
)


class ValidationError(RuntimeError):
    """Raised when an input or output violates an audit invariant."""


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent

    parser = argparse.ArgumentParser(
        description="Prepare and validate the Ophiostoma gene-level DESeq2 matrix."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=repo_root
        / "results"
        / "ophiostoma"
        / "gene_counts"
        / "ophiostoma_gene_featurecounts.txt",
        help="featureCounts output table.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=repo_root / "results" / "ophiostoma" / "gene_counts",
        help="Output directory.",
    )
    parser.add_argument(
        "--expected-genes",
        type=int,
        default=EXPECTED_GENES,
        help="Expected number of gene rows.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_id_from_column(column_name: str) -> str:
    """Extract sample 149-157 from an absolute or relative BAM path."""
    basename = Path(column_name).name
    match = re.fullmatch(r"(149|150|151|152|153|154|155|156|157)_ophiostoma\.sorted\.bam", basename)
    if match is None:
        raise ValidationError(
            f"Could not map featureCounts column to a canonical sample ID: {column_name!r}"
        )
    return match.group(1)


def read_featurecounts(
    path: Path,
) -> Tuple[List[str], List[List[int]], List[str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValidationError(f"Input is missing or empty: {path}")

    header: List[str] | None = None
    rows: List[List[int]] = []
    gene_ids: List[str] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t",
        )

        for record_number, record in enumerate(reader, start=1):
            if header is None:
                header = record
                continue

            if len(record) != len(header):
                raise ValidationError(
                    f"Row {record_number} has {len(record)} fields; expected {len(header)}."
                )

            gene_id = record[0].strip()
            if not gene_id:
                raise ValidationError(f"Row {record_number} has an empty Geneid.")

            try:
                counts = [int(value) for value in record[6:]]
            except ValueError as exc:
                raise ValidationError(
                    f"Non-integer count at row {record_number}, gene {gene_id!r}."
                ) from exc

            if any(value < 0 for value in counts):
                raise ValidationError(
                    f"Negative count at row {record_number}, gene {gene_id!r}."
                )

            gene_ids.append(gene_id)
            rows.append(counts)

    if header is None:
        raise ValidationError("No featureCounts header was found.")

    expected_fixed = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    if header[:6] != expected_fixed:
        raise ValidationError(
            "Unexpected featureCounts fixed columns. "
            f"Observed {header[:6]!r}; expected {expected_fixed!r}."
        )

    raw_sample_columns = header[6:]
    sample_ids = [sample_id_from_column(column) for column in raw_sample_columns]

    if len(set(sample_ids)) != len(sample_ids):
        raise ValidationError(f"Duplicated sample columns: {sample_ids!r}")

    expected_samples = [row[0] for row in SAMPLE_METADATA]
    if set(sample_ids) != set(expected_samples):
        raise ValidationError(
            f"Unexpected sample set. Observed {sample_ids!r}; expected {expected_samples!r}."
        )

    if len(set(gene_ids)) != len(gene_ids):
        seen = set()
        duplicates = []
        for gene_id in gene_ids:
            if gene_id in seen:
                duplicates.append(gene_id)
            seen.add(gene_id)
        raise ValidationError(
            "Duplicated gene IDs detected: " + ", ".join(sorted(set(duplicates))[:20])
        )

    # Reorder count columns into canonical numeric sample order.
    position = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    reordered_rows = [
        [counts[position[sample_id]] for sample_id in expected_samples]
        for counts in rows
    ]

    return gene_ids, reordered_rows, expected_samples


def write_matrix(
    path: Path,
    gene_ids: Sequence[str],
    counts: Sequence[Sequence[int]],
    sample_ids: Sequence[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_ids])
        for gene_id, row in zip(gene_ids, counts):
            writer.writerow([gene_id, *row])


def write_metadata(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "condition", "replicate"])
        writer.writerows(SAMPLE_METADATA)


def column_sums(
    rows: Sequence[Sequence[int]], sample_ids: Sequence[str]
) -> Dict[str, int]:
    return {
        sample_id: sum(row[index] for row in rows)
        for index, sample_id in enumerate(sample_ids)
    }


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_path = outdir / "ophiostoma_gene_counts_matrix.tsv"
    metadata_path = outdir / "sample_metadata.tsv"
    report_path = outdir / "05_prepare_deseq2_matrix.validation.json"
    checksum_path = outdir / "05_prepare_deseq2_matrix.sha256"

    try:
        gene_ids, counts, sample_ids = read_featurecounts(input_path)

        if len(gene_ids) != args.expected_genes:
            raise ValidationError(
                f"Expected {args.expected_genes} genes, observed {len(gene_ids)}."
            )

        write_matrix(matrix_path, gene_ids, counts, sample_ids)
        write_metadata(metadata_path)

        sums = column_sums(counts, sample_ids)
        zero_total_genes = sum(1 for row in counts if sum(row) == 0)

        report = {
            "status": "PASS",
            "input_featurecounts": str(input_path),
            "input_sha256": sha256(input_path),
            "expected_gene_rows": args.expected_genes,
            "observed_gene_rows": len(gene_ids),
            "sample_ids": sample_ids,
            "sample_metadata": [
                {
                    "sample_id": sample_id,
                    "condition": condition,
                    "replicate": replicate,
                }
                for sample_id, condition, replicate in SAMPLE_METADATA
            ],
            "library_sizes": sums,
            "zero_total_gene_rows": zero_total_genes,
            "matrix": str(matrix_path),
            "metadata": str(metadata_path),
        }

        with report_path.open("w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")

        with checksum_path.open("w", encoding="utf-8") as handle:
            for path in (input_path, matrix_path, metadata_path, report_path):
                handle.write(f"{sha256(path)}  {path}\n")

    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("PASS: Ophiostoma gene-level DESeq2 inputs were created.")
    print(f"Genes: {len(gene_ids):,}")
    print(f"Samples: {', '.join(sample_ids)}")
    print(f"Matrix: {matrix_path}")
    print(f"Metadata: {metadata_path}")
    print(f"Validation report: {report_path}")
    print(f"Checksums: {checksum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
