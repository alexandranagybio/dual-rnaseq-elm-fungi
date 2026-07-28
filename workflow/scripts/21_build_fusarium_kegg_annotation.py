#!/usr/bin/env python3

"""
Build a validated gene-level KEGG pathway annotation for Fusarium cf. salinense.

This script consumes the validated functional-annotation tables produced by
18_build_fusarium_functional_annotation.py. It does not reparse the raw eggNOG
file or reconstruct the Trinity protein-to-gene mapping.

Input
-----
1. Validated protein-level functional annotation table.
2. Validated gene-level functional annotation table.

KEGG pathway normalization
--------------------------
eggNOG-mapper reports both `ko#####` and `map#####` forms for the same generic
KEGG pathway. These are collapsed to one canonical identifier:

    ko00010  -> map00010
    map00010 -> map00010

Only tokens matching `ko` or `map` followed by exactly five digits are usable.

Outputs
-------
- Protein-level canonical KEGG pathway annotation
- Gene-level canonical KEGG pathway annotation
- TERM2GENE table for clusterProfiler::enricher()
- Canonical pathway metadata
- Raw-to-canonical pathway normalization table
- Genes without a usable KEGG pathway
- Invalid-token diagnostic table
- Validation JSON
- SHA256 checksum file
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_PROTEIN_ANNOTATION = Path(
    "results/fusarium/functional_annotation/tables/"
    "fusarium_protein_functional_annotation.tsv"
)

DEFAULT_GENE_ANNOTATION = Path(
    "results/fusarium/functional_annotation/tables/"
    "fusarium_gene_functional_annotation.tsv"
)

DEFAULT_OUT_DIR = Path("results/fusarium/kegg_annotation")

MISSING_TOKENS = {"", "-", r"\N", "NA", "NaN", "nan", "None"}
PATHWAY_PATTERN = re.compile(r"^(ko|map)(\d{5})$")

EXPECTED_PROTEIN_ROWS = 35_327
EXPECTED_UNIQUE_PROTEINS = 35_327
EXPECTED_PROTEIN_TABLE_GENES = 9_197
EXPECTED_STRUCTURAL_GENES = 15_192
EXPECTED_PROTEINS_WITH_PATHWAY_FIELD = 8_282
EXPECTED_UNIQUE_RAW_PATHWAY_STRINGS = 914
EXPECTED_UNIQUE_RAW_PATHWAY_TOKENS = 766
EXPECTED_CANONICAL_PATHWAYS = 383
EXPECTED_PROTEIN_PATHWAY_PAIRS = 33_884
EXPECTED_GENE_PATHWAY_PAIRS = 10_897
EXPECTED_GENES_WITH_PATHWAY = 2_503
EXPECTED_INVALID_PATHWAY_TOKENS = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build validated canonical KEGG pathway annotations from the "
            "Fusarium functional-annotation tables."
        )
    )
    parser.add_argument(
        "--protein-annotation",
        type=Path,
        default=DEFAULT_PROTEIN_ANNOTATION,
        help="Validated protein-level functional annotation TSV.",
    )
    parser.add_argument(
        "--gene-annotation",
        type=Path,
        default=DEFAULT_GENE_ANNOTATION,
        help="Validated gene-level functional annotation TSV.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory.",
    )
    parser.add_argument(
        "--skip-expected-counts",
        action="store_true",
        help=(
            "Skip dataset-specific expected-count checks. Structural and "
            "relational validation checks are still enforced."
        ),
    )
    return parser.parse_args()


def stop(message: str) -> None:
    raise RuntimeError(message)


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_csv_field_limit() -> None:
    limit = 2**31 - 1
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def read_header(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration:
            stop(f"Input file is empty: {path}")


def require_columns(
    path: Path,
    observed: Iterable[str],
    required: set[str],
) -> None:
    missing = sorted(required - set(observed))
    if missing:
        stop(f"{path} is missing required columns: " + ", ".join(missing))


def write_tsv(
    rows: Iterable[dict[str, object]],
    path: Path,
    fieldnames: list[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def canonicalize_pathway(raw_token: str) -> tuple[str, str, str] | None:
    match = PATHWAY_PATTERN.fullmatch(raw_token)
    if match is None:
        return None
    prefix, numeric_id = match.groups()
    return prefix, numeric_id, f"map{numeric_id}"


def main() -> None:
    args = parse_args()
    configure_csv_field_limit()

    protein_path = args.protein_annotation.expanduser().resolve()
    gene_path = args.gene_annotation.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    diagnostics_dir = out_dir / "diagnostics"

    require_file(protein_path, "Protein annotation table")
    require_file(gene_path, "Gene annotation table")

    out_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    protein_header = read_header(protein_path)
    gene_header = read_header(gene_path)

    require_columns(
        protein_path,
        protein_header,
        {"gene_id", "protein_id", "eggnog_kegg_pathway"},
    )
    require_columns(gene_path, gene_header, {"gene_id"})

    all_gene_ids: list[str] = []
    seen_gene_ids: set[str] = set()

    with gene_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=2):
            gene_id = row["gene_id"].strip()
            if not gene_id:
                stop(
                    f"Blank gene_id in gene annotation table at line "
                    f"{line_number}."
                )
            if gene_id in seen_gene_ids:
                stop(f"Duplicated gene_id in gene annotation table: {gene_id}")
            seen_gene_ids.add(gene_id)
            all_gene_ids.append(gene_id)

    protein_rows = 0
    protein_ids: set[str] = set()
    protein_table_genes: set[str] = set()
    proteins_with_pathway_field = 0
    raw_pathway_strings: set[str] = set()
    raw_token_counts: Counter[str] = Counter()
    invalid_token_counts: Counter[str] = Counter()
    normalization_counts: Counter[tuple[str, str]] = Counter()

    protein_pathway_pairs: set[tuple[str, str, str]] = set()
    gene_pathway_pairs: set[tuple[str, str]] = set()
    pathway_to_prefixes: defaultdict[str, set[str]] = defaultdict(set)

    with protein_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=2):
            protein_rows += 1
            protein_id = row["protein_id"].strip()
            gene_id = row["gene_id"].strip()
            raw_value = row["eggnog_kegg_pathway"].strip()

            if not protein_id:
                stop(
                    f"Blank protein_id in protein annotation table at line "
                    f"{line_number}."
                )
            if protein_id in protein_ids:
                stop(
                    f"Duplicated protein_id in protein annotation table: "
                    f"{protein_id}"
                )
            if not gene_id:
                stop(
                    f"Blank gene_id in protein annotation table at line "
                    f"{line_number}."
                )
            if gene_id not in seen_gene_ids:
                stop(
                    f"Protein-table gene absent from validated gene table: "
                    f"{gene_id}"
                )

            protein_ids.add(protein_id)
            protein_table_genes.add(gene_id)

            if raw_value in MISSING_TOKENS:
                continue

            proteins_with_pathway_field += 1
            raw_pathway_strings.add(raw_value)

            for token in raw_value.split(","):
                token = token.strip()
                if token in MISSING_TOKENS:
                    continue

                raw_token_counts[token] += 1
                parsed = canonicalize_pathway(token)
                if parsed is None:
                    invalid_token_counts[token] += 1
                    continue

                prefix, _numeric_id, canonical = parsed
                normalization_counts[(token, canonical)] += 1
                pathway_to_prefixes[canonical].add(prefix)
                protein_pathway_pairs.add((protein_id, gene_id, canonical))
                gene_pathway_pairs.add((gene_id, canonical))

    observed_counts = {
        "protein_rows": protein_rows,
        "unique_proteins": len(protein_ids),
        "unique_genes_in_protein_table": len(protein_table_genes),
        "structural_genes": len(all_gene_ids),
        "proteins_with_nonempty_kegg_pathway_field": proteins_with_pathway_field,
        "unique_raw_kegg_pathway_strings": len(raw_pathway_strings),
        "unique_raw_kegg_pathway_tokens": len(raw_token_counts),
        "canonical_kegg_pathways": len(pathway_to_prefixes),
        "unique_protein_pathway_pairs": len(protein_pathway_pairs),
        "unique_gene_pathway_pairs": len(gene_pathway_pairs),
        "genes_with_canonical_kegg_pathway": len(
            {gene_id for gene_id, _ in gene_pathway_pairs}
        ),
        "invalid_kegg_pathway_tokens": len(invalid_token_counts),
    }

    expected_counts = {
        "protein_rows": EXPECTED_PROTEIN_ROWS,
        "unique_proteins": EXPECTED_UNIQUE_PROTEINS,
        "unique_genes_in_protein_table": EXPECTED_PROTEIN_TABLE_GENES,
        "structural_genes": EXPECTED_STRUCTURAL_GENES,
        "proteins_with_nonempty_kegg_pathway_field": EXPECTED_PROTEINS_WITH_PATHWAY_FIELD,
        "unique_raw_kegg_pathway_strings": EXPECTED_UNIQUE_RAW_PATHWAY_STRINGS,
        "unique_raw_kegg_pathway_tokens": EXPECTED_UNIQUE_RAW_PATHWAY_TOKENS,
        "canonical_kegg_pathways": EXPECTED_CANONICAL_PATHWAYS,
        "unique_protein_pathway_pairs": EXPECTED_PROTEIN_PATHWAY_PAIRS,
        "unique_gene_pathway_pairs": EXPECTED_GENE_PATHWAY_PAIRS,
        "genes_with_canonical_kegg_pathway": EXPECTED_GENES_WITH_PATHWAY,
        "invalid_kegg_pathway_tokens": EXPECTED_INVALID_PATHWAY_TOKENS,
    }

    if not args.skip_expected_counts:
        for label, expected in expected_counts.items():
            observed = observed_counts[label]
            if observed != expected:
                stop(
                    f"Expected {expected:,} for {label}, "
                    f"observed {observed:,}."
                )

    if invalid_token_counts:
        examples = ", ".join(
            token for token, _ in invalid_token_counts.most_common(20)
        )
        stop(
            f"Observed {len(invalid_token_counts):,} invalid KEGG pathway "
            f"tokens. Examples: {examples}"
        )

    incomplete_prefixes = {
        pathway: prefixes
        for pathway, prefixes in pathway_to_prefixes.items()
        if prefixes != {"ko", "map"}
    }
    if incomplete_prefixes:
        preview = ", ".join(
            f"{pathway}={','.join(sorted(prefixes))}"
            for pathway, prefixes in list(sorted(incomplete_prefixes.items()))[:20]
        )
        stop(
            "Some canonical pathways were not represented by both ko and "
            f"map prefixes: {preview}"
        )

    protein_output = out_dir / "fusarium_protein_kegg_pathways.tsv"
    gene_output = out_dir / "fusarium_gene_kegg_pathways.tsv"
    term2gene_output = out_dir / "fusarium_kegg_term2gene.tsv"
    metadata_output = out_dir / "fusarium_kegg_pathway_metadata.tsv"
    normalization_output = out_dir / "pathway_normalization.tsv"
    missing_genes_output = out_dir / "genes_without_kegg_pathway.tsv"
    invalid_output = diagnostics_dir / "invalid_kegg_pathway_tokens.tsv"
    validation_output = out_dir / "21_build_fusarium_kegg_annotation.validation.json"
    checksum_output = out_dir / "21_build_fusarium_kegg_annotation.sha256"

    write_tsv(
        [
            {"protein_id": protein_id, "gene_id": gene_id, "pathway": pathway}
            for protein_id, gene_id, pathway in sorted(
                protein_pathway_pairs,
                key=lambda x: (x[2], x[1], x[0]),
            )
        ],
        protein_output,
        ["protein_id", "gene_id", "pathway"],
    )

    write_tsv(
        [
            {"gene_id": gene_id, "pathway": pathway}
            for gene_id, pathway in sorted(
                gene_pathway_pairs,
                key=lambda x: (x[1], x[0]),
            )
        ],
        gene_output,
        ["gene_id", "pathway"],
    )

    write_tsv(
        [
            {"pathway": pathway, "gene_id": gene_id}
            for gene_id, pathway in sorted(
                gene_pathway_pairs,
                key=lambda x: (x[1], x[0]),
            )
        ],
        term2gene_output,
        ["pathway", "gene_id"],
    )

    write_tsv(
        [
            {
                "pathway": pathway,
                "pathway_numeric_id": pathway.removeprefix("map"),
                "pathway_name": pathway,
                "raw_prefixes_observed": ",".join(
                    sorted(pathway_to_prefixes[pathway])
                ),
            }
            for pathway in sorted(pathway_to_prefixes)
        ],
        metadata_output,
        [
            "pathway",
            "pathway_numeric_id",
            "pathway_name",
            "raw_prefixes_observed",
        ],
    )

    write_tsv(
        [
            {
                "raw_pathway": raw_pathway,
                "canonical_pathway": canonical_pathway,
                "occurrences_in_protein_table": count,
            }
            for (raw_pathway, canonical_pathway), count in sorted(
                normalization_counts.items(),
                key=lambda x: (x[0][1], x[0][0]),
            )
        ],
        normalization_output,
        [
            "raw_pathway",
            "canonical_pathway",
            "occurrences_in_protein_table",
        ],
    )

    genes_with_pathway = {gene_id for gene_id, _ in gene_pathway_pairs}
    genes_without_pathway = [
        {"gene_id": gene_id}
        for gene_id in sorted(seen_gene_ids - genes_with_pathway)
    ]
    write_tsv(genes_without_pathway, missing_genes_output, ["gene_id"])

    write_tsv(
        [
            {"raw_pathway_token": token, "occurrences": count}
            for token, count in sorted(invalid_token_counts.items())
        ],
        invalid_output,
        ["raw_pathway_token", "occurrences"],
    )

    validation = {
        "script": "21_build_fusarium_kegg_annotation.py",
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "normalization": {
            "accepted_raw_pattern": r"^(ko|map)(\d{5})$",
            "canonical_form": "map#####",
            "example": {"ko00010": "map00010", "map00010": "map00010"},
            "missing_tokens": sorted(MISSING_TOKENS),
        },
        "inputs": {
            "protein_annotation": str(protein_path),
            "gene_annotation": str(gene_path),
        },
        "input_sha256": {
            "protein_annotation": sha256sum(protein_path),
            "gene_annotation": sha256sum(gene_path),
        },
        "counts": observed_counts,
        "expected_counts": None if args.skip_expected_counts else expected_counts,
        "derived_counts": {
            "genes_without_canonical_kegg_pathway": len(genes_without_pathway),
            "canonical_pathways_with_both_ko_and_map_prefixes": sum(
                prefixes == {"ko", "map"}
                for prefixes in pathway_to_prefixes.values()
            ),
            "canonical_pathways_missing_one_prefix": len(incomplete_prefixes),
        },
        "outputs": {
            "protein_kegg_pathways": str(protein_output),
            "gene_kegg_pathways": str(gene_output),
            "term2gene": str(term2gene_output),
            "pathway_metadata": str(metadata_output),
            "pathway_normalization": str(normalization_output),
            "genes_without_kegg_pathway": str(missing_genes_output),
            "invalid_kegg_pathway_tokens": str(invalid_output),
        },
        "validation_passed": True,
    }

    with validation_output.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")

    output_files = [
        protein_output,
        gene_output,
        term2gene_output,
        metadata_output,
        normalization_output,
        missing_genes_output,
        invalid_output,
        validation_output,
    ]

    with checksum_output.open("w", encoding="utf-8") as handle:
        for path in output_files:
            relative_path = path.relative_to(out_dir)
            handle.write(
                f"{sha256sum(path)}  {relative_path.as_posix()}\n"
            )

    print("Fusarium KEGG annotation completed successfully.")
    print(f"Protein rows:                         {protein_rows:,}")
    print(f"Unique proteins:                     {len(protein_ids):,}")
    print(
        "Genes represented in protein table: "
        f"{len(protein_table_genes):,}"
    )
    print(f"Structural Trinity genes:            {len(all_gene_ids):,}")
    print(
        "Proteins with KEGG pathway field:   "
        f"{proteins_with_pathway_field:,}"
    )
    print(f"Unique raw pathway tokens:          {len(raw_token_counts):,}")
    print(f"Canonical KEGG pathways:            {len(pathway_to_prefixes):,}")
    print(f"Unique protein-pathway pairs:       {len(protein_pathway_pairs):,}")
    print(f"Unique gene-pathway pairs:          {len(gene_pathway_pairs):,}")
    print(f"Genes with canonical pathway:       {len(genes_with_pathway):,}")
    print(f"Genes without canonical pathway:    {len(genes_without_pathway):,}")
    print(f"Invalid pathway tokens:             {len(invalid_token_counts):,}")
    print(f"Results written to:                 {out_dir}")


if __name__ == "__main__":
    main()
