#!/usr/bin/env python3

"""
Build a validated gene-level KEGG pathway annotation for Ophiostoma novo-ulmi.

Inputs
------
1. Raw JGI/MycoCosm KEGG annotation table.
2. Validated gene_id-to-protein_id annotation map produced by script 09.

Outputs
-------
- Full row-level gene–KEGG annotation
- TERM2GENE table for clusterProfiler::enricher()
- Pathway metadata table
- Genes without a usable KEGG pathway
- Validation JSON
- SHA256 checksum file

A usable KEGG pathway is a row where the pathway field is not empty or '\\N'.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


DEFAULT_KEGG = Path(
    "data/annotation/ophiostoma/KEGG/"
    "Ophnu1_GeneCatalog_proteins_20170425_KEGG.tab"
)

DEFAULT_GENE_MAP = Path(
    "results/ophiostoma/annotation/"
    "ophiostoma_gene_annotation_map.tsv"
)

DEFAULT_OUT_DIR = Path(
    "results/ophiostoma/kegg_annotation"
)

MISSING_TOKENS = {"", r"\N", "NA", "NaN", "nan", "None"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build validated gene-level KEGG pathway annotation."
    )
    parser.add_argument("--kegg", type=Path, default=DEFAULT_KEGG)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--expected-genes", type=int, default=8640)
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="")


def main() -> None:
    args = parse_args()

    args.kegg = args.kegg.expanduser().resolve()
    args.gene_map = args.gene_map.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()

    require_file(args.kegg, "KEGG annotation table")
    require_file(args.gene_map, "Gene annotation map")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    kegg = pd.read_csv(
        args.kegg,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    gene_map = pd.read_csv(
        args.gene_map,
        sep="\t",
        dtype=str,
        usecols=["gene_id", "protein_id"],
        keep_default_na=False,
    )

    if "#proteinId" in kegg.columns:
        kegg = kegg.rename(columns={"#proteinId": "protein_id"})

    required_kegg_columns = {
        "protein_id",
        "ecNum",
        "definition",
        "catalyticActivity",
        "cofactors",
        "associatedDiseases",
        "pathway",
        "pathway_class",
        "pathway_type",
    }

    missing_columns = required_kegg_columns - set(kegg.columns)
    if missing_columns:
        raise ValueError(
            "KEGG table is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    if len(gene_map) != args.expected_genes:
        raise ValueError(
            f"Expected {args.expected_genes} gene-map rows, "
            f"observed {len(gene_map)}."
        )

    if gene_map["gene_id"].duplicated().any():
        duplicated = gene_map.loc[
            gene_map["gene_id"].duplicated(keep=False), "gene_id"
        ].unique()
        raise ValueError(
            "Duplicated gene IDs in gene map: "
            + ", ".join(map(str, duplicated[:10]))
        )

    if gene_map["protein_id"].duplicated().any():
        duplicated = gene_map.loc[
            gene_map["protein_id"].duplicated(keep=False), "protein_id"
        ].unique()
        raise ValueError(
            "Protein IDs mapping to multiple genes: "
            + ", ".join(map(str, duplicated[:10]))
        )

    for column in kegg.columns:
        kegg[column] = kegg[column].astype(str).str.strip()

    valid = kegg[
        ~kegg["pathway"].isin(MISSING_TOKENS)
    ].copy()

    raw_duplicate_pairs = int(
        valid.duplicated(subset=["protein_id", "pathway"]).sum()
    )

    valid = valid.drop_duplicates(
        subset=["protein_id", "pathway"]
    ).copy()

    annotation = valid.merge(
        gene_map,
        on="protein_id",
        how="left",
        validate="many_to_one",
    )

    missing_gene_ids = annotation["gene_id"].isna() | (
        annotation["gene_id"].str.strip() == ""
    )

    if missing_gene_ids.any():
        missing_proteins = sorted(
            annotation.loc[missing_gene_ids, "protein_id"].unique()
        )
        raise ValueError(
            "KEGG pathway proteins missing from gene map: "
            + ", ".join(missing_proteins[:20])
        )

    annotation = annotation[
        [
            "gene_id",
            "protein_id",
            "ecNum",
            "definition",
            "catalyticActivity",
            "cofactors",
            "associatedDiseases",
            "pathway",
            "pathway_class",
            "pathway_type",
        ]
    ].rename(
        columns={
            "ecNum": "ec_num",
            "catalyticActivity": "catalytic_activity",
            "associatedDiseases": "associated_diseases",
        }
    )

    annotation = annotation.sort_values(
        ["pathway", "gene_id", "protein_id"]
    ).reset_index(drop=True)

    term2gene = (
        annotation[["pathway", "gene_id"]]
        .drop_duplicates()
        .sort_values(["pathway", "gene_id"])
        .reset_index(drop=True)
    )

    pathway_metadata = (
        annotation[
            ["pathway", "pathway_class", "pathway_type"]
        ]
        .drop_duplicates()
        .sort_values(["pathway", "pathway_class", "pathway_type"])
        .reset_index(drop=True)
    )

    class_counts = (
        pathway_metadata.groupby("pathway")["pathway_class"]
        .nunique()
    )
    type_counts = (
        pathway_metadata.groupby("pathway")["pathway_type"]
        .nunique()
    )

    class_collisions = class_counts[class_counts > 1]
    type_collisions = type_counts[type_counts > 1]

    if not class_collisions.empty:
        raise ValueError(
            "Pathway names associated with multiple pathway classes: "
            + ", ".join(class_collisions.index[:20])
        )

    if not type_collisions.empty:
        raise ValueError(
            "Pathway names associated with multiple pathway types: "
            + ", ".join(type_collisions.index[:20])
        )

    annotated_gene_ids = set(term2gene["gene_id"])
    all_gene_ids = set(gene_map["gene_id"])

    genes_without_pathway = (
        gene_map.loc[
            ~gene_map["gene_id"].isin(annotated_gene_ids),
            ["gene_id", "protein_id"],
        ]
        .sort_values("gene_id")
        .reset_index(drop=True)
    )

    annotation_path = (
        args.out_dir / "ophiostoma_gene_kegg_pathways.tsv"
    )
    term2gene_path = (
        args.out_dir / "ophiostoma_kegg_term2gene.tsv"
    )
    metadata_path = (
        args.out_dir / "ophiostoma_kegg_pathway_metadata.tsv"
    )
    missing_path = (
        args.out_dir / "genes_without_kegg_pathway.tsv"
    )
    validation_path = (
        args.out_dir
        / "12_build_ophiostoma_kegg_annotation.validation.json"
    )
    checksum_path = (
        args.out_dir
        / "12_build_ophiostoma_kegg_annotation.sha256"
    )

    write_tsv(annotation, annotation_path)
    write_tsv(term2gene, term2gene_path)
    write_tsv(pathway_metadata, metadata_path)
    write_tsv(genes_without_pathway, missing_path)

    validation = {
        "run_timestamp": datetime.now().astimezone().isoformat(
            timespec="seconds"
        ),
        "inputs": {
            "kegg": str(args.kegg),
            "gene_map": str(args.gene_map),
        },
        "input_sha256": {
            "kegg": sha256sum(args.kegg),
            "gene_map": sha256sum(args.gene_map),
        },
        "counts": {
            "total_structural_genes": len(gene_map),
            "raw_kegg_rows": len(kegg),
            "unique_proteins_anywhere_in_kegg": int(
                kegg["protein_id"].nunique()
            ),
            "raw_rows_with_valid_pathway": len(
                kegg[~kegg["pathway"].isin(MISSING_TOKENS)]
            ),
            "unique_proteins_with_valid_pathway": int(
                valid["protein_id"].nunique()
            ),
            "unique_genes_with_valid_pathway": int(
                term2gene["gene_id"].nunique()
            ),
            "unique_gene_pathway_pairs": len(term2gene),
            "unique_pathways": int(term2gene["pathway"].nunique()),
            "genes_without_valid_pathway": len(genes_without_pathway),
            "duplicate_raw_protein_pathway_pairs": raw_duplicate_pairs,
            "pathway_class_collisions": int(len(class_collisions)),
            "pathway_type_collisions": int(len(type_collisions)),
            "unmapped_pathway_proteins": int(missing_gene_ids.sum()),
        },
        "coverage": {
            "fraction_of_structural_genes_with_valid_pathway": (
                term2gene["gene_id"].nunique() / len(gene_map)
            )
        },
        "outputs": {
            "gene_kegg_pathways": str(annotation_path),
            "term2gene": str(term2gene_path),
            "pathway_metadata": str(metadata_path),
            "genes_without_kegg_pathway": str(missing_path),
        },
        "validation_passed": True,
    }

    with validation_path.open("w", encoding="utf-8") as handle:
        json.dump(validation, handle, indent=2, sort_keys=True)
        handle.write("\n")

    output_files = [
        annotation_path,
        term2gene_path,
        metadata_path,
        missing_path,
        validation_path,
    ]

    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in output_files:
            handle.write(f"{sha256sum(path)}  {path.name}\n")

    print("KEGG pathway annotation completed successfully.")
    print(f"Structural genes:                 {len(gene_map):,}")
    print(
        "Genes with usable KEGG pathway:  "
        f"{term2gene['gene_id'].nunique():,}"
    )
    print(
        "Genes without usable pathway:    "
        f"{len(genes_without_pathway):,}"
    )
    print(f"Unique gene-pathway pairs:        {len(term2gene):,}")
    print(
        "Unique KEGG pathways:             "
        f"{term2gene['pathway'].nunique():,}"
    )
    print(f"Results written to:               {args.out_dir}")


if __name__ == "__main__":
    main()

