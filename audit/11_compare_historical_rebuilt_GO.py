#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

HISTORICAL_DIR = ROOT / "audit/go_historical_comparison/input"
REBUILT_FILE = (
    ROOT
    / "results/ophiostoma/go_enrichment/all_significant_GO_terms.tsv"
)
OUTPUT_DIR = ROOT / "audit/go_historical_comparison/output"


CONTRASTS = {
    "interaction_vs_self": {
        "historical_file": "GO_interaction_vs_self.csv",
        "historical_label": "interaction_vs_self",
    },
    "interaction_vs_onu": {
        "historical_file": "GO_interaction_vs_ONU.csv",
        "historical_label": "interaction_vs_onu",
    },
    "onu_vs_self": {
        "historical_file": "GO_ONU_vs_self.csv",
        "historical_label": "onu_vs_self",
    },
}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def load_historical(path: Path, contrast: str) -> pd.DataFrame:
    if not path.is_file():
        fail(f"Historical file not found: {path}")

    df = pd.read_csv(path)

    required = {
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "FoldEnrichment",
        "pvalue",
        "p.adjust",
        "Count",
    }

    missing = required.difference(df.columns)
    if missing:
        fail(
            f"{path.name} is missing columns: "
            + ", ".join(sorted(missing))
        )

    result = df.copy()
    result["contrast"] = contrast
    result["historical_present"] = True

    return result


def load_rebuilt(path: Path) -> pd.DataFrame:
    if not path.is_file():
        fail(f"Rebuilt GO file not found: {path}")

    df = pd.read_csv(path, sep="\t")

    required = {
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "FoldEnrichment",
        "pvalue",
        "p.adjust",
        "Count",
        "contrast",
        "gene_set",
        "direction",
        "ontology",
    }

    missing = required.difference(df.columns)
    if missing:
        fail(
            f"{path.name} is missing columns: "
            + ", ".join(sorted(missing))
        )

    # Main comparison uses the padj < 0.05 gene sets,
    # not the stronger |log2FC|-filtered sensitivity sets.
    result = df.loc[df["gene_set"].eq("significant")].copy()

    duplicated = result.duplicated(
        subset=["contrast", "direction", "ontology", "ID"]
    )

    if duplicated.any():
        fail(
            "Duplicate rebuilt GO entries detected for the same "
            "contrast, direction, ontology and GO ID."
        )

    return result


def compare_contrast(
    historical: pd.DataFrame,
    rebuilt: pd.DataFrame,
    contrast: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rebuilt_contrast = rebuilt.loc[
        rebuilt["contrast"].eq(contrast)
    ].copy()

    historical_ids = set(historical["ID"].dropna())
    rebuilt_ids = set(rebuilt_contrast["ID"].dropna())

    shared_ids = historical_ids & rebuilt_ids
    historical_only_ids = historical_ids - rebuilt_ids
    rebuilt_only_ids = rebuilt_ids - historical_ids

    shared = rebuilt_contrast.loc[
        rebuilt_contrast["ID"].isin(shared_ids)
    ].copy()

    historical_lookup = (
        historical[
            [
                "ID",
                "Description",
                "GeneRatio",
                "BgRatio",
                "FoldEnrichment",
                "pvalue",
                "p.adjust",
                "Count",
            ]
        ]
        .rename(
            columns={
                "Description": "historical_description",
                "GeneRatio": "historical_gene_ratio",
                "BgRatio": "historical_bg_ratio",
                "FoldEnrichment": "historical_fold_enrichment",
                "pvalue": "historical_pvalue",
                "p.adjust": "historical_p_adjust",
                "Count": "historical_count",
            }
        )
    )

    shared = shared.merge(
        historical_lookup,
        on="ID",
        how="left",
        validate="many_to_one",
    )

    shared = shared.rename(
        columns={
            "Description": "rebuilt_description",
            "GeneRatio": "rebuilt_gene_ratio",
            "BgRatio": "rebuilt_bg_ratio",
            "FoldEnrichment": "rebuilt_fold_enrichment",
            "pvalue": "rebuilt_pvalue",
            "p.adjust": "rebuilt_p_adjust",
            "Count": "rebuilt_count",
        }
    )

    shared["comparison_status"] = "exact_GO_ID_reproduced"

    shared_columns = [
        "contrast",
        "ID",
        "historical_description",
        "rebuilt_description",
        "direction",
        "ontology",
        "historical_gene_ratio",
        "rebuilt_gene_ratio",
        "historical_bg_ratio",
        "rebuilt_bg_ratio",
        "historical_fold_enrichment",
        "rebuilt_fold_enrichment",
        "historical_p_adjust",
        "rebuilt_p_adjust",
        "historical_count",
        "rebuilt_count",
        "comparison_status",
    ]

    shared = shared[shared_columns].sort_values(
        ["direction", "ontology", "rebuilt_p_adjust", "ID"]
    )

    historical_only = historical.loc[
        historical["ID"].isin(historical_only_ids)
    ].copy()

    historical_only["comparison_status"] = (
        "historical_only_exact_GO_ID"
    )

    historical_only = historical_only[
        [
            "contrast",
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "FoldEnrichment",
            "pvalue",
            "p.adjust",
            "Count",
            "comparison_status",
        ]
    ].sort_values(["p.adjust", "ID"])

    rebuilt_only = rebuilt_contrast.loc[
        rebuilt_contrast["ID"].isin(rebuilt_only_ids)
    ].copy()

    rebuilt_only["comparison_status"] = "rebuilt_only_exact_GO_ID"

    rebuilt_only = rebuilt_only[
        [
            "contrast",
            "direction",
            "ontology",
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "FoldEnrichment",
            "pvalue",
            "p.adjust",
            "Count",
            "comparison_status",
        ]
    ].sort_values(
        ["direction", "ontology", "p.adjust", "ID"]
    )

    return shared, historical_only, rebuilt_only


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rebuilt = load_rebuilt(REBUILT_FILE)

    all_shared = []
    all_historical_only = []
    all_rebuilt_only = []
    summary_rows = []

    for contrast, config in CONTRASTS.items():
        historical_path = (
            HISTORICAL_DIR / config["historical_file"]
        )

        historical = load_historical(
            historical_path,
            contrast=config["historical_label"],
        )

        rebuilt_contrast = rebuilt.loc[
            rebuilt["contrast"].eq(contrast)
        ]

        shared, historical_only, rebuilt_only = compare_contrast(
            historical=historical,
            rebuilt=rebuilt,
            contrast=contrast,
        )

        contrast_dir = OUTPUT_DIR / contrast
        contrast_dir.mkdir(parents=True, exist_ok=True)

        shared.to_csv(
            contrast_dir / "shared_exact_GO_terms.tsv",
            sep="\t",
            index=False,
        )

        historical_only.to_csv(
            contrast_dir / "historical_only_GO_terms.tsv",
            sep="\t",
            index=False,
        )

        rebuilt_only.to_csv(
            contrast_dir / "rebuilt_only_GO_terms.tsv",
            sep="\t",
            index=False,
        )

        all_shared.append(shared)
        all_historical_only.append(historical_only)
        all_rebuilt_only.append(rebuilt_only)

        summary_rows.append(
            {
                "contrast": contrast,
                "historical_significant_terms": len(historical),
                "rebuilt_significant_rows": len(rebuilt_contrast),
                "rebuilt_unique_GO_IDs": rebuilt_contrast["ID"].nunique(),
                "shared_unique_GO_IDs": shared["ID"].nunique(),
                "historical_only_unique_GO_IDs": (
                    historical_only["ID"].nunique()
                ),
                "rebuilt_only_unique_GO_IDs": (
                    rebuilt_only["ID"].nunique()
                ),
                "shared_up_rows": int(
                    shared["direction"].eq("up").sum()
                ),
                "shared_down_rows": int(
                    shared["direction"].eq("down").sum()
                ),
            }
        )

    pd.concat(all_shared, ignore_index=True).to_csv(
        OUTPUT_DIR / "all_shared_exact_GO_terms.tsv",
        sep="\t",
        index=False,
    )

    pd.concat(all_historical_only, ignore_index=True).to_csv(
        OUTPUT_DIR / "all_historical_only_GO_terms.tsv",
        sep="\t",
        index=False,
    )

    pd.concat(all_rebuilt_only, ignore_index=True).to_csv(
        OUTPUT_DIR / "all_rebuilt_only_GO_terms.tsv",
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        OUTPUT_DIR / "GO_exact_overlap_summary.tsv",
        sep="\t",
        index=False,
    )

    ontology_summary = (
        rebuilt.groupby(
            ["contrast", "direction", "ontology"],
            observed=True,
        )
        .size()
        .reset_index(name="significant_term_rows")
    )

    ontology_summary.to_csv(
        OUTPUT_DIR / "rebuilt_significant_term_counts.tsv",
        sep="\t",
        index=False,
    )

    print("PASS: historical versus rebuilt GO comparison completed")
    print(f"Rebuilt significant rows: {len(rebuilt)}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
