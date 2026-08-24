#!/usr/bin/env python3
"""
34_build_cazyme_plotting_tables.py

Create cleaned CAZyme family- and substrate-level plotting tables from the
validated comparative gene-family table produced by script 33.

Input
-----
results/publication/cazyme_comparison/
    comparative_cazyme_gene_family_long.tsv

Outputs
-------
results/publication/cazyme_comparison/plotting_tables/
    cazyme_family_plotting_long.tsv
    cazyme_family_plotting_wide.tsv
    cazyme_substrate_gene_long.tsv
    cazyme_substrate_summary_long.tsv
    cazyme_substrate_summary_wide.tsv
    cazyme_class_summary_long.tsv
    cazyme_class_summary_wide.tsv
    cazyme_plotting_audit.tsv
    run_info.tsv

Counting rules
--------------
1. Family summaries count one gene once per distinct CAZyme family.
2. Substrate annotations are split into standardized biological categories.
3. A gene may contribute to multiple substrate categories when its dbCAN
   annotation contains multiple substrates.
4. Within one family × gene combination, duplicate substrate labels are removed.
5. Class summaries count each gene once per CAZyme class.
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path.cwd()

INPUT = (
    ROOT
    / "results/publication/cazyme_comparison/"
      "comparative_cazyme_gene_family_long.tsv"
)

OUTDIR = (
    ROOT
    / "results/publication/cazyme_comparison/"
      "plotting_tables"
)

FAMILY_LONG_OUT = OUTDIR / "cazyme_family_plotting_long.tsv"
FAMILY_WIDE_OUT = OUTDIR / "cazyme_family_plotting_wide.tsv"

SUBSTRATE_GENE_OUT = OUTDIR / "cazyme_substrate_gene_long.tsv"
SUBSTRATE_LONG_OUT = OUTDIR / "cazyme_substrate_summary_long.tsv"
SUBSTRATE_WIDE_OUT = OUTDIR / "cazyme_substrate_summary_wide.tsv"

CLASS_LONG_OUT = OUTDIR / "cazyme_class_summary_long.tsv"
CLASS_WIDE_OUT = OUTDIR / "cazyme_class_summary_wide.tsv"

AUDIT_OUT = OUTDIR / "cazyme_plotting_audit.tsv"
RUN_INFO_OUT = OUTDIR / "run_info.tsv"

EXPECTED_SPECIES = {
    "Fusarium",
    "Ophiostoma",
}

EXPECTED_GENE_COUNTS = {
    "Fusarium": 142,
    "Ophiostoma": 206,
}

REQUIRED_COLUMNS = {
    "species",
    "gene_id",
    "family",
    "class",
    "substrate",
    "direction",
    "secreted",
    "raw_log2FoldChange",
    "shrunk_log2FoldChange",
    "padj",
}

EMPTY_VALUES = {
    "",
    "-",
    "NA",
    "NaN",
    "nan",
    "None",
    "none",
    "NULL",
    "null",
    "Unassigned",
    "unassigned",
}

CLASS_ORDER = {
    "AA": 0,
    "CBM": 1,
    "CE": 2,
    "GH": 3,
    "GT": 4,
    "PL": 5,
    "Other": 6,
}

SUBSTRATE_ORDER = {
    "Chitin": 0,
    "Chitosan": 1,
    "Beta-glucan": 2,
    "Alpha-glucan": 3,
    "Cellulose": 4,
    "Hemicellulose": 5,
    "Xylan": 6,
    "Xyloglucan": 7,
    "Mannan": 8,
    "Pectin": 9,
    "Arabinan": 10,
    "Arabinogalactan": 11,
    "Fructan": 12,
    "Trehalose": 13,
    "Host glycan": 14,
    "Peptidoglycan": 15,
    "Other": 98,
    "Unassigned": 99,
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
    rows: Iterable[dict[str, object]],
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
    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid {field} for {gene_id}: {value}"
        ) from exc


def parse_bool(value: str, gene_id: str) -> bool:
    normalized = value.strip().upper()

    if normalized == "TRUE":
        return True

    if normalized == "FALSE":
        return False

    fail(
        f"invalid secreted Boolean for {gene_id}: {value}"
    )

    return False


def family_sort_key(family: str) -> tuple[int, int, int, str]:
    cls_match = re.match(
        r"^(AA|CBM|CE|GH|GT|PL)",
        family,
    )

    cls = (
        cls_match.group(1)
        if cls_match
        else "Other"
    )

    number_match = re.match(
        r"^[A-Z]+(\d+)(?:_(\d+))?$",
        family,
    )

    if number_match:
        main_number = int(number_match.group(1))
        subfamily_number = (
            int(number_match.group(2))
            if number_match.group(2)
            else -1
        )
    else:
        main_number = 10**9
        subfamily_number = 10**9

    return (
        CLASS_ORDER.get(cls, 99),
        main_number,
        subfamily_number,
        family,
    )


def substrate_sort_key(
    substrate: str,
) -> tuple[int, str]:
    return (
        SUBSTRATE_ORDER.get(substrate, 98),
        substrate,
    )


def normalize_substrate_token(token: str) -> str:
    text = token.strip().lower()

    text = text.replace("β", "beta")
    text = text.replace("α", "alpha")
    text = re.sub(r"\s+", " ", text)
    text = text.strip(" ._-")

    if text in EMPTY_VALUES:
        return "Unassigned"

    mappings = [
        (r"\bchitosan\b", "Chitosan"),
        (r"\bchitin\b", "Chitin"),
        (r"\bbeta[- ]?glucan\b", "Beta-glucan"),
        (r"\balpha[- ]?glucan\b", "Alpha-glucan"),
        (r"\bcellulose\b", "Cellulose"),
        (r"\bhemicellulose\b", "Hemicellulose"),
        (r"\bxyloglucan\b", "Xyloglucan"),
        (r"\bxylan\b", "Xylan"),
        (r"\barabinogalactan\b", "Arabinogalactan"),
        (r"\barabinan\b", "Arabinan"),
        (r"\bmannan\b", "Mannan"),
        (r"\bpectin\b", "Pectin"),
        (r"\bfructan\b", "Fructan"),
        (r"\btrehalose\b", "Trehalose"),
        (r"\bhost glycan\b", "Host glycan"),
        (r"\bpeptidoglycan\b", "Peptidoglycan"),
    ]

    for pattern, label in mappings:
        if re.search(pattern, text):
            return label

    return "Other"


def parse_substrates(value: str) -> list[str]:
    if value.strip() in EMPTY_VALUES:
        return ["Unassigned"]

    tokens = re.split(
        r"[;|,+/]+",
        value,
    )

    normalized = {
        normalize_substrate_token(token)
        for token in tokens
        if token.strip()
    }

    if not normalized:
        return ["Unassigned"]

    if len(normalized) > 1 and "Unassigned" in normalized:
        normalized.remove("Unassigned")

    return sorted(
        normalized,
        key=substrate_sort_key,
    )


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)

    if n == 0:
        return float("nan")

    midpoint = n // 2

    if n % 2 == 1:
        return ordered[midpoint]

    return (
        ordered[midpoint - 1]
        + ordered[midpoint]
    ) / 2


def summarize_group(
    rows: list[dict[str, object]],
) -> dict[str, object]:
    gene_ids = {
        str(row["gene_id"])
        for row in rows
    }

    up_gene_ids = {
        str(row["gene_id"])
        for row in rows
        if row["direction"] == "up"
    }

    down_gene_ids = {
        str(row["gene_id"])
        for row in rows
        if row["direction"] == "down"
    }

    secreted_gene_ids = {
        str(row["gene_id"])
        for row in rows
        if row["secreted"] is True
    }

    secreted_up_gene_ids = {
        str(row["gene_id"])
        for row in rows
        if (
            row["secreted"] is True
            and row["direction"] == "up"
        )
    }

    secreted_down_gene_ids = {
        str(row["gene_id"])
        for row in rows
        if (
            row["secreted"] is True
            and row["direction"] == "down"
        )
    }

    shrunk_by_gene: dict[str, float] = {}

    for row in rows:
        gene_id = str(row["gene_id"])
        shrunk_by_gene[gene_id] = float(
            row["shrunk_log2FoldChange"]
        )

    shrunk_values = list(shrunk_by_gene.values())

    return {
        "significant_genes": len(gene_ids),
        "up": len(up_gene_ids),
        "down": len(down_gene_ids),
        "secreted": len(secreted_gene_ids),
        "secreted_up": len(secreted_up_gene_ids),
        "secreted_down": len(secreted_down_gene_ids),
        "mean_shrunk_lfc": (
            sum(shrunk_values) / len(shrunk_values)
        ),
        "median_shrunk_lfc": median(shrunk_values),
        "max_abs_shrunk_lfc": max(
            abs(value)
            for value in shrunk_values
        ),
    }


def build_family_summaries(
    rows: list[dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
]:
    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in rows:
        grouped[
            (
                str(row["family"]),
                str(row["species"]),
            )
        ].append(row)

    long_rows: list[dict[str, object]] = []

    for (family, species), group_rows in grouped.items():
        summary = summarize_group(group_rows)

        substrates = sorted(
            {
                substrate
                for row in group_rows
                for substrate in parse_substrates(
                    str(row["substrate"])
                )
            },
            key=substrate_sort_key,
        )

        long_rows.append(
            {
                "family": family,
                "class": str(group_rows[0]["class"]),
                "substrates": " | ".join(substrates),
                "species": species,
                **summary,
            }
        )

    long_rows.sort(
        key=lambda row: (
            family_sort_key(str(row["family"])),
            str(row["species"]),
        )
    )

    index = {
        (
            str(row["family"]),
            str(row["species"]),
        ): row
        for row in long_rows
    }

    families = sorted(
        {
            str(row["family"])
            for row in long_rows
        },
        key=family_sort_key,
    )

    wide_rows: list[dict[str, object]] = []

    for family in families:
        fusarium = index.get(
            (family, "Fusarium")
        )

        ophiostoma = index.get(
            (family, "Ophiostoma")
        )

        source = fusarium or ophiostoma

        row: dict[str, object] = {
            "family": family,
            "class": (
                source["class"]
                if source is not None
                else ""
            ),
            "substrates": " | ".join(
                sorted(
                    {
                        substrate
                        for species_row in (
                            fusarium,
                            ophiostoma,
                        )
                        if species_row is not None
                        for substrate in str(
                            species_row["substrates"]
                        ).split(" | ")
                    },
                    key=substrate_sort_key,
                )
            ),
        }

        for prefix, species_row in (
            ("fusarium", fusarium),
            ("ophiostoma", ophiostoma),
        ):
            for field in (
                "significant_genes",
                "up",
                "down",
                "secreted",
                "secreted_up",
                "secreted_down",
                "mean_shrunk_lfc",
                "median_shrunk_lfc",
                "max_abs_shrunk_lfc",
            ):
                row[f"{prefix}_{field}"] = (
                    species_row[field]
                    if species_row is not None
                    else (
                        ""
                        if "lfc" in field
                        else 0
                    )
                )

        row["total_significant"] = (
            int(row["fusarium_significant_genes"])
            + int(row["ophiostoma_significant_genes"])
        )

        row["total_secreted"] = (
            int(row["fusarium_secreted"])
            + int(row["ophiostoma_secreted"])
        )

        row["total_secreted_up"] = (
            int(row["fusarium_secreted_up"])
            + int(row["ophiostoma_secreted_up"])
        )

        row["total_secreted_down"] = (
            int(row["fusarium_secreted_down"])
            + int(row["ophiostoma_secreted_down"])
        )

        wide_rows.append(row)

    wide_rows.sort(
        key=lambda row: (
            -int(row["total_secreted"]),
            -int(row["total_significant"]),
            family_sort_key(str(row["family"])),
        )
    )

    return long_rows, wide_rows


def build_substrate_gene_rows(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    expanded: dict[
        tuple[str, str, str],
        dict[str, object],
    ] = {}

    for row in rows:
        species = str(row["species"])
        gene_id = str(row["gene_id"])

        for substrate in parse_substrates(
            str(row["substrate"])
        ):
            key = (
                species,
                gene_id,
                substrate,
            )

            if key not in expanded:
                expanded[key] = {
                    "species": species,
                    "gene_id": gene_id,
                    "substrate": substrate,
                    "direction": row["direction"],
                    "secreted": row["secreted"],
                    "raw_log2FoldChange":
                        row["raw_log2FoldChange"],
                    "shrunk_log2FoldChange":
                        row["shrunk_log2FoldChange"],
                    "padj": row["padj"],
                    "families": set(),
                    "classes": set(),
                }

            expanded[key]["families"].add(
                str(row["family"])
            )

            expanded[key]["classes"].add(
                str(row["class"])
            )

    output: list[dict[str, object]] = []

    for row in expanded.values():
        output.append(
            {
                **row,
                "families": "|".join(
                    sorted(
                        row["families"],
                        key=family_sort_key,
                    )
                ),
                "classes": "|".join(
                    sorted(
                        row["classes"],
                        key=lambda cls: (
                            CLASS_ORDER.get(cls, 99),
                            cls,
                        ),
                    )
                ),
            }
        )

    output.sort(
        key=lambda row: (
            substrate_sort_key(
                str(row["substrate"])
            ),
            str(row["species"]),
            str(row["direction"]),
            str(row["gene_id"]),
        )
    )

    return output


def build_substrate_summaries(
    substrate_gene_rows: list[dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
]:
    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in substrate_gene_rows:
        grouped[
            (
                str(row["substrate"]),
                str(row["species"]),
            )
        ].append(row)

    long_rows: list[dict[str, object]] = []

    for (substrate, species), group_rows in grouped.items():
        long_rows.append(
            {
                "substrate": substrate,
                "species": species,
                **summarize_group(group_rows),
            }
        )

    long_rows.sort(
        key=lambda row: (
            substrate_sort_key(
                str(row["substrate"])
            ),
            str(row["species"]),
        )
    )

    index = {
        (
            str(row["substrate"]),
            str(row["species"]),
        ): row
        for row in long_rows
    }

    substrates = sorted(
        {
            str(row["substrate"])
            for row in long_rows
        },
        key=substrate_sort_key,
    )

    wide_rows: list[dict[str, object]] = []

    for substrate in substrates:
        fusarium = index.get(
            (substrate, "Fusarium")
        )

        ophiostoma = index.get(
            (substrate, "Ophiostoma")
        )

        row: dict[str, object] = {
            "substrate": substrate,
        }

        for prefix, species_row in (
            ("fusarium", fusarium),
            ("ophiostoma", ophiostoma),
        ):
            for field in (
                "significant_genes",
                "up",
                "down",
                "secreted",
                "secreted_up",
                "secreted_down",
                "mean_shrunk_lfc",
                "median_shrunk_lfc",
                "max_abs_shrunk_lfc",
            ):
                row[f"{prefix}_{field}"] = (
                    species_row[field]
                    if species_row is not None
                    else (
                        ""
                        if "lfc" in field
                        else 0
                    )
                )

        row["total_significant"] = (
            int(row["fusarium_significant_genes"])
            + int(row["ophiostoma_significant_genes"])
        )

        row["total_secreted"] = (
            int(row["fusarium_secreted"])
            + int(row["ophiostoma_secreted"])
        )

        wide_rows.append(row)

    wide_rows.sort(
        key=lambda row: (
            -int(row["total_secreted"]),
            -int(row["total_significant"]),
            substrate_sort_key(
                str(row["substrate"])
            ),
        )
    )

    return long_rows, wide_rows


def build_class_summaries(
    rows: list[dict[str, object]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
]:
    unique_rows: dict[
        tuple[str, str, str],
        dict[str, object],
    ] = {}

    for row in rows:
        key = (
            str(row["species"]),
            str(row["gene_id"]),
            str(row["class"]),
        )

        unique_rows[key] = row

    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in unique_rows.values():
        grouped[
            (
                str(row["class"]),
                str(row["species"]),
            )
        ].append(row)

    long_rows: list[dict[str, object]] = []

    for (cls, species), group_rows in grouped.items():
        long_rows.append(
            {
                "class": cls,
                "species": species,
                **summarize_group(group_rows),
            }
        )

    long_rows.sort(
        key=lambda row: (
            CLASS_ORDER.get(
                str(row["class"]),
                99,
            ),
            str(row["species"]),
        )
    )

    index = {
        (
            str(row["class"]),
            str(row["species"]),
        ): row
        for row in long_rows
    }

    classes = sorted(
        {
            str(row["class"])
            for row in long_rows
        },
        key=lambda cls: (
            CLASS_ORDER.get(cls, 99),
            cls,
        ),
    )

    wide_rows: list[dict[str, object]] = []

    for cls in classes:
        fusarium = index.get(
            (cls, "Fusarium")
        )

        ophiostoma = index.get(
            (cls, "Ophiostoma")
        )

        row: dict[str, object] = {
            "class": cls,
        }

        for prefix, species_row in (
            ("fusarium", fusarium),
            ("ophiostoma", ophiostoma),
        ):
            for field in (
                "significant_genes",
                "up",
                "down",
                "secreted",
                "secreted_up",
                "secreted_down",
                "mean_shrunk_lfc",
                "median_shrunk_lfc",
                "max_abs_shrunk_lfc",
            ):
                row[f"{prefix}_{field}"] = (
                    species_row[field]
                    if species_row is not None
                    else (
                        ""
                        if "lfc" in field
                        else 0
                    )
                )

        wide_rows.append(row)

    return long_rows, wide_rows


def main() -> None:
    fields, input_rows = read_tsv(INPUT)

    missing = REQUIRED_COLUMNS - set(fields)

    if missing:
        fail(
            f"input is missing columns: {sorted(missing)}"
        )

    parsed_rows: list[dict[str, object]] = []

    observed_species: set[str] = set()

    for row in input_rows:
        gene_id = row["gene_id"]
        species = row["species"]
        direction = row["direction"].lower()

        observed_species.add(species)

        if direction not in {"up", "down"}:
            fail(
                f"unexpected direction for {gene_id}: "
                f"{row['direction']}"
            )

        parsed_rows.append(
            {
                "species": species,
                "gene_id": gene_id,
                "family": row["family"],
                "class": row["class"],
                "substrate": row["substrate"],
                "direction": direction,
                "secreted": parse_bool(
                    row["secreted"],
                    gene_id,
                ),
                "raw_log2FoldChange": parse_float(
                    row["raw_log2FoldChange"],
                    "raw_log2FoldChange",
                    gene_id,
                ),
                "shrunk_log2FoldChange": parse_float(
                    row["shrunk_log2FoldChange"],
                    "shrunk_log2FoldChange",
                    gene_id,
                ),
                "padj": parse_float(
                    row["padj"],
                    "padj",
                    gene_id,
                ),
            }
        )

    if observed_species != EXPECTED_SPECIES:
        fail(
            f"expected species {sorted(EXPECTED_SPECIES)}, "
            f"observed {sorted(observed_species)}"
        )

    observed_gene_counts = {
        species: len(
            {
                str(row["gene_id"])
                for row in parsed_rows
                if row["species"] == species
            }
        )
        for species in EXPECTED_SPECIES
    }

    for species, expected in EXPECTED_GENE_COUNTS.items():
        observed = observed_gene_counts[species]

        if observed != expected:
            fail(
                f"{species}: expected {expected} unique genes, "
                f"observed {observed}"
            )

    family_long, family_wide = build_family_summaries(
        parsed_rows
    )

    substrate_gene_rows = build_substrate_gene_rows(
        parsed_rows
    )

    substrate_long, substrate_wide = (
        build_substrate_summaries(
            substrate_gene_rows
        )
    )

    class_long, class_wide = build_class_summaries(
        parsed_rows
    )

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_tsv(
        FAMILY_LONG_OUT,
        [
            "family",
            "class",
            "substrates",
            "species",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
            "mean_shrunk_lfc",
            "median_shrunk_lfc",
            "max_abs_shrunk_lfc",
        ],
        family_long,
    )

    write_tsv(
        FAMILY_WIDE_OUT,
        list(family_wide[0]),
        family_wide,
    )

    write_tsv(
        SUBSTRATE_GENE_OUT,
        [
            "species",
            "gene_id",
            "substrate",
            "direction",
            "secreted",
            "raw_log2FoldChange",
            "shrunk_log2FoldChange",
            "padj",
            "families",
            "classes",
        ],
        substrate_gene_rows,
    )

    write_tsv(
        SUBSTRATE_LONG_OUT,
        [
            "substrate",
            "species",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
            "mean_shrunk_lfc",
            "median_shrunk_lfc",
            "max_abs_shrunk_lfc",
        ],
        substrate_long,
    )

    write_tsv(
        SUBSTRATE_WIDE_OUT,
        list(substrate_wide[0]),
        substrate_wide,
    )

    write_tsv(
        CLASS_LONG_OUT,
        [
            "class",
            "species",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
            "mean_shrunk_lfc",
            "median_shrunk_lfc",
            "max_abs_shrunk_lfc",
        ],
        class_long,
    )

    write_tsv(
        CLASS_WIDE_OUT,
        list(class_wide[0]),
        class_wide,
    )

    audit_rows: list[dict[str, object]] = []

    for species in sorted(EXPECTED_SPECIES):
        species_family_rows = [
            row
            for row in parsed_rows
            if row["species"] == species
        ]

        species_substrate_rows = [
            row
            for row in substrate_gene_rows
            if row["species"] == species
        ]

        unique_genes = {
            str(row["gene_id"])
            for row in species_family_rows
        }

        secreted_genes = {
            str(row["gene_id"])
            for row in species_family_rows
            if row["secreted"] is True
        }

        unique_families = {
            str(row["family"])
            for row in species_family_rows
        }

        unique_substrates = {
            str(row["substrate"])
            for row in species_substrate_rows
        }

        audit_rows.append(
            {
                "species": species,
                "unique_significant_cazyme_genes":
                    len(unique_genes),
                "family_assignment_rows":
                    len(species_family_rows),
                "unique_families":
                    len(unique_families),
                "substrate_assignment_rows":
                    len(species_substrate_rows),
                "unique_substrates":
                    len(unique_substrates),
                "secreted_genes":
                    len(secreted_genes),
                "status": "PASS",
            }
        )

    write_tsv(
        AUDIT_OUT,
        [
            "species",
            "unique_significant_cazyme_genes",
            "family_assignment_rows",
            "unique_families",
            "substrate_assignment_rows",
            "unique_substrates",
            "secreted_genes",
            "status",
        ],
        audit_rows,
    )

    run_info = [
        {
            "field": "script",
            "value":
                "34_build_cazyme_plotting_tables.py",
        },
        {
            "field": "run_timestamp_utc",
            "value":
                datetime.now(timezone.utc).isoformat(),
        },
        {
            "field": "input_file",
            "value": str(INPUT),
        },
        {
            "field": "input_sha256",
            "value": sha256(INPUT),
        },
        {
            "field": "family_counting_rule",
            "value":
                "one gene once per distinct CAZyme family",
        },
        {
            "field": "class_counting_rule",
            "value":
                "one gene once per CAZyme class",
        },
        {
            "field": "substrate_counting_rule",
            "value": (
                "one gene once per standardized substrate; "
                "multi-substrate annotations are expanded"
            ),
        },
        {
            "field": "direction",
            "value":
                "sign of raw log2 fold change",
        },
        {
            "field": "effect_size",
            "value":
                "apeglm-shrunken log2 fold change",
        },
    ]

    write_tsv(
        RUN_INFO_OUT,
        ["field", "value"],
        run_info,
    )

    print()
    print("============================================================")
    print("CAZYME PLOTTING TABLES COMPLETE")
    print("============================================================")
    print()

    for audit in audit_rows:
        print(audit["species"])
        print(
            "  Significant CAZyme genes: "
            f"{audit['unique_significant_cazyme_genes']}"
        )
        print(
            "  Unique families: "
            f"{audit['unique_families']}"
        )
        print(
            "  Standardized substrates: "
            f"{audit['unique_substrates']}"
        )
        print(
            "  Secreted genes: "
            f"{audit['secreted_genes']}"
        )
        print()

    print(f"Family table:    {FAMILY_WIDE_OUT}")
    print(f"Substrate table: {SUBSTRATE_WIDE_OUT}")
    print(f"Class table:     {CLASS_WIDE_OUT}")
    print(f"Audit:           {AUDIT_OUT}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
