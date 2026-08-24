#!/usr/bin/env python3
"""
33_build_comparative_cazyme_family_table.py

Build gene-level and family-level comparative CAZyme tables for the
biologically matched interaction-versus-self contrasts in:

    Fusarium cf. salinense
    Ophiostoma novo-ulmi

Primary input subset
--------------------
High-confidence CAZymes:
    dbCAN support from at least two tools

Differential-expression definition:
    padj < 0.05

Direction:
    sign of the raw log2 fold change

Effect size:
    shrunken log2 fold change retained for ranking and interpretation

Secreted status:
    confident SignalP-positive gene annotation

Multi-family proteins
---------------------
A gene assigned to multiple CAZyme families contributes once to each distinct
family. Duplicate representations of the same family within a gene are removed.

Outputs
-------
results/publication/cazyme_comparison/
    comparative_cazyme_gene_family_long.tsv
    comparative_cazyme_family_summary_long.tsv
    comparative_cazyme_family_summary_wide.tsv
    comparative_cazyme_class_summary.tsv
    comparative_cazyme_substrate_summary.tsv
    comparative_cazyme_audit.tsv
    run_info.tsv
"""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path.cwd()

FUSARIUM_INPUT = (
    ROOT
    / "results/fusarium/cazyme_analysis/tables/"
      "fusarium_cazymes_high_confidence_significant.tsv"
)

OPHIOSTOMA_INPUT = (
    ROOT
    / "results/ophiostoma/cazyme_analysis/tables/"
      "interaction_vs_onu_high_confidence_cazyme_significant.tsv"
)

OUTDIR = ROOT / "results/publication/cazyme_comparison"

GENE_FAMILY_OUT = OUTDIR / "comparative_cazyme_gene_family_long.tsv"
FAMILY_LONG_OUT = OUTDIR / "comparative_cazyme_family_summary_long.tsv"
FAMILY_WIDE_OUT = OUTDIR / "comparative_cazyme_family_summary_wide.tsv"
CLASS_OUT = OUTDIR / "comparative_cazyme_class_summary.tsv"
SUBSTRATE_OUT = OUTDIR / "comparative_cazyme_substrate_summary.tsv"
AUDIT_OUT = OUTDIR / "comparative_cazyme_audit.tsv"
RUN_INFO_OUT = OUTDIR / "run_info.tsv"

EXPECTED_SPECIES_COUNTS = {
    "Fusarium": 142,
    "Ophiostoma": 206,
}

FAMILY_PATTERN = re.compile(
    r"^(AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?$",
    flags=re.IGNORECASE,
)

FAMILY_TOKEN_PATTERN = re.compile(
    r"(AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?",
    flags=re.IGNORECASE,
)

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


def choose_column(
    fields: list[str],
    candidates: tuple[str, ...],
    label: str,
    path: Path,
    required: bool = True,
) -> str | None:
    for candidate in candidates:
        if candidate in fields:
            return candidate

    if required:
        fail(
            f"could not identify {label} column in {path}\n"
            f"Expected one of: {', '.join(candidates)}\n"
            f"Observed columns: {', '.join(fields)}"
        )

    return None


def parse_float(
    value: str,
    field: str,
    gene_id: str,
) -> float:
    if value in EMPTY_VALUES:
        return float("nan")

    try:
        return float(value)
    except ValueError as exc:
        raise SystemExit(
            f"ERROR: invalid {field} for {gene_id}: {value}"
        ) from exc


def is_nan(value: float) -> bool:
    return value != value


def parse_bool(value: str) -> bool:
    normalized = value.strip().upper()

    if normalized in {"TRUE", "T", "1", "YES", "Y"}:
        return True

    if normalized in {"FALSE", "F", "0", "NO", "N", ""}:
        return False

    fail(f"could not parse Boolean value: {value}")
    return False


def normalize_direction(
    supplied_direction: str,
    raw_lfc: float,
    gene_id: str,
) -> str:
    if is_nan(raw_lfc):
        fail(f"missing raw log2 fold change for {gene_id}")

    expected = (
        "up"
        if raw_lfc > 0
        else "down"
        if raw_lfc < 0
        else "unchanged"
    )

    supplied = supplied_direction.strip().lower()

    aliases = {
        "upregulated": "up",
        "induced": "up",
        "positive": "up",
        "downregulated": "down",
        "repressed": "down",
        "negative": "down",
        "stable": "unchanged",
    }

    supplied = aliases.get(supplied, supplied)

    if supplied and supplied not in {"na", "nan"}:
        if supplied != expected:
            fail(
                f"direction mismatch for {gene_id}: "
                f"column says {supplied}, raw LFC implies {expected}"
            )

    return expected


def canonical_family(token: str) -> str | None:
    token = token.strip().upper()

    if token in EMPTY_VALUES:
        return None

    token = re.sub(r"_E\d+$", "", token)
    token = re.sub(r"\([^)]*\)", "", token)
    token = token.strip()

    match = FAMILY_TOKEN_PATTERN.search(token)

    if not match:
        return None

    family = match.group(0).upper()

    if not FAMILY_PATTERN.match(family):
        return None

    return family


def parse_families(value: str) -> list[str]:
    if value.strip() in EMPTY_VALUES:
        return []

    raw_tokens = re.split(
        r"[|+,;/\s]+",
        value.strip(),
    )

    families: set[str] = set()

    for token in raw_tokens:
        family = canonical_family(token)

        if family is not None:
            families.add(family)

    return sorted(
        families,
        key=family_sort_key,
    )


def cazyme_class(family: str) -> str:
    match = re.match(r"^(AA|CBM|CE|GH|GT|PL)", family)

    if match is None:
        return "Other"

    return match.group(1)


def family_sort_key(family: str) -> tuple[int, int, int, str]:
    class_order = {
        "AA": 0,
        "CBM": 1,
        "CE": 2,
        "GH": 3,
        "GT": 4,
        "PL": 5,
        "Other": 6,
    }

    cls = cazyme_class(family)

    match = re.match(
        r"^[A-Z]+(\d+)(?:_(\d+))?$",
        family,
    )

    if match:
        major = int(match.group(1))
        subfamily = (
            int(match.group(2))
            if match.group(2) is not None
            else -1
        )
    else:
        major = 10**9
        subfamily = 10**9

    return (
        class_order.get(cls, 99),
        major,
        subfamily,
        family,
    )


def normalize_substrate(value: str) -> str:
    value = value.strip()

    if value in EMPTY_VALUES:
        return "Unassigned"

    return value


def detect_schema(
    species: str,
    path: Path,
    fields: list[str],
) -> dict[str, str | None]:
    common = {
        "gene_id": choose_column(
            fields,
            ("gene_id",),
            "gene ID",
            path,
        ),
        "family": choose_column(
            fields,
            (
                "dbcan_recommended",
                "dbcan_family",
                "cazyme_family",
            ),
            "recommended CAZyme family",
            path,
        ),
        "substrate": choose_column(
            fields,
            (
                "dbcan_substrate",
                "substrate",
            ),
            "CAZyme substrate",
            path,
            required=False,
        ),
        "direction": choose_column(
            fields,
            (
                "de_direction",
                "direction",
                "regulation",
            ),
            "DE direction",
            path,
        ),
    }

    if species == "Fusarium":
        common.update(
            {
                "raw_lfc": choose_column(
                    fields,
                    (
                        "raw_log2FoldChange",
                        "log2FoldChange",
                    ),
                    "raw log2 fold change",
                    path,
                ),
                "shrunk_lfc": choose_column(
                    fields,
                    (
                        "shrunk_log2FoldChange",
                        "plot_lfc",
                    ),
                    "shrunken log2 fold change",
                    path,
                ),
                "padj": choose_column(
                    fields,
                    (
                        "raw_padj",
                        "padj",
                    ),
                    "adjusted P value",
                    path,
                ),
                "secreted": choose_column(
                    fields,
                    (
                        "signalp_confident_positive",
                        "signalp_is_sp",
                        "signalp_raw_positive",
                    ),
                    "confident SignalP status",
                    path,
                ),
            }
        )

    elif species == "Ophiostoma":
        common.update(
            {
                "raw_lfc": choose_column(
                    fields,
                    ("log2FoldChange",),
                    "raw log2 fold change",
                    path,
                ),
                "shrunk_lfc": choose_column(
                    fields,
                    ("shrunk_log2FoldChange",),
                    "shrunken log2 fold change",
                    path,
                ),
                "padj": choose_column(
                    fields,
                    ("padj",),
                    "adjusted P value",
                    path,
                ),
                "secreted": choose_column(
                    fields,
                    (
                        "signalp_is_sp",
                        "signalp_confident_positive",
                    ),
                    "confident SignalP status",
                    path,
                ),
            }
        )

    else:
        fail(f"unsupported species: {species}")

    return common


def process_species(
    species: str,
    path: Path,
) -> tuple[
    list[dict[str, object]],
    dict[str, object],
]:
    fields, rows = read_tsv(path)
    schema = detect_schema(species, path, fields)

    expected_count = EXPECTED_SPECIES_COUNTS[species]

    if len(rows) != expected_count:
        fail(
            f"{species}: expected {expected_count} significant "
            f"high-confidence CAZyme genes, observed {len(rows)}"
        )

    gene_ids: set[str] = set()
    long_rows: list[dict[str, object]] = []

    genes_without_family: list[str] = []
    genes_with_multiple_families = 0
    secreted_gene_ids: set[str] = set()

    direction_gene_counts = Counter()
    family_assignment_counts = Counter()

    for row in rows:
        gene_id_column = str(schema["gene_id"])
        gene_id = row[gene_id_column]

        if gene_id in gene_ids:
            fail(
                f"{species}: duplicated gene ID in input: {gene_id}"
            )

        gene_ids.add(gene_id)

        raw_lfc = parse_float(
            row[str(schema["raw_lfc"])],
            str(schema["raw_lfc"]),
            gene_id,
        )

        shrunk_lfc = parse_float(
            row[str(schema["shrunk_lfc"])],
            str(schema["shrunk_lfc"]),
            gene_id,
        )

        padj = parse_float(
            row[str(schema["padj"])],
            str(schema["padj"]),
            gene_id,
        )

        if is_nan(padj) or padj >= 0.05:
            fail(
                f"{species}: non-significant gene in significant "
                f"CAZyme input: {gene_id}, padj={padj}"
            )

        direction = normalize_direction(
            row[str(schema["direction"])],
            raw_lfc,
            gene_id,
        )

        if direction not in {"up", "down"}:
            fail(
                f"{species}: significant CAZyme has direction "
                f"{direction}: {gene_id}"
            )

        secreted = parse_bool(
            row[str(schema["secreted"])]
        )

        if secreted:
            secreted_gene_ids.add(gene_id)

        family_value = row[str(schema["family"])]
        families = parse_families(family_value)

        if not families:
            genes_without_family.append(gene_id)
            continue

        if len(families) > 1:
            genes_with_multiple_families += 1

        substrate = (
            normalize_substrate(
                row[str(schema["substrate"])]
            )
            if schema["substrate"] is not None
            else "Unassigned"
        )

        direction_gene_counts[direction] += 1

        for family in families:
            cls = cazyme_class(family)

            long_rows.append(
                {
                    "species": species,
                    "gene_id": gene_id,
                    "family": family,
                    "class": cls,
                    "substrate": substrate,
                    "direction": direction,
                    "secreted": str(secreted).upper(),
                    "raw_log2FoldChange": raw_lfc,
                    "shrunk_log2FoldChange": shrunk_lfc,
                    "padj": padj,
                    "source_family_annotation": family_value,
                    "source_file": str(path),
                }
            )

            family_assignment_counts[family] += 1

    if genes_without_family:
        fail(
            f"{species}: {len(genes_without_family)} high-confidence "
            "CAZyme genes lacked a parsable recommended family. "
            "First examples: "
            + ", ".join(genes_without_family[:10])
        )

    represented_gene_ids = {
        str(row["gene_id"])
        for row in long_rows
    }

    if represented_gene_ids != gene_ids:
        fail(
            f"{species}: gene-family expansion did not retain the "
            "complete input gene set"
        )

    audit = {
        "species": species,
        "input_file": str(path),
        "input_gene_count": len(rows),
        "unique_gene_count": len(gene_ids),
        "family_assignment_rows": len(long_rows),
        "unique_families": len(family_assignment_counts),
        "genes_with_multiple_families":
            genes_with_multiple_families,
        "secreted_genes": len(secreted_gene_ids),
        "upregulated_genes":
            direction_gene_counts["up"],
        "downregulated_genes":
            direction_gene_counts["down"],
        "input_sha256": sha256(path),
        "status": "PASS",
    }

    return long_rows, audit


def summarize_family_long(
    gene_family_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in gene_family_rows:
        key = (
            str(row["species"]),
            str(row["family"]),
        )

        grouped[key].append(row)

    output: list[dict[str, object]] = []

    for (species, family), rows in grouped.items():
        up_rows = [
            row
            for row in rows
            if row["direction"] == "up"
        ]

        down_rows = [
            row
            for row in rows
            if row["direction"] == "down"
        ]

        secreted_rows = [
            row
            for row in rows
            if row["secreted"] == "TRUE"
        ]

        secreted_up = [
            row
            for row in secreted_rows
            if row["direction"] == "up"
        ]

        secreted_down = [
            row
            for row in secreted_rows
            if row["direction"] == "down"
        ]

        substrates = sorted(
            {
                str(row["substrate"])
                for row in rows
                if str(row["substrate"]) != "Unassigned"
            }
        )

        shrunk_values = [
            float(row["shrunk_log2FoldChange"])
            for row in rows
        ]

        mean_shrunk = (
            sum(shrunk_values) / len(shrunk_values)
        )

        output.append(
            {
                "species": species,
                "family": family,
                "class": cazyme_class(family),
                "substrates": (
                    " | ".join(substrates)
                    if substrates
                    else "Unassigned"
                ),
                "significant_genes": len(rows),
                "up": len(up_rows),
                "down": len(down_rows),
                "secreted": len(secreted_rows),
                "secreted_up": len(secreted_up),
                "secreted_down": len(secreted_down),
                "mean_shrunk_log2FoldChange":
                    mean_shrunk,
                "max_abs_shrunk_log2FoldChange":
                    max(abs(value) for value in shrunk_values),
            }
        )

    output.sort(
        key=lambda row: (
            family_sort_key(str(row["family"])),
            str(row["species"]),
        )
    )

    return output


def summarize_family_wide(
    family_long: list[dict[str, object]],
) -> list[dict[str, object]]:
    indexed: dict[
        tuple[str, str],
        dict[str, object],
    ] = {
        (
            str(row["family"]),
            str(row["species"]),
        ): row
        for row in family_long
    }

    families = sorted(
        {
            str(row["family"])
            for row in family_long
        },
        key=family_sort_key,
    )

    output: list[dict[str, object]] = []

    for family in families:
        f_row = indexed.get((family, "Fusarium"))
        o_row = indexed.get((family, "Ophiostoma"))

        substrates = sorted(
            {
                value
                for source in (f_row, o_row)
                if source is not None
                for value in str(
                    source["substrates"]
                ).split(" | ")
                if value != "Unassigned"
            }
        )

        row: dict[str, object] = {
            "family": family,
            "class": cazyme_class(family),
            "substrates": (
                " | ".join(substrates)
                if substrates
                else "Unassigned"
            ),
        }

        for prefix, source in (
            ("fusarium", f_row),
            ("ophiostoma", o_row),
        ):
            row[f"{prefix}_significant"] = (
                source["significant_genes"]
                if source is not None
                else 0
            )

            row[f"{prefix}_up"] = (
                source["up"]
                if source is not None
                else 0
            )

            row[f"{prefix}_down"] = (
                source["down"]
                if source is not None
                else 0
            )

            row[f"{prefix}_secreted"] = (
                source["secreted"]
                if source is not None
                else 0
            )

            row[f"{prefix}_secreted_up"] = (
                source["secreted_up"]
                if source is not None
                else 0
            )

            row[f"{prefix}_secreted_down"] = (
                source["secreted_down"]
                if source is not None
                else 0
            )

            row[f"{prefix}_mean_shrunk_lfc"] = (
                source[
                    "mean_shrunk_log2FoldChange"
                ]
                if source is not None
                else ""
            )

            row[
                f"{prefix}_max_abs_shrunk_lfc"
            ] = (
                source[
                    "max_abs_shrunk_log2FoldChange"
                ]
                if source is not None
                else ""
            )

        row["total_significant"] = (
            int(row["fusarium_significant"])
            + int(row["ophiostoma_significant"])
        )

        row["total_secreted"] = (
            int(row["fusarium_secreted"])
            + int(row["ophiostoma_secreted"])
        )

        output.append(row)

    output.sort(
        key=lambda row: (
            -int(row["total_secreted"]),
            -int(row["total_significant"]),
            family_sort_key(str(row["family"])),
        )
    )

    return output


def summarize_category(
    gene_family_rows: list[dict[str, object]],
    category_field: str,
) -> list[dict[str, object]]:
    grouped: dict[
        tuple[str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in gene_family_rows:
        key = (
            str(row["species"]),
            str(row[category_field]),
        )

        grouped[key].append(row)

    output: list[dict[str, object]] = []

    for (species, category), rows in grouped.items():
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
            if row["secreted"] == "TRUE"
        }

        secreted_up_gene_ids = {
            str(row["gene_id"])
            for row in rows
            if (
                row["secreted"] == "TRUE"
                and row["direction"] == "up"
            )
        }

        secreted_down_gene_ids = {
            str(row["gene_id"])
            for row in rows
            if (
                row["secreted"] == "TRUE"
                and row["direction"] == "down"
            )
        }

        output.append(
            {
                category_field: category,
                "species": species,
                "significant_genes": len(gene_ids),
                "up": len(up_gene_ids),
                "down": len(down_gene_ids),
                "secreted": len(secreted_gene_ids),
                "secreted_up":
                    len(secreted_up_gene_ids),
                "secreted_down":
                    len(secreted_down_gene_ids),
            }
        )

    output.sort(
        key=lambda row: (
            str(row[category_field]),
            str(row["species"]),
        )
    )

    return output


def main() -> None:
    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fusarium_rows, fusarium_audit = process_species(
        "Fusarium",
        FUSARIUM_INPUT,
    )

    ophiostoma_rows, ophiostoma_audit = process_species(
        "Ophiostoma",
        OPHIOSTOMA_INPUT,
    )

    gene_family_rows = (
        fusarium_rows
        + ophiostoma_rows
    )

    gene_family_rows.sort(
        key=lambda row: (
            family_sort_key(str(row["family"])),
            str(row["species"]),
            str(row["direction"]),
            str(row["gene_id"]),
        )
    )

    write_tsv(
        GENE_FAMILY_OUT,
        [
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
            "source_family_annotation",
            "source_file",
        ],
        gene_family_rows,
    )

    family_long = summarize_family_long(
        gene_family_rows
    )

    write_tsv(
        FAMILY_LONG_OUT,
        [
            "species",
            "family",
            "class",
            "substrates",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
            "mean_shrunk_log2FoldChange",
            "max_abs_shrunk_log2FoldChange",
        ],
        family_long,
    )

    family_wide = summarize_family_wide(
        family_long
    )

    write_tsv(
        FAMILY_WIDE_OUT,
        [
            "family",
            "class",
            "substrates",
            "fusarium_significant",
            "fusarium_up",
            "fusarium_down",
            "fusarium_secreted",
            "fusarium_secreted_up",
            "fusarium_secreted_down",
            "fusarium_mean_shrunk_lfc",
            "fusarium_max_abs_shrunk_lfc",
            "ophiostoma_significant",
            "ophiostoma_up",
            "ophiostoma_down",
            "ophiostoma_secreted",
            "ophiostoma_secreted_up",
            "ophiostoma_secreted_down",
            "ophiostoma_mean_shrunk_lfc",
            "ophiostoma_max_abs_shrunk_lfc",
            "total_significant",
            "total_secreted",
        ],
        family_wide,
    )

    class_summary = summarize_category(
        gene_family_rows,
        "class",
    )

    write_tsv(
        CLASS_OUT,
        [
            "class",
            "species",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
        ],
        class_summary,
    )

    substrate_summary = summarize_category(
        gene_family_rows,
        "substrate",
    )

    write_tsv(
        SUBSTRATE_OUT,
        [
            "substrate",
            "species",
            "significant_genes",
            "up",
            "down",
            "secreted",
            "secreted_up",
            "secreted_down",
        ],
        substrate_summary,
    )

    audits = [
        fusarium_audit,
        ophiostoma_audit,
    ]

    write_tsv(
        AUDIT_OUT,
        [
            "species",
            "input_file",
            "input_gene_count",
            "unique_gene_count",
            "family_assignment_rows",
            "unique_families",
            "genes_with_multiple_families",
            "secreted_genes",
            "upregulated_genes",
            "downregulated_genes",
            "input_sha256",
            "status",
        ],
        audits,
    )

    run_info = [
        {
            "field": "script",
            "value":
                "33_build_comparative_cazyme_family_table.py",
        },
        {
            "field": "run_timestamp_utc",
            "value":
                datetime.now(timezone.utc).isoformat(),
        },
        {
            "field": "repository_root",
            "value": str(ROOT),
        },
        {
            "field": "contrast",
            "value": "matched interaction_vs_control; "
                     "Fusarium internal contrast=interaction_vs_self; "
                     "Ophiostoma internal contrast=interaction_vs_onu",
        },
        {
            "field": "cazyme_definition",
            "value":
                "dbCAN high confidence; at least two tools",
        },
        {
            "field": "significance_definition",
            "value": "padj < 0.05",
        },
        {
            "field": "direction_definition",
            "value":
                "sign of raw log2 fold change",
        },
        {
            "field": "effect_size",
            "value":
                "apeglm-shrunken log2 fold change",
        },
        {
            "field": "family_counting",
            "value": (
                "one gene counted once per distinct "
                "recommended CAZyme family"
            ),
        },
        {
            "field": "fusarium_input_sha256",
            "value": sha256(FUSARIUM_INPUT),
        },
        {
            "field": "ophiostoma_input_sha256",
            "value": sha256(OPHIOSTOMA_INPUT),
        },
    ]

    write_tsv(
        RUN_INFO_OUT,
        ["field", "value"],
        run_info,
    )

    print()
    print("============================================================")
    print("COMPARATIVE CAZYME FAMILY TABLE COMPLETE")
    print("============================================================")
    print()

    for audit in audits:
        print(audit["species"])
        print(
            "  Significant high-confidence genes: "
            f"{audit['unique_gene_count']}"
        )
        print(
            "  Family-assignment rows: "
            f"{audit['family_assignment_rows']}"
        )
        print(
            "  Unique families: "
            f"{audit['unique_families']}"
        )
        print(
            "  Multi-family genes: "
            f"{audit['genes_with_multiple_families']}"
        )
        print(
            "  Secreted significant CAZyme genes: "
            f"{audit['secreted_genes']}"
        )
        print(
            "  Direction: "
            f"up {audit['upregulated_genes']}, "
            f"down {audit['downregulated_genes']}"
        )
        print()

    print(f"Family summary: {FAMILY_WIDE_OUT}")
    print(f"Gene-family table: {GENE_FAMILY_OUT}")
    print(f"Audit: {AUDIT_OUT}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
