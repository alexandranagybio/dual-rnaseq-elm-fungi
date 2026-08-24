#!/usr/bin/env python3

"""
Build the canonical manuscript-facing results source table.

This script does NOT perform statistical analyses.
It extracts manuscript-facing values from previously audited canonical outputs.

Outputs
-------
results/publication/manuscript_source/MANUSCRIPT_RESULTS_SOURCE.tsv
results/publication/manuscript_source/MANUSCRIPT_RESULTS_SOURCE_run_info.tsv

Design principles
-----------------
- No manuscript number is hard-coded if it can be extracted from a canonical table.
- Hard-coded expected values are used only as explicit audit locks.
- Biological labels are publication-facing:
    Ophiostoma interaction = reaction zone
    Ophiostoma self        = non-contact region
    Ophiostoma onu         = control colony
- Cross-species comparison uses:
    Fusarium interaction_vs_self
    Ophiostoma interaction_vs_onu
- Percentages are rounded only for manuscript-facing display.
- Every canonical input is SHA256-hashed.
- The generating script itself is SHA256-hashed.
- Source row counts, Python version, and Git commit (when available) are recorded.
- Missing, duplicated, unexpected, or biologically inconsistent rows cause failure.
- The manuscript-facing result manifest has a locked expected row count.
"""

from __future__ import annotations

import csv
import hashlib
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


# ==========================================================================
# Paths and fixed audit configuration
# ==========================================================================

ROOT = Path(__file__).resolve().parents[2]

MANUSCRIPT_SOURCE_DIR = (
    ROOT / "results/publication/manuscript_source"
)

OUT = (
    MANUSCRIPT_SOURCE_DIR
    / "MANUSCRIPT_RESULTS_SOURCE.tsv"
)

RUN_INFO = (
    MANUSCRIPT_SOURCE_DIR
    / "MANUSCRIPT_RESULTS_SOURCE_run_info.tsv"
)

EXPECTED_RESULT_COUNT = 43


SOURCES = {
    "fus_pca":
        ROOT / "results/publication/figure3/pca/fusarium_pca_variance.tsv",

    "oph_pca":
        ROOT / "results/publication/figure3/pca/ophiostoma_pca_variance.tsv",

    "deg":
        ROOT / "results/publication/figure3/extent_direction/"
        "figure_3c_extent_direction_differential_expression_data.tsv",

    "cog_coverage":
        ROOT / "results/publication/cog_enrichment/cog_annotation_coverage.tsv",

    "cog_enrichment":
        ROOT / "results/publication/cog_enrichment/"
        "cog_enrichment_significant.tsv",

    "cog_interaction":
        ROOT / "results/publication/cog_species_interaction/"
        "species_cog_interaction_significant.tsv",

    "extracellular":
        ROOT / "results/publication/extracellular_response/"
        "extracellular_summary.tsv",

    "cazyme_audit":
        ROOT / "results/publication/cazyme_comparison/"
        "comparative_cazyme_audit.tsv",

    "spatial_summary":
        ROOT / "results/publication/ophiostoma_spatial_heatmap/"
        "ophiostoma_spatial_response_summary.tsv",

    "spatial_function":
        ROOT / "results/publication/ophiostoma_spatial_functional_enrichment/"
        "spatial_secretome_cazyme_enrichment.tsv",

    "spatial_run_info":
        ROOT / "results/publication/ophiostoma_spatial_functional_enrichment/"
        "run_info.tsv",
}


FIELDS = [
    "section",
    "result_id",
    "organism_or_group",
    "contrast_or_context",
    "metric",
    "value",
    "unit",
    "canonical_source",
    "source_sha256",
    "reporting_note",
]


# ==========================================================================
# Utilities
# ==========================================================================

def fail(message: str) -> None:
    raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        fail(f"Missing canonical source: {relative(path)}")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        fail(f"Canonical source is empty: {relative(path)}")

    return rows


def one(rows, description: str, **conditions):
    matches = [
        row for row in rows
        if all(
            row.get(key) == str(value)
            for key, value in conditions.items()
        )
    ]

    if len(matches) != 1:
        fail(
            f"{description}: expected exactly 1 matching row; "
            f"observed {len(matches)}. Conditions={conditions}"
        )

    return matches[0]


def fmt1(value: str | float) -> str:
    return f"{float(value):.1f}"


def fmt2(value: str | float) -> str:
    return f"{float(value):.2f}"


def fmt3(value: str | float) -> str:
    return f"{float(value):.3f}"


def fmt_sci(value: str | float) -> str:
    return f"{float(value):.3g}"


def git_commit() -> str:
    """
    Return the current Git commit when ROOT is inside a Git repository.
    This is provenance metadata only; inability to recover a commit does
    not invalidate the manuscript results source table.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def git_status() -> str:
    """
    Record whether the working tree was clean at generation time.
    """
    try:
        result = subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        )
        return "clean" if not result.stdout.strip() else "modified"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


# ==========================================================================
# Validate canonical inputs and compute provenance
# ==========================================================================

for name, path in SOURCES.items():
    if not path.exists():
        fail(
            f"Required source missing [{name}]: "
            f"{relative(path)}"
        )


HASHES = {
    key: sha256(path)
    for key, path in SOURCES.items()
}


SOURCE_ROWS = {
    key: len(read_tsv(path))
    for key, path in SOURCES.items()
}


SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_SHA256 = sha256(SCRIPT_PATH)


# ==========================================================================
# Manuscript output collector
# ==========================================================================

rows_out: list[dict[str, str]] = []


def add(
    section: str,
    result_id: str,
    organism: str,
    context: str,
    metric: str,
    value,
    unit: str,
    source_key: str,
    note: str,
) -> None:

    source = SOURCES[source_key]

    rows_out.append({
        "section": section,
        "result_id": result_id,
        "organism_or_group": organism,
        "contrast_or_context": context,
        "metric": metric,
        "value": str(value),
        "unit": unit,
        "canonical_source": relative(source),
        "source_sha256": HASHES[source_key],
        "reporting_note": note,
    })


# ==========================================================================
# Figure 3 — PCA
# ==========================================================================

fus_pca = read_tsv(SOURCES["fus_pca"])
oph_pca = read_tsv(SOURCES["oph_pca"])

for component, result_id, note in [
    ("PC1", "F3_PCA_FUS_PC1", "Report in Results"),
    ("PC2", "F3_PCA_FUS_PC2", "Usually caption/supporting text"),
]:
    r = one(
        fus_pca,
        f"Fusarium {component}",
        component=component,
    )

    add(
        "Figure3",
        result_id,
        "Fusarium",
        "confrontation vs control",
        f"{component} variance",
        fmt1(r["percent_variance"]),
        "percent",
        "fus_pca",
        note,
    )


for component, result_id in [
    ("PC1", "F3_PCA_OPH_PC1"),
    ("PC2", "F3_PCA_OPH_PC2"),
]:
    r = one(
        oph_pca,
        f"Ophiostoma {component}",
        component=component,
    )

    add(
        "Figure3",
        result_id,
        "Ophiostoma",
        "three-condition design",
        f"{component} variance",
        fmt1(r["percent_variance"]),
        "percent",
        "oph_pca",
        "Report in Results",
    )


# ==========================================================================
# Figure 3 — DEG extent
# ==========================================================================

deg = read_tsv(SOURCES["deg"])

# Biological mapping audit lock:
# Ophiostoma interaction = reaction zone
# Ophiostoma self        = non-contact region
# Ophiostoma onu         = control colony

fus = one(
    deg,
    "Fusarium confrontation vs control DEG row",
    organism="Fusarium",
    contrast="interaction_vs_self",
)

oph_control = one(
    deg,
    "Ophiostoma reaction zone vs control DEG row",
    organism="Ophiostoma",
    contrast="interaction_vs_onu",
)

oph_reaction_noncontact = one(
    deg,
    "Ophiostoma reaction zone vs non-contact DEG row",
    organism="Ophiostoma",
    contrast="interaction_vs_self",
)

oph_control_noncontact = one(
    deg,
    "Ophiostoma control vs non-contact DEG row",
    organism="Ophiostoma",
    contrast="onu_vs_self",
)


def add_deg_rows(
    row,
    prefix,
    organism,
    context,
    primary_note,
):
    add(
        "Figure3",
        f"{prefix}_TOTAL",
        organism,
        context,
        "significant DE genes",
        row["significant_total"],
        "genes",
        "deg",
        primary_note,
    )

    add(
        "Figure3",
        f"{prefix}_PERCENT",
        organism,
        context,
        "significant DE fraction",
        fmt1(row["percent_total"]),
        "percent",
        "deg",
        primary_note,
    )


add_deg_rows(
    fus,
    "F3_DE_FUS",
    "Fusarium",
    "confrontation vs control",
    "Primary cross-species comparison",
)

add(
    "Figure3",
    "F3_DE_FUS_UP",
    "Fusarium",
    "confrontation vs control",
    "induced genes",
    fus["significant_up"],
    "genes",
    "deg",
    "Primary",
)

add(
    "Figure3",
    "F3_DE_FUS_DOWN",
    "Fusarium",
    "confrontation vs control",
    "repressed genes",
    fus["significant_down"],
    "genes",
    "deg",
    "Primary",
)


add_deg_rows(
    oph_control,
    "F3_DE_OPH_CONTROL",
    "Ophiostoma",
    "reaction zone vs control",
    "Primary cross-species comparison",
)

add(
    "Figure3",
    "F3_DE_OPH_CONTROL_UP",
    "Ophiostoma",
    "reaction zone vs control",
    "induced genes",
    oph_control["significant_up"],
    "genes",
    "deg",
    "Primary",
)

add(
    "Figure3",
    "F3_DE_OPH_CONTROL_DOWN",
    "Ophiostoma",
    "reaction zone vs control",
    "repressed genes",
    oph_control["significant_down"],
    "genes",
    "deg",
    "Primary",
)


add_deg_rows(
    oph_reaction_noncontact,
    "F3_DE_OPH_NONCONTACT",
    "Ophiostoma",
    "reaction zone vs non-contact region",
    (
        "Spatial-design result; do not use as "
        "cross-species control comparison"
    ),
)

add_deg_rows(
    oph_control_noncontact,
    "F3_DE_OPH_CONTROL_NONCONTACT",
    "Ophiostoma",
    "control vs non-contact region",
    "Spatial-design result",
)


# ==========================================================================
# Figure 4 — COG coverage
# ==========================================================================

coverage = read_tsv(SOURCES["cog_coverage"])

coverage_specs = [
    ("Fusarium", "Induced", "F4_COG_FUS_INDUCED_COVERAGE"),
    ("Fusarium", "Repressed", "F4_COG_FUS_REPRESSED_COVERAGE"),
    ("Ophiostoma", "Induced", "F4_COG_OPH_INDUCED_COVERAGE"),
    ("Ophiostoma", "Repressed", "F4_COG_OPH_REPRESSED_COVERAGE"),
]

for organism, direction, result_id in coverage_specs:

    r = one(
        coverage,
        f"{organism} {direction} significant COG coverage",
        organism=organism,
        direction=direction,
        significant="TRUE",
    )

    note = (
        f"{r['genes_with_cog']} of {r['genes']} genes"
    )

    add(
        "Figure4",
        result_id,
        organism,
        f"{direction.lower()} significant genes",
        "COG annotation coverage",
        fmt1(r["cog_coverage_percent"]),
        "percent",
        "cog_coverage",
        note,
    )


# ==========================================================================
# Figure 4 — COG enrichment summary
# ==========================================================================

cog = read_tsv(SOURCES["cog_enrichment"])

# Audit locks: expected significant COG sets.
# These are assertions against the canonical enrichment table, not
# independent manuscript values. Any change upstream causes failure.
expected_cog_sets = {
    ("Fusarium", "Induced"): {"G"},
    ("Fusarium", "Repressed"): {"J", "C", "E", "H"},
    ("Ophiostoma", "Induced"): {"I", "C", "O", "Z", "U"},
    ("Ophiostoma", "Repressed"): {"S", "P", "K", "B", "Q", "G"},
}

for key, expected in expected_cog_sets.items():
    organism, direction = key

    observed = {
        row["cog"]
        for row in cog
        if row["organism"] == organism
        and row["direction"] == direction
    }

    if observed != expected:
        fail(
            f"Unexpected COG set for {organism} {direction}. "
            f"Expected={sorted(expected)}, "
            f"observed={sorted(observed)}"
        )


fus_cog_value = (
    "G induced; "
    + ",".join(
        sorted(
            expected_cog_sets[
                ("Fusarium", "Repressed")
            ]
        )
    )
    + " repressed"
)

oph_cog_value = (
    ",".join(
        sorted(
            expected_cog_sets[
                ("Ophiostoma", "Induced")
            ]
        )
    )
    + " induced; "
    + ",".join(
        sorted(
            expected_cog_sets[
                ("Ophiostoma", "Repressed")
            ]
        )
    )
    + " repressed"
)

add(
    "Figure4",
    "F4_COG_FUS_ENRICHED",
    "Fusarium",
    "confrontation vs control",
    "significant enriched COG categories",
    fus_cog_value,
    "categories",
    "cog_enrichment",
    "Functional summary",
)

add(
    "Figure4",
    "F4_COG_OPH_ENRICHED",
    "Ophiostoma",
    "reaction zone vs control",
    "significant enriched COG categories",
    oph_cog_value,
    "categories",
    "cog_enrichment",
    "Functional summary",
)


interaction_rows = read_tsv(SOURCES["cog_interaction"])

EXPECTED_COG_INTERACTION_COUNT = 23

if len(interaction_rows) != EXPECTED_COG_INTERACTION_COUNT:
    fail(
        "Unexpected number of significant species × COG interactions: "
        f"expected {EXPECTED_COG_INTERACTION_COUNT}, "
        f"observed {len(interaction_rows)}"
    )

add(
    "Figure4",
    "F4_COG_INTERACTIONS",
    "Both",
    "matched confrontation vs control",
    "significant species x COG interactions",
    len(interaction_rows),
    "interactions",
    "cog_interaction",
    "Report selected biologically informative examples only",
)


# ==========================================================================
# Figure 5 — extracellular response
# ==========================================================================

extra = read_tsv(SOURCES["extracellular"])

extra_fus = one(
    extra,
    "Fusarium extracellular summary",
    species="Fusarium",
)

extra_oph = one(
    extra,
    "Ophiostoma extracellular summary",
    species="Ophiostoma",
)

for row, species, context, prefix in [
    (
        extra_fus,
        "Fusarium",
        "confrontation vs control",
        "FUS",
    ),
    (
        extra_oph,
        "Ophiostoma",
        "reaction zone vs control",
        "OPH",
    ),
]:
    add(
        "Figure5",
        f"F5_SECRETED_{prefix}",
        species,
        context,
        "significant SignalP-positive genes",
        row["significant_secreted_genes"],
        "genes",
        "extracellular",
        "Classical secretion candidates",
    )

    add(
        "Figure5",
        f"F5_SECRETED_CAZYME_{prefix}",
        species,
        context,
        "significant secreted CAZyme genes",
        row["significant_secreted_cazymes"],
        "genes",
        "extracellular",
        (
            "Subset of significant CAZymes and "
            "SignalP-positive genes"
        ),
    )


cazy = read_tsv(SOURCES["cazyme_audit"])

cazy_fus = one(
    cazy,
    "Fusarium CAZyme audit",
    species="Fusarium",
)

cazy_oph = one(
    cazy,
    "Ophiostoma CAZyme audit",
    species="Ophiostoma",
)

for row, species, context, prefix in [
    (
        cazy_fus,
        "Fusarium",
        "confrontation vs control",
        "FUS",
    ),
    (
        cazy_oph,
        "Ophiostoma",
        "reaction zone vs control",
        "OPH",
    ),
]:
    if row["status"] != "PASS":
        fail(
            f"{species} CAZyme audit status is not PASS"
        )

    add(
        "Figure5",
        f"F5_CAZYME_{prefix}",
        species,
        context,
        "significant high-confidence CAZyme genes",
        row["unique_gene_count"],
        "genes",
        "cazyme_audit",
        "dbCAN >=2 tools",
    )

    add(
        "Figure5",
        f"F5_CAZYME_DIRECTION_{prefix}",
        species,
        context,
        "CAZyme direction",
        (
            f"{row['upregulated_genes']} induced; "
            f"{row['downregulated_genes']} repressed"
        ),
        "genes",
        "cazyme_audit",
        "Unique genes",
    )


# ==========================================================================
# Figure 6 — spatial response programs
# ==========================================================================

spatial = read_tsv(SOURCES["spatial_summary"])

spatial_specs = [
    (
        "Reaction-zone specific",
        "F6_REACTION_SPECIFIC",
    ),
    (
        "Plate-wide confrontation response",
        "F6_PLATE_WIDE",
    ),
    (
        "Complex spatial response",
        "F6_COMPLEX",
    ),
    (
        "Non-contact-region specific",
        "F6_NONCONTACT_SPECIFIC",
    ),
]

classified_total = 0

for group, result_id in spatial_specs:

    r = one(
        spatial,
        f"Spatial response group {group}",
        response_group=group,
    )

    available = int(r["available_genes"])
    classified_total += available

    add(
        "Figure6",
        result_id,
        "Ophiostoma",
        "spatial response",
        f"{group} genes",
        available,
        "genes",
        "spatial_summary",
        "Canonical spatial class",
    )


run_info_rows = read_tsv(SOURCES["spatial_run_info"])

run_info_lookup = {
    row["field"]: row["value"]
    for row in run_info_rows
}

if "classified_gene_count" not in run_info_lookup:
    fail(
        "Spatial run_info is missing classified_gene_count"
    )

reported_classified = int(
    run_info_lookup["classified_gene_count"]
)

if classified_total != reported_classified:
    fail(
        "Spatial-class sum does not equal run_info "
        "classified_gene_count: "
        f"{classified_total} != {reported_classified}"
    )

add(
    "Figure6",
    "F6_CLASSIFIED_TOTAL",
    "Ophiostoma",
    "spatial response",
    "classified genes",
    classified_total,
    "genes",
    "spatial_run_info",
    "Sum of four canonical spatial classes",
)


# ==========================================================================
# Figure 6 — secretome / CAZyme enrichment
# ==========================================================================

spatial_fun = read_tsv(SOURCES["spatial_function"])


def add_spatial_enrichment(
    result_id,
    group,
    feature,
    metric,
):
    r = one(
        spatial_fun,
        f"{group} / {feature}",
        response_group=group,
        feature=feature,
    )

    significant = r["significant"] == "TRUE"

    if significant:
        value = r["enrichment_direction"].lower()
    else:
        value = "not significant"

    note = (
        f"{r['group_positive']} of {r['group_total']} genes; "
        f"{fmt2(r['group_percent'])}%; "
        f"OR {fmt3(r['odds_ratio'])}; "
        f"padj {fmt_sci(r['padj'])}"
    )

    add(
        "Figure6",
        result_id,
        "Ophiostoma",
        group.lower(),
        metric,
        value,
        "result",
        "spatial_function",
        note,
    )


# Complex group: retain counts as manuscript-facing values.
complex_sec = one(
    spatial_fun,
    "Complex spatial / secreted",
    response_group="Complex spatial response",
    feature="Secreted protein",
)

complex_caz = one(
    spatial_fun,
    "Complex spatial / CAZyme",
    response_group="Complex spatial response",
    feature="High-confidence CAZyme",
)

# Audit lock for the Figure 6 biological interpretation:
# the complex spatial-response class must remain significantly enriched
# for both secreted proteins and high-confidence CAZymes.
for r, description in [
    (
        complex_sec,
        "Complex spatial / secreted",
    ),
    (
        complex_caz,
        "Complex spatial / CAZyme",
    ),
]:
    if r["significant"] != "TRUE":
        fail(
            f"{description} is no longer statistically significant"
        )

    if r["enrichment_direction"].lower() != "enriched":
        fail(
            f"{description} is no longer enriched"
        )


add(
    "Figure6",
    "F6_COMPLEX_SECRETED",
    "Ophiostoma",
    "complex spatial response",
    "secreted proteins",
    complex_sec["group_positive"],
    "genes",
    "spatial_function",
    (
        f"{fmt2(complex_sec['group_percent'])}%; "
        f"OR {fmt3(complex_sec['odds_ratio'])}; "
        f"padj {fmt_sci(complex_sec['padj'])}"
    ),
)

add(
    "Figure6",
    "F6_COMPLEX_CAZYME",
    "Ophiostoma",
    "complex spatial response",
    "high-confidence CAZymes",
    complex_caz["group_positive"],
    "genes",
    "spatial_function",
    (
        f"{fmt2(complex_caz['group_percent'])}%; "
        f"OR {fmt3(complex_caz['odds_ratio'])}; "
        f"padj {fmt_sci(complex_caz['padj'])}"
    ),
)


add_spatial_enrichment(
    "F6_REACTION_SECRETED_ENRICHMENT",
    "Reaction-zone specific",
    "Secreted protein",
    "secreted-protein enrichment",
)

add_spatial_enrichment(
    "F6_REACTION_CAZYME_ENRICHMENT",
    "Reaction-zone specific",
    "High-confidence CAZyme",
    "CAZyme enrichment",
)

add_spatial_enrichment(
    "F6_PLATE_WIDE_CAZYME",
    "Plate-wide confrontation response",
    "High-confidence CAZyme",
    "CAZyme enrichment",
)

add_spatial_enrichment(
    "F6_NONCONTACT_SECRETED",
    "Non-contact-region specific",
    "Secreted protein",
    "secreted-protein enrichment",
)

add_spatial_enrichment(
    "F6_NONCONTACT_CAZYME",
    "Non-contact-region specific",
    "High-confidence CAZyme",
    "CAZyme enrichment",
)


# ==========================================================================
# Final validation
# ==========================================================================

result_ids = [
    row["result_id"]
    for row in rows_out
]

if len(result_ids) != len(set(result_ids)):
    fail(
        "Duplicate manuscript result_id detected"
    )

if len(rows_out) != EXPECTED_RESULT_COUNT:
    fail(
        "Unexpected manuscript result count: "
        f"expected {EXPECTED_RESULT_COUNT}, "
        f"observed {len(rows_out)}"
    )


# Cross-table invariants.
if int(fus["significant_total"]) != int(
    extra_fus["significant_de_genes"]
):
    fail(
        "Fusarium DEG count disagrees with "
        "extracellular summary"
    )

if int(oph_control["significant_total"]) != int(
    extra_oph["significant_de_genes"]
):
    fail(
        "Ophiostoma DEG count disagrees with "
        "extracellular summary"
    )

if int(cazy_fus["unique_gene_count"]) != int(
    extra_fus["significant_high_confidence_cazymes"]
):
    fail(
        "Fusarium CAZyme count disagrees "
        "across canonical tables"
    )

if int(cazy_oph["unique_gene_count"]) != int(
    extra_oph["significant_high_confidence_cazymes"]
):
    fail(
        "Ophiostoma CAZyme count disagrees "
        "across canonical tables"
    )

if int(cazy_fus["secreted_genes"]) != int(
    extra_fus["significant_secreted_cazymes"]
):
    fail(
        "Fusarium secreted CAZyme count disagrees "
        "across canonical tables"
    )

if int(cazy_oph["secreted_genes"]) != int(
    extra_oph["significant_secreted_cazymes"]
):
    fail(
        "Ophiostoma secreted CAZyme count disagrees "
        "across canonical tables"
    )


# ==========================================================================
# Write manuscript source table
# ==========================================================================

MANUSCRIPT_SOURCE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

with OUT.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows_out)


# ==========================================================================
# Write provenance / run info
# ==========================================================================

run_info = [
    (
        "script",
        relative(SCRIPT_PATH),
    ),
    (
        "script_sha256",
        SCRIPT_SHA256,
    ),
    (
        "python_version",
        platform.python_version(),
    ),
    (
        "git_commit",
        git_commit(),
    ),
    (
        "git_working_tree",
        git_status(),
    ),
    (
        "run_timestamp_utc",
        datetime.now(timezone.utc).isoformat(),
    ),
    (
        "output",
        relative(OUT),
    ),
    (
        "output_sha256",
        sha256(OUT),
    ),
    (
        "row_count",
        str(len(rows_out)),
    ),
    (
        "expected_row_count",
        str(EXPECTED_RESULT_COUNT),
    ),
    (
        "biological_mapping",
        (
            "Ophiostoma: interaction=reaction zone; "
            "self=non-contact region; "
            "onu=control colony"
        ),
    ),
    (
        "cross_species_contrast",
        (
            "Fusarium interaction_vs_self; "
            "Ophiostoma interaction_vs_onu"
        ),
    ),
    (
        "significant_de_definition",
        "padj < 0.05",
    ),
    (
        "cazyme_definition",
        "high-confidence dbCAN; >=2 tools",
    ),
]

for key, path in SOURCES.items():
    run_info.append(
        (
            f"source_{key}",
            relative(path),
        )
    )
    run_info.append(
        (
            f"source_rows_{key}",
            str(SOURCE_ROWS[key]),
        )
    )
    run_info.append(
        (
            f"sha256_{key}",
            HASHES[key],
        )
    )


with RUN_INFO.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        ["field", "value"]
    )
    writer.writerows(run_info)


print("=" * 60)
print("MANUSCRIPT RESULTS SOURCE COMPLETE")
print("=" * 60)
print()
print(f"Rows: {len(rows_out)}")
print(f"Output:   {relative(OUT)}")
print(f"Run info: {relative(RUN_INFO)}")
print(f"Output SHA256: {sha256(OUT)}")
print()
print(
    "PASS: canonical values extracted, provenance recorded, "
    "and cross-table invariants validated."
)
