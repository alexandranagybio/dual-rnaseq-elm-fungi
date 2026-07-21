#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path.cwd()
FASTA = ROOT / "data/reference/ophiostoma/ophiostoma.proteins.clean.fa"
SIGNALP = ROOT / "results/ophiostoma/signalp_raw/prediction_results.txt"
DBCAN = ROOT / "results/ophiostoma/dbcan_raw/overview.tsv"
OUTDIR = ROOT / "results/ophiostoma/functional_annotation"
OUT_TSV = OUTDIR / "ophiostoma_functional_annotation.tsv"
OUT_SUMMARY = OUTDIR / "validation_summary.txt"

EXPECTED = 8640
MRNA_RE = re.compile(r"^mRNA_(\d+)$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def fasta_ids(path: Path) -> list[str]:
    ids = [line[1:].split()[0] for line in path.read_text().splitlines()
           if line.startswith(">")]
    if len(ids) != len(set(ids)):
        fail("duplicated FASTA identifiers")
    return ids


def gene_id(mrna_id: str) -> str:
    match = MRNA_RE.fullmatch(mrna_id)
    if not match:
        fail(f"unexpected FASTA ID: {mrna_id}")
    return f"gene_{match.group(1)}"


def read_signalp(path: Path) -> dict[str, dict[str, str]]:
    out = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            pid, prediction, other, sp = row[:4]
            cs = row[4] if len(row) > 4 else ""
            if pid in out:
                fail(f"duplicated SignalP ID: {pid}")
            if prediction not in {"OTHER", "SP"}:
                fail(f"unexpected SignalP prediction: {prediction}")
            out[pid] = {
                "signalp_prediction": prediction,
                "signalp_is_sp": str(prediction == "SP").upper(),
                "signalp_other_score": other,
                "signalp_sp_score": sp,
                "signalp_cs_position": cs,
            }
    return out


def read_dbcan(path: Path) -> dict[str, dict[str, str]]:
    out = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pid = row["Gene ID"]
            if pid in out:
                fail(f"duplicated dbCAN ID: {pid}")
            n = int(row["#ofTools"])
            out[pid] = {
                "dbcan_any_hit": "TRUE",
                "dbcan_high_confidence": str(n >= 2).upper(),
                "dbcan_n_tools": str(n),
                "dbcan_ec": row["EC#"],
                "dbcan_hmm": row["dbCAN_hmm"],
                "dbcan_sub": row["dbCAN_sub"],
                "dbcan_diamond": row["DIAMOND"],
                "dbcan_recommended": row["Recommend Results"],
                "dbcan_substrate": row["Substrate"],
            }
    return out


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


for path in (FASTA, SIGNALP, DBCAN):
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing or empty file: {path}")

ids = fasta_ids(FASTA)
signalp = read_signalp(SIGNALP)
dbcan = read_dbcan(DBCAN)

if len(ids) != EXPECTED:
    fail(f"expected {EXPECTED} proteins, observed {len(ids)}")
if set(signalp) != set(ids):
    fail("SignalP IDs do not exactly match FASTA IDs")
if not set(dbcan).issubset(set(ids)):
    fail("dbCAN contains IDs absent from FASTA")

rows = []
for pid in sorted(ids, key=lambda x: int(MRNA_RE.fullmatch(x).group(1))):
    row = {"gene_id": gene_id(pid), "mrna_id": pid}
    row.update(signalp[pid])
    row.update(dbcan.get(pid, blank_dbcan()))
    rows.append(row)

OUTDIR.mkdir(parents=True, exist_ok=True)
fieldnames = list(rows[0])
with OUT_TSV.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

sp_count = sum(r["signalp_is_sp"] == "TRUE" for r in rows)
any_count = sum(r["dbcan_any_hit"] == "TRUE" for r in rows)
hc_count = sum(r["dbcan_high_confidence"] == "TRUE" for r in rows)
support = Counter(r["dbcan_n_tools"] for r in rows)

summary = [
    "Ophiostoma functional annotation validation",
    "",
    f"Authoritative FASTA proteins: {len(ids)}",
    f"SignalP predictions: {len(signalp)}",
    f"SignalP SP predictions: {sp_count}",
    f"SignalP OTHER predictions: {len(rows) - sp_count}",
    f"dbCAN proteins with >=1 method hit: {any_count}",
    f"dbCAN high-confidence proteins (#ofTools >=2): {hc_count}",
    f"Proteins with no dbCAN hit: {len(rows) - any_count}",
    "",
    "dbCAN support distribution across complete proteome:",
]
for n in sorted(support, key=int):
    summary.append(f"  {n} tool(s): {support[n]}")
summary += ["", "Validation status: PASS", f"Output: {OUT_TSV}"]

OUT_SUMMARY.write_text("\n".join(summary) + "\n")
print("\n".join(summary))
