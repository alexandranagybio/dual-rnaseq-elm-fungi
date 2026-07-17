# Analysis audit

This document records every analytical decision used in the manuscript pipeline.

## Core principles

- The historical thesis analysis is preserved separately and is not edited.
- All manuscript analyses are rebuilt from verified inputs.
- Thresholds are fixed before results are inspected.
- Genes, transcripts, exons, proteins, and annotation entries are reported separately.
- Every manuscript number must be reproducible from a named script and output file.
- Corrected manuscript results may differ from the thesis results.

## Differential-expression thresholds

Primary statistical threshold:

- adjusted p-value < 0.05

Effect-size threshold:

- absolute log2 fold change > 1

These thresholds are fixed across organisms and contrasts.

Dynamic significance thresholds are not permitted.

## Threshold audit

| Analysis | Historical rule | Manuscript rule | Decision | Status |
|---|---|---|---|---|
| Ophiostoma DESeq2 | To verify | padj < 0.05 | Pending | Not audited |
| Ophiostoma effect-size filter | To verify | padj < 0.05 and abs(log2FC) > 1 | Pending | Not audited |
| Ophiostoma GO | Dynamic threshold found | Fixed threshold | Remove dynamic rule | Not rerun |
| Ophiostoma KEGG | Dynamic threshold found | Fixed threshold | Remove dynamic rule | Not rerun |
| Ophiostoma dbCAN | padj < 0.05 and direction only | padj < 0.05 and abs(log2FC) > 1 | Change | Not rerun |
| Fusarium analyses | To verify | Fixed threshold | Pending | Not audited |

## Reporting units

These units must not be treated as interchangeable:

- genes
- transcripts
- exons
- proteins
- unique protein-family combinations
- CAZyme annotation entries

## Previously verified observations

- The historical Ophiostoma DESeq2 analysis reproduced in a clean R session.
- The historical count matrix contains exon-level identifiers.
- The reference annotation contains 8,640 genes or mRNAs and 17,638 exons.
- Significant proteins in the historical protein table were unique after filtering.
- Duplicate protein rows occurred only among rows with missing adjusted p-values.
- Historical downstream scripts contained inconsistent and dynamic significance rules.

## Next action

Identify and copy only verified manuscript inputs into this repository.
