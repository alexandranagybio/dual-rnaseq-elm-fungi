# Ophiostoma–Fusarium comparability checklist

Use this checklist before interpreting or comparing the Fusarium and Ophiostoma datasets.

The two organisms do **not** need identical upstream quantification software when their
reference resources differ. Ophiostoma uses reference-guided gene counting, whereas
Fusarium currently uses a de novo transcriptome and Salmon/tximport. They do need the
same biological definitions, statistical thresholds, annotation rules, and reporting
logic at the final gene-level comparison stage.

## Locked definitions from the validated Ophiostoma workflow

| Analysis component | Validated Ophiostoma definition | Requirement for Fusarium comparison |
|---|---|---|
| Statistical unit | Gene-level DESeq2 analysis | Aggregate transcript estimates to a validated gene-level unit before comparison; do not compare Ophiostoma genes with uncollapsed Fusarium isoforms |
| DESeq2 design | `~ condition` | Use the equivalent condition-only design unless a documented design difference is biologically required |
| DESeq2 condition coding / biological reference | Internal factor labels must be mapped explicitly to biological sample classes. For Ophiostoma, `interaction` = reaction zone, `self` = non-contact region of the confrontation plate, and `onu` = control colony | Define the biologically corresponding Fusarium conditions explicitly and document every contrast numerator and denominator; do not assume that internal factor labels denote equivalent biological conditions across organisms |
| DE test | DESeq2 Wald test | Use the same test unless a justified model change is documented |
| Multiple-testing threshold | `padj < 0.05` | Use `padj < 0.05` |
| “Significant DE” | `padj < 0.05` | Use exactly the same definition |
| “Strong DE” | `padj < 0.05` and `abs(log2FoldChange) > 1` | Use exactly the same descriptive subset |
| Fold-change boundary | Strictly `> 1`, not `>= 1` | Preserve the strict boundary |
| Independent filtering | Enabled in DESeq2 | Keep enabled and report tested/unfiltered gene counts |
| All-zero genes | Removed before inference | Remove all-zero gene rows before DESeq2 and report how many were removed |
| PCA transformation | Blind VST | Use blind VST for comparable unsupervised PCA |
| PCA feature selection | Top 500 most variable genes | Use the same number unless Fusarium has fewer eligible genes; document any deviation |
| SignalP software | SignalP 6.0 fast model | Run the same SignalP generation/model mode where technically possible |
| SignalP-positive definition | `signalp_is_sp == TRUE`, corresponding to classical Sec/SPI signal peptide prediction | Use exactly the same positive class |
| “Secretome” wording | SignalP-positive candidates for classical secretion; not a complete membrane-filtered secretome | Use the same cautious terminology |
| dbCAN software | dbCAN 5.2.8 | Use the same version/database snapshot where possible |
| Any dbCAN annotation | At least one dbCAN method hit; `dbcan_any_hit == TRUE` | Preserve as a secondary, broad annotation set |
| High-confidence CAZyme | At least two dbCAN methods; `dbcan_high_confidence == TRUE`; `#ofTools >= 2` | Use as the primary CAZyme definition |
| GO enrichment universe | All DESeq2-tested genes with valid gene-level mapping; ontology-specific annotated universe used by enrichment | Build the Fusarium universe from its own tested genes, not from all assembled transcripts or only significant genes |
| GO selected sets | Significant and strong DE sets, separated into up and down directions and BP/MF/CC ontologies | Use the same gene sets, directions, and ontologies |
| GO significance | Enrichment-adjusted `p.adjust < 0.05` | Use the same threshold and multiple-testing procedure |
| KEGG enrichment universe | All DESeq2-tested genes with KEGG annotation | Use the Fusarium tested gene universe with valid KEGG mapping |
| KEGG selected sets | Significant and strong DE sets, separated into up and down directions | Use the same definitions |
| KEGG significance | Enrichment-adjusted `p.adjust < 0.05` | Use the same threshold |
| Annotation denominator | Report complete annotation set separately from DESeq2-tested set | Always distinguish complete Fusarium annotation, tested genes, significant genes, and strong genes |
| Direction | Positive log2FC = up in the named numerator; negative = down | Verify every contrast orientation before comparing counts |
| Reproducibility | Save scripts, versions, hashes, validation tables, and run metadata | Mirror the same audit structure for Fusarium |

## Manuscript-level functional comparison

For the current manuscript, GO and KEGG analyses are retained as validated
exploratory/supporting analyses but are not used as the primary cross-species
functional comparison.

The primary manuscript-level functional comparison uses gene-level COG
enrichment with organism-specific tested-gene backgrounds. Differences in
functional response between Fusarium and Ophiostoma are additionally evaluated
using the formal species × COG interaction analysis.

The same comparability principle applies: upstream organism-specific resources
may differ, but biological definitions, gene-level statistical units,
significance thresholds, annotation logic, contrast orientation, and reporting
rules must be explicit and reproducible.

## Required denominator reporting for each organism

Always report these as separate quantities:

1. Complete gene/protein annotation universe
2. Genes retained for DESeq2 testing
3. SignalP-positive genes in the complete and tested universes
4. Any-hit and high-confidence CAZymes in the complete and tested universes
5. Significant DE genes
6. Strong DE genes
7. Significant and strong SignalP-positive subsets
8. Significant and strong high-confidence CAZyme subsets
9. GO-annotated and KEGG-annotated tested-gene universes

## Critical warning for the Fusarium de novo assembly

A Trinity “gene” and a Trinity transcript/isoform are not interchangeable. Before
comparison, confirm that:

- Salmon estimates are imported with tximport using a validated transcript-to-gene map.
- Differential expression is performed at Trinity gene level, not transcript-isoform level.
- SignalP and dbCAN protein predictions map uniquely and reproducibly back to the same
  gene-level identifiers used by DESeq2.
- When several proteins or isoforms map to one Trinity gene, the rule for collapsing
  annotation to gene level is defined before examining biological results.
- Counts and percentages are compared using organism-specific tested-gene denominators;
  raw numbers alone are not directly comparable because the annotation universes differ.

## Interpretation rule

Matching thresholds make the analyses methodologically comparable, but they do not make
raw gene counts directly equivalent. Compare proportions, directions, functional
patterns, and enrichment results alongside absolute counts, while retaining each
organism's own valid tested-gene universe.
