# Fungal confrontation reveals distinct and spatially structured transcriptional responses in the elm endophyte *Fusarium salinense* and the Dutch elm disease pathogen *Ophiostoma novo-ulmi*

## Overview

This repository contains the computational workflow and validated analysis outputs used for a dual RNA-seq study of fungal confrontation between the Dutch elm disease pathogen *Ophiostoma novo-ulmi* and the fungal endophyte *Fusarium salinense*.

The two fungi were analyzed using complementary transcriptomic strategies:

* *Ophiostoma novo-ulmi*: reference-guided alignment and gene-level quantification
* *Fusarium salinense*: de novo transcriptome assembly followed by transcript quantification and gene-level analysis

The repository focuses on the reproducible manuscript-analysis workflow from validated processed inputs. Large sequencing files and computationally expensive upstream products are not stored in Git. Instead, repository-relative paths under `data/external/` provide local access to these inputs when available.

Downstream analyses include differential expression, functional annotation, Gene Ontology and KEGG enrichment, secretome and CAZyme analyses, comparative functional analyses, spatial transcriptional responses, and generation of publication figures and manuscript-facing result tables.

---

## Experimental design

Dual RNA-seq libraries were generated from fungal interaction cultures and species-specific controls, with biological triplicates for each experimental condition.

For *O. novo-ulmi*, three sampled conditions were analyzed:

* reaction zone during interaction
* non-contact region of the interacting colony
* control colony

For *Fusarium salinense*, confrontation samples were compared with the corresponding control condition.

The authoritative sample definitions and biological contrast mappings used by the computational workflow are stored in:

```text
config/samples.tsv
config/biological_contrasts.tsv
```

---

## Repository structure

```text
.
├── config/
│   ├── samples.tsv
│   └── biological_contrasts.tsv
│
├── data/
│   ├── annotation/
│   ├── external/
│   ├── metadata/
│   └── reference/
│
├── figures/
│   ├── figure4_cog_enrichment/
│   ├── figure5_secretome_cazymes/
│   ├── figure6_ophiostoma_spatial_response/
│   └── publication/
│
├── results/
│   ├── fusarium/
│   ├── ophiostoma/
│   └── publication/
│
├── workflow/
│   ├── rules/
│   ├── scripts/
│   └── setup/
│
├── CHANGELOG.md
├── CITATION.cff
├── environment.yml
├── LICENSE
└── README.md
```

---

## Data organization

### Version-controlled data

Small reference, annotation, metadata, validation, and derived analysis files required by the manuscript workflow are stored under `data/`, `results/`, and `config/`.

### External data

Large sequencing and intermediate files are not version-controlled.

Local external inputs are accessed through:

```text
data/external/
```

On the original analysis workstation these paths are symbolic links to data stored outside the repository.

External inputs include:

```text
data/external/raw_reads/
data/external/trimmed_reads/
data/external/ophiostoma_bams/
data/external/fusarium_assembly/
data/external/fusarium_salmon/
data/external/fusarium_annotation/
data/external/fusarium_eggnog/
```

Examples of canonical external inputs include:

* paired-end sequencing reads
* trimmed sequencing reads
* coordinate-sorted *O. novo-ulmi* BAM files
* the final *Fusarium* Trinity assembly
* Salmon quantifications
* TransDecoder protein predictions
* eggNOG-mapper annotations
* dbCAN results

These files and local symbolic links are excluded from Git.

See:

```text
data/external/README.md
```

for details.

---

## Analysis workflow

### *Ophiostoma novo-ulmi*

The *O. novo-ulmi* analysis uses a reference-guided workflow.

```text
validated alignment inputs
        │
        ▼
gene-level counting
(featureCounts)
        │
        ▼
DESeq2 dataset construction and QC
        │
        ▼
differential expression
        │
        ├── functional annotation
        ├── GO enrichment
        ├── KEGG enrichment
        ├── secretome analysis
        ├── CAZyme analysis
        └── spatial transcriptional analysis
```

The repository includes validation of reference annotations, sample identities, gene-count matrices, contrast definitions, and downstream result consistency.

### *Fusarium salinense*

The *Fusarium* analysis uses a de novo transcriptome-based workflow.

```text
Trinity assembly + Salmon quantifications
        │
        ▼
gene-level abundance reconstruction
        │
        ▼
DESeq2
        │
        ▼
differential expression
        │
        ├── eggNOG functional annotation
        ├── GO enrichment
        ├── KEGG enrichment
        ├── SignalP-based secretome analysis
        └── dbCAN-based CAZyme analysis
```

The Trinity assembly, Salmon quantifications, TransDecoder proteins, and large annotation outputs are treated as external upstream inputs.

---

## Comparative and publication analyses

Cross-species manuscript analyses are stored primarily under:

```text
results/publication/
```

These include:

* global transcriptional response summaries
* PCA analyses
* differential-expression extent and direction
* COG annotation coverage and enrichment
* species × COG interaction tests
* comparative CAZyme analyses
* extracellular-response summaries
* *O. novo-ulmi* spatial response analyses
* audit and validation outputs

Publication figure scripts are located in:

```text
workflow/scripts/
```

Final publication figures and previews are stored under:

```text
figures/publication/
```

---

## Canonical manuscript result source

Manuscript-facing numerical results are consolidated in:

```text
results/publication/MANUSCRIPT_RESULTS_SOURCE.tsv
```

with provenance information in:

```text
results/publication/MANUSCRIPT_RESULTS_SOURCE_run_info.tsv
```

The table is generated by:

```text
workflow/scripts/40_build_manuscript_results_source.py
```

The generating workflow extracts values from validated analysis tables rather than manually entering manuscript numbers. Canonical input files are SHA-256 hashed, and expected biological mappings and result structures are checked before the manuscript source table is produced.

This table should be treated as the primary computational source for numerical values reported in the manuscript.

---

## Reproducibility and validation

The repository was organized to minimize silent inconsistencies between analyses.

Safeguards include:

* canonical sample metadata
* explicit biological contrast definitions
* repository-relative paths
* deterministic sample ordering
* input validation
* expected gene and transcript count checks
* contrast consistency checks
* annotation coverage checks
* cross-species comparability checks
* SHA-256 checksums
* machine-readable audit tables
* explicit failure when expected analytical constraints are violated

Supporting validation material is stored under locations including:

```text
data/metadata/
results/publication/audit/
workflow/rules/
```

---

## Software environment

A reconstructed Conda environment for running the released analysis workflow is provided in:

```text
environment.yml
```

Create it with:

```bash
conda env create -f environment.yml
conda activate dual-rnaseq-elm-fungi
```

An additional R dependency used for spatial heatmap generation can be installed with:

```bash
Rscript workflow/setup/install_r_extras.R
```

Core software versions verified on the analysis workstation included:

```text
R              4.4.3
samtools       1.19.2
featureCounts  2.0.6
```

The R analysis stack included DESeq2, apeglm, clusterProfiler, ComplexHeatmap, tximport, and associated CRAN/Bioconductor packages.

Python scripts use primarily the Python standard library, with:

* `pandas` for construction of the *O. novo-ulmi* KEGG annotation table
* PyMuPDF for assembly of composite publication figures

The repository evolved across more than one computational environment. The supplied `environment.yml` is therefore a curated environment for reproducing the released workflow rather than an exact export of every historical software environment used during development.

---

## Upstream software

Several upstream products used by this repository were generated with software that is not invoked directly by the current manuscript-analysis scripts.

These include tools such as:

* fastp
* HISAT2
* Trinity
* Salmon
* TransDecoder
* eggNOG-mapper
* SignalP
* dbCAN

Their processed outputs are supplied to the manuscript workflow as validated external inputs.

---

## Running the analysis

Clone the repository and create the software environment:

```bash
git clone <repository-url>
cd dual-rnaseq-elm-fungi

conda env create -f environment.yml
conda activate dual-rnaseq-elm-fungi

Rscript workflow/setup/install_r_extras.R
```

Configure the external data paths described in:

```text
data/external/README.md
```

Analysis scripts are located in:

```text
workflow/scripts/
```

Scripts should be executed from the repository root so that repository-relative paths resolve consistently.

The workflow is currently represented as a set of individually validated analysis scripts rather than a single workflow-manager entry point. Dependencies between major stages are reflected by their numbered filenames and validated input/output tables.

---

## Citation

Citation information for this repository is provided in:

```text
CITATION.cff
```

The manuscript citation should be used once the associated article is published.

---

## License

See `LICENSE` for the terms governing reuse of the repository code.

