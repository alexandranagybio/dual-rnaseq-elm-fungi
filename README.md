# Dual RNA-seq workflow for transcriptomic analysis of fungal antagonism between *Ophiostoma novo-ulmi* and *Fusarium cf. salinense*

## Overview

This repository contains the complete computational workflow used to analyze dual RNA-seq data from interactions between the Dutch elm disease pathogen *Ophiostoma novo-ulmi* and the fungal endophyte *Fusarium cf. salinense*. The workflow accompanies a transcriptomic study investigating molecular responses associated with fungal antagonism.

Because the two organisms required different analytical strategies, the repository contains two complementary pipelines. *Ophiostoma novo-ulmi* gene expression is quantified using reference-guided alignment and gene-level counting, whereas *Fusarium cf. salinense* gene expression is quantified from a de novo assembled transcriptome. Downstream analyses include differential expression, functional annotation, Gene Ontology (GO) enrichment, and generation of publication-ready tables and figures.

The repository has been systematically audited to maximize computational reproducibility. Analysis scripts perform extensive validation of inputs and intermediate results, generate checksum manifests, and produce machine-readable validation reports to ensure transparent and repeatable analyses.

---

## Experimental design

Dual RNA-seq libraries were generated from dual cultures of *Ophiostoma novo-ulmi* and *Fusarium cf. salinense*, together with species-specific control cultures. Biological triplicates were analyzed for each experimental condition.

Following sequencing, reads originating from *Ophiostoma novo-ulmi* were analyzed using a reference-guided workflow based on genome alignment and gene-level counting. Reads not assigned to *Ophiostoma* were assembled de novo to reconstruct the *Fusarium cf. salinense* transcriptome, which was subsequently used for transcript quantification and downstream differential expression analysis.

---

## Repository structure

```text
.
├── config/              Canonical sample metadata
├── data/                Reference data, annotations and supporting resources
├── figures/             Publication figures
├── logs/                Workflow execution logs
├── manuscript/          Manuscript source files
├── results/             Generated analysis results
├── workflow/            Analysis scripts, workflow rules and environments
├── CHANGELOG.md         Repository development history
├── CITATION.cff         Citation metadata
├── environment.yml      Conda software environment
├── LICENSE              Repository license
└── README.md
```

---

## Workflow overview

The repository implements two complementary transcriptomic analysis pipelines.

```text
                           Raw RNA-seq reads
                                   │
                    ┌──────────────┴──────────────┐
                    │                             │
                    ▼                             ▼
          Ophiostoma novo-ulmi        Fusarium cf. salinense
                    │                             │
          Reference genome            De novo transcriptome
               alignment                   assembly
               (HISAT2)                   (Trinity)
                    │                             │
           Gene-level counting         Transcript quantification
           (featureCounts)                 (Salmon)
                    │                             │
                    ▼                             ▼
                  DESeq2                       DESeq2
                    │                             │
        Functional annotation       Functional annotation
    (eggNOG-mapper, SignalP, dbCAN for both pipelines)
                    └──────────────┬──────────────┘
                                   ▼
                       Gene Ontology enrichment
                         (clusterProfiler)
                                   │
                                   ▼
                           Tables and figures
```

Although the analytical strategies differ between the two organisms, both pipelines follow the same principles of reproducible data processing, validation, and downstream functional interpretation.

---

## Repository organization

### `config/`

Contains the canonical sample metadata used throughout the repository.

The file

```text
config/samples.tsv
```

is the single authoritative description of the experimental design. Workflow scripts derive analysis-specific metadata directly from this file while validating consistency with the audited study design.

---

### `workflow/`

Contains the complete computational workflow, including:

* shell scripts
* Python utilities
* R analysis scripts
* workflow rules
* software environments

All scripts use repository-relative paths to maximize portability and reproducibility.

---

### `data/`

Contains reference genomes, transcriptomes, annotations, metadata, and supplementary resources required by the workflow.

---

### `results/`

Stores all generated analysis outputs.

```text
results/
├── ophiostoma/
│   ├── gene counts
│   ├── DESeq2 analyses
│   ├── functional annotation
│   └── GO enrichment
│
└── fusarium/
    ├── transcript quantification
    ├── DESeq2 analyses
    ├── functional annotation
    └── GO enrichment
```

Results are generated automatically by the workflow and should not be edited manually.

---

### `figures/`

Contains publication figures generated from the workflow outputs.

---

### `manuscript/`

Contains manuscript source files associated with this repository.

---

## Reproducibility

Computational reproducibility is a central design objective of this repository.

The workflow incorporates multiple safeguards, including:

* repository-relative file paths
* canonical sample metadata
* deterministic sample ordering
* strict input validation
* automatic integrity checks
* SHA-256 checksum manifests
* machine-readable validation reports
* explicit failure when expected biological or technical constraints are violated

These validation steps prevent silent propagation of errors and ensure that downstream analyses are performed only on verified inputs.

---

## Software requirements

The primary computational environment is defined in:

```text
environment.yml
```

The workflow uses widely adopted open-source software, including:

### Read processing and quantification

* fastp
* HISAT2
* samtools
* featureCounts (Subread)
* Trinity
* Salmon

### Statistical analysis

* R
* DESeq2

### Functional annotation

* eggNOG-mapper
* SignalP
* dbCAN

### Functional enrichment

* clusterProfiler

---

## Running the workflow

Clone the repository:

```bash
git clone <repository-url>
cd dual-rnaseq-elm-fungi
```

Workflow components are executed from the repository root. Individual analysis scripts are located in:

```text
workflow/scripts/
```

Intermediate and final outputs are written to the `results/` directory, where validated outputs from one stage serve as inputs for subsequent analyses. Validation reports and checksum manifests are generated automatically throughout the workflow to document computational integrity.

---

## Main outputs

The repository produces publication-ready outputs, including:

* validated count matrices and transcript quantifications
* DESeq2 datasets
* differential expression tables
* functional annotation tables
* Gene Ontology enrichment analyses
* quality-control summaries
* validation reports
* SHA-256 checksum manifests
* publication figures

---

## Citation

If you use this repository in your research, please cite the associated publication when available.

Citation metadata are provided in:

```text
CITATION.cff
```

---

## License

This repository is distributed under the terms of the accompanying `LICENSE` file.

