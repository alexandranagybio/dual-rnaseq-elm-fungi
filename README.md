# Ophiostoma gene-level recount scripts

These files are designed for the split layout discovered on 18 July 2026:

- Audit repository: `~/dual-rnaseq-elm-fungi`
- Original RNA-seq data: `~/rnaseq`

Place the three scripts in:

```text
~/dual-rnaseq-elm-fungi/workflow/scripts/
```

## Important prerequisite

`04_count_ophiostoma_genes.sh` must use the **validated GTF produced yesterday**,
not an arbitrary GFF3. The GTF must contain:

- exon features;
- a `gene_id` attribute on every counted exon;
- exactly 8,640 unique genes represented by the exon features.

The script tries to find a unique plausible GTF. If several candidates exist,
it stops rather than guessing. In that case, supply the exact file:

```bash
OPN_GTF=/absolute/path/to/the/validated_ophiostoma.gtf \
bash workflow/scripts/04_count_ophiostoma_genes.sh
```

## Installation

From a terminal:

```bash
cd ~/dual-rnaseq-elm-fungi

mkdir -p workflow/scripts

cp ~/Downloads/04_count_ophiostoma_genes.sh workflow/scripts/
cp ~/Downloads/05_prepare_deseq2_matrix.py workflow/scripts/
cp ~/Downloads/06_build_ophiostoma_deseq2_dataset.R workflow/scripts/

chmod +x workflow/scripts/04_count_ophiostoma_genes.sh
chmod +x workflow/scripts/05_prepare_deseq2_matrix.py
chmod +x workflow/scripts/06_build_ophiostoma_deseq2_dataset.R
```

Depending on browser download settings, the files may be somewhere other than
`~/Downloads`; use their actual downloaded locations.

## Run sequence

### 1. Count genes

```bash
cd ~/dual-rnaseq-elm-fungi
bash workflow/scripts/04_count_ophiostoma_genes.sh
```

This performs paired-end, forward-stranded counting:

```text
featureCounts -p --countReadPairs -B -C -s 1 -t exon -g gene_id
```

It stops unless exactly 8,640 gene rows are produced.

### 2. Prepare the clean DESeq2 inputs

```bash
python3 workflow/scripts/05_prepare_deseq2_matrix.py
```

This creates:

```text
results/ophiostoma/gene_counts/ophiostoma_gene_counts_matrix.tsv
results/ophiostoma/gene_counts/sample_metadata.tsv
results/ophiostoma/gene_counts/05_prepare_deseq2_matrix.validation.json
results/ophiostoma/gene_counts/05_prepare_deseq2_matrix.sha256
```

Canonical sample design:

| Sample | Condition | Replicate |
|---|---|---|
| 149 | interaction | 1 |
| 150 | interaction | 2 |
| 151 | interaction | 3 |
| 152 | self | 1 |
| 153 | self | 2 |
| 154 | self | 3 |
| 155 | onu | 1 |
| 156 | onu | 2 |
| 157 | onu | 3 |

### 3. Build the DESeq2 object

Activate the R environment containing DESeq2, then run:

```bash
Rscript workflow/scripts/06_build_ophiostoma_deseq2_dataset.R
```

This creates an **unfiltered, unanalysed** DESeq2 object with:

```r
design = ~ condition
```

and `self` as the reference level. It does not yet run `DESeq()` or produce
differential-expression results.

## Outputs from the counting step

```text
results/ophiostoma/gene_counts/ophiostoma_gene_featurecounts.txt
results/ophiostoma/gene_counts/ophiostoma_gene_featurecounts.txt.summary
results/ophiostoma/gene_counts/04_count_ophiostoma_genes.run_info.txt
results/ophiostoma/gene_counts/04_count_ophiostoma_genes.sha256
logs/04_count_ophiostoma_genes.log
```

## Why the scripts stop aggressively

This audit is intended to be publication-grade. The scripts therefore refuse to
continue when:

- the annotation is ambiguous;
- BAMs or BAM indices are missing;
- `samtools quickcheck` fails;
- the gene universe is not 8,640;
- sample names do not match 149-157;
- genes are duplicated;
- counts are negative, missing, or non-integer;
- count columns and metadata rows are misaligned.
