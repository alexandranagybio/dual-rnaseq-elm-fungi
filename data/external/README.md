# External data inputs

This directory provides repository-relative access to large sequencing files,
assemblies, quantifications, and annotation outputs that are not stored in Git.

The local workflow expects the following paths:

| Repository path | Contents |
|---|---|
| `data/external/raw_reads/` | Original paired-end FASTQ files |
| `data/external/trimmed_reads/` | fastp-trimmed paired-end FASTQ files |
| `data/external/ophiostoma_bams/` | Coordinate-sorted *Ophiostoma novo-ulmi* BAM files |
| `data/external/fusarium_assembly/` | Final *Fusarium salinense* Trinity assembly |
| `data/external/fusarium_salmon/` | *Fusarium* Salmon index and quantifications |
| `data/external/fusarium_annotation/` | *Fusarium* annotation inputs used by the released workflow |
| `data/external/fusarium_eggnog/` | Complete *Fusarium* eggNOG-mapper output directory |

On the original analysis workstation, these paths may be symbolic links to files
stored outside the repository.

## Example setup

```bash
ln -s /path/to/raw_reads data/external/raw_reads
ln -s /path/to/trimmed_reads data/external/trimmed_reads
ln -s /path/to/ophiostoma_alignments data/external/ophiostoma_bams
ln -s /path/to/trinity_final data/external/fusarium_assembly
ln -s /path/to/fusarium_salmon data/external/fusarium_salmon
ln -s /path/to/fusarium_annotation data/external/fusarium_annotation
ln -s /path/to/fusarium_eggnog data/external/fusarium_eggnog
```

Canonical external inputs include:

- `data/external/ophiostoma_bams/*.bam`
- `data/external/fusarium_assembly/Fusarium_pure.Trinity.fasta`
- `data/external/fusarium_salmon/quants/*/quant.sf`
- `data/external/fusarium_annotation/`
- `data/external/fusarium_eggnog/`

The sequencing data, large upstream products, and local symbolic links are
excluded from Git. Only this README is version-controlled.
