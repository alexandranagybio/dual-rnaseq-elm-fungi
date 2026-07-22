# External data inputs

This directory provides repository-relative access to large sequencing files and
intermediate data that are not stored in Git.

The local workflow expects the following paths:

| Repository path | Contents |
|---|---|
| `data/external/raw_reads/` | Original paired-end FASTQ files |
| `data/external/trimmed_reads/` | fastp-trimmed paired-end FASTQ files |
| `data/external/ophiostoma_bams/` | Coordinate-sorted Ophiostoma BAM files |
| `data/external/fusarium_assembly/` | Final Fusarium Trinity assembly |
| `data/external/fusarium_salmon/` | Fusarium Salmon index and quantifications |

On the original analysis workstation, these paths are symbolic links to files
stored outside the repository.

## Example setup

```bash
ln -s /path/to/raw_reads data/external/raw_reads
ln -s /path/to/trimmed_reads data/external/trimmed_reads
ln -s /path/to/ophiostoma_alignments data/external/ophiostoma_bams
ln -s /path/to/trinity_final data/external/fusarium_assembly
ln -s /path/to/fusarium_salmon data/external/fusarium_salmon
```

Canonical publication inputs include:

- `data/external/ophiostoma_bams/*.bam`
- `data/external/fusarium_assembly/Fusarium_pure.Trinity.fasta`
- `data/external/fusarium_salmon/quants/*/quant.sf`

The sequencing data and symbolic links are excluded from Git. Only this README
is version-controlled.
