# dRNA Methylation Pipeline

A Snakemake pipeline for RNA base modification calling from Oxford Nanopore direct RNA sequencing (dRNA-seq) data. Based on the [DOGME](https://github.com/mortazavilab/dogme) workflow and the approach of Laurens Lambrechts, adapted for SGE cluster execution.

---
## Repository structure

dRNA_methylation_pipeline/
├── README.md               — this file
├── Snakefile               — main pipeline (edit MODS and modkit_pileup here)
├── cluster.yaml            — SGE resource requests per rule
├── submit.sh               — cluster submission script
├── cluster_qsub.sh         — SGE jobscript template
├── config/
│   └── config.yaml         — paths, containers, samples, filterbed thresholds
└── env/
    └── minimap.yaml        — conda environment (samtools, minimap2, bedtools)


## Workflow

```
pod5 files (per sample)
    └── 0. dorado_models_download  — download models once (skipped if present)
    └── 1. dorado_basecall         — basecall each pod5 → unaligned BAM
    └── 2. merge_unaligned_bam     — merge per-pod5 BAMs into one per sample
    └── 3. minimap2_align          — splice-aware alignment (preserves MM/ML tags)
            ├── 4. sample_probs    — QC: modification probability distributions
            ├── 5. modkit_pileup   — call modifications → raw BED
            └── 6. modkit_summary  — QC: global modification summary
                    └── 7. filterbed   — filter by coverage and mod percentage
                            └── 8. splitbed    — split BED by modification type
                                    └── 9. qc_report   — final QC report
```

### Supported modifications

| Modification | Description | SAM code | grep code |
|---|---|---|---|
| m5C | 5-Methylcytosine | C+m | `m` |
| m6A | N6-Methyladenosine | A+a | `a` |
| inosine | Inosine (A→I editing) | A+17596 | `17596` |
| pseU | Pseudouridine | T+17802 | `17802` |
| 2OmeC | 2′-O-methylcytidine | C+19228 | `19228` |
| 2OmeA | 2′-O-methyladenosine | A+69426 | `69426` |
| 2OmeG | 2′-O-methylguanosine | G+19229 | `19229` |
| 2OmeU | 2′-O-methyluridine | T+19227 | `19227` |

---

## Requirements

### Software

- [Snakemake](https://snakemake.readthedocs.io) ≥ 7.0
- [Singularity](https://sylabs.io/singularity/) ≥ 3.8
- [Conda](https://docs.conda.io) or [Mamba](https://github.com/mamba-org/mamba)
- SGE cluster with GPU nodes

### Singularity containers

```bash
# Dorado basecaller
singularity pull ontresearch-dorado-1.4.0.sif docker://ontresearch/dorado:1.4.0

# modkit modification calling
singularity pull ontresearch-modkit-0.6.1.sif docker://ontresearch/modkit:0.6.1
```

### Conda environments

```bash
conda env create -f env/minimap.yaml   # samtools, minimap2, bedtools
```

### Reference files

Download from Gencode:

```bash
# Reference genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz

# GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
```

Generate the junction BED file once from the GTF (required for splice-aware alignment):

```bash
paftools.js gff2bed gencode.v49.primary_assembly.annotation.gtf > gencode.v49.basic.annotation.junc.bed
```

`paftools.js` is included with minimap2.

---

## Installation

```bash
git clone https://github.com/your-username/dRNA_methylation_pipeline.git
cd dRNA_methylation_pipeline
```

---

## Configuration

### 1. Edit `config/config.yaml`

Set all file paths and thresholds:

```yaml
# Singularity bind mount — adjust to your cluster mount point
singularity_bind: "/your/data/mount:/your/data/mount"

# Container paths (dorado v.1.4.0; modkit v.0.6.1)
dorado_sif: "/path/to/ontresearch-dorado-1.4.0.sif"
modkit_sif: "/path/to/ontresearch-modkit-0.6.1.sif"

# Dorado model — full names required with --models-directory (the names of each model can be cheched on https://software-docs.nanoporetech.com/dorado/latest/models/list/)
dorado_model: "rna004_130bps_sup@v5.3.0,rna004_130bps_sup@v5.3.0_inosine_m6A_2OmeA@v1,..."

# Local directory for pre-downloaded dorado models
dorado_models_dir: "/path/to/dorado_models"

# Reference files
genome:  "/path/to/GRCh38.p14.genome.fa"
gtf:     "/path/to/gencode.v49.primary_assembly.annotation.gtf"
gtf_bed: "/path/to/gencode.v49.basic.annotation.junc.bed"

# Filterbed thresholds (can be changed and re-applied without re-running pileup)
min_coverage: 5   # minimum reads covering a site (col 5 of modkit BED)
mod_pct:      5   # minimum % modified reads per site (col 11 of modkit BED)

# Samples — one entry per sample, pointing to its pod5 directory
samples:
  sample_01: "/path/to/sample_01/pod5/"
  sample_02: "/path/to/sample_02/pod5/"
```

Pod5 files must be organised in **separate directories per sample**.

### 2. Edit `Snakefile` — modkit pileup arguments

Open the Snakefile and find `rule modkit_pileup` (marked `<-- EDIT`). Adjust the modkit arguments directly in the shell command:

```python
shell:
    """
    modkit pileup \
        -t {threads} \
        --filter-threshold 0.7 \ # one line per base
        --mod-thresholds m:0.9 \
        {input.bam} {output.bed} \
        --log-filepath {log}
    """
```

**`--filter-threshold`** — global minimum base modification probability. Reads below this are excluded from the pileup entirely

**`--mod-thresholds BASE:FLOAT`** — per-canonical-base minimum probability for a modification call to be counted.



---

## Running the pipeline


## Output structure

```
results/
├── basecalled/
│   └── {sample}/
│       └── {sample}.merged.unaligned.bam     ← merged unaligned BAM
├── bams/
│   └── {sample}/
│       ├── {sample}.bam                       ← aligned BAM
│       └── {sample}.bam.bai
├── modkit/
│   ├── {sample}.raw.bed                       ← all sites, unfiltered
│   └── {sample}.summary.tsv                  ← modkit summary QC
├── bedMethyl/
│   ├── {sample}.filtered.bed                 ← all sites, filtered
│   ├── {sample}.m5C.filtered.bed
│   ├── {sample}.m6A.filtered.bed
│   ├── {sample}.inosine.filtered.bed
│   ├── {sample}.pseU.filtered.bed
│   ├── {sample}.2OmeC.filtered.bed
│   ├── {sample}.2OmeA.filtered.bed
│   ├── {sample}.2OmeG.filtered.bed
│   └── {sample}.2OmeU.filtered.bed
├── sample_probs/
│   └── {sample}/
│       ├── probabilities.tsv
│       ├── thresholds.tsv
│       ├── counts.html                        ← open in browser
│       └── proportion.html                    ← open in browser
└── qc/
    └── {sample}.qc_summary.txt               ← final QC report
```

If a modification is absent in the data, the corresponding BED file will be empty — the pipeline does not fail.

---

