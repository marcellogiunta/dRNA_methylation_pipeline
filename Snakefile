import os
import re
from pathlib import Path

# ================================================================
#
#   dRNA METHYLATION PIPELINE — Snakemake / SGE cluster
#   Based on DOGME (nanoporeModule.nf v1.2.3)
#
#   WORKFLOW:
#   0.  dorado_models_download — download dorado models once
#   1.  dorado_basecall        — basecall each pod5 → unaligned BAM
#   2.  merge_unaligned_bam    — merge BAMs per sample
#   FROM THIS STEP ONWARD THE PIPELINE PROCESSES IN PARALLEL GENOME AND TRANSCRIPTOME ALIGNED BAM FILES
#   3.  minimap2_align         — splice-aware alignment
#   4.  sample_probs           — modkit sample-probs QC
#   5.  modkit_pileup          — call modifications
#   6.  modkit_summary         — modkit summary QC
#   7.  filterbed              — filter by coverage and mod percentage
#   8.  splitbed               — split BED by modification type
#   9.  qc_report              — final QC report
#
# ================================================================
#
#   USER CONFIGURATION — edit sections marked with  <-- EDIT
#
#   1. Paths and thresholds → config/config.yaml
#   2. modkit pileup arguments → rule modkit_pileup below
#
# ================================================================

 configfile: "config/config.yaml"

DORADO_SIF = config["dorado_sif"]
MODKIT = config["modkit_bin"]
BIND       = config["singularity_bind"]


# =========================================================
# FILE DISCOVERY
# =========================================================
def collect_jobs():
    jobs = []
    for sample, path in config["samples"].items():
        p = Path(path).expanduser().resolve()
        for f in sorted(p.glob("*.pod5")):
            jobs.append({
                "sample"  : sample,
                "pod5"    : str(f),
                "basename": f.stem
            })
    return jobs

JOBS    = collect_jobs()
SAMPLES = list(config["samples"].keys())

JOB_POD5 = {
    (j["sample"], j["basename"]): j["pod5"]
    for j in JOBS
}


# =========================================================
# FINAL OUTPUT
# =========================================================
rule all:
    input:
        # BED filtrati e separati per tipo di modifica
        expand("results/bedMethyl_genome/{sample}/{sample}.m5C.filtered.bed.gz",      sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.m6A.filtered.bed.gz",      sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.pseU.filtered.bed.gz",     sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.inosine.filtered.bed.gz",  sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.2OmeC.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.2OmeA.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.2OmeG.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_genome/{sample}/{sample}.2OmeU.filtered.bed.gz",       sample=SAMPLES),
        # Sample probs (QC modificazioni)
        expand("results/sample_probs_genome/{sample}/probabilities.tsv",  sample=SAMPLES),
        # QC report finale
        expand("results/qc_genome/{sample}.qc_summary.txt",               sample=SAMPLES),
        # BED filtrati e separati per tipo di modifica — trascrittoma
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.m5C.filtered.bed.gz",      sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.m6A.filtered.bed.gz",      sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.pseU.filtered.bed.gz",     sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.inosine.filtered.bed.gz",  sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.2OmeC.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.2OmeA.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.2OmeG.filtered.bed.gz",       sample=SAMPLES),
        expand("results/bedMethyl_transcriptome/{sample}/{sample}.2OmeU.filtered.bed.gz",       sample=SAMPLES),
        # Sample probs (QC modificazioni) — trascrittoma
        expand("results/sample_probs_transcriptome/{sample}/probabilities.tsv",  sample=SAMPLES),
        # QC report finale — trascrittoma
        expand("results/qc_transcriptome/{sample}.qc_summary.txt",               sample=SAMPLES)


# =========================================================
# 0. DORADO MODELS DOWNLOAD
#
# Downloads all required models to dorado_models_dir.
# Skipped automatically by Snakemake if model directories
# already exist on disk.
# =========================================================
rule dorado_models_download:
    output:
        directory(config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0"),
        directory(config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0_m5C_2OmeC@v1"),
        directory(config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0_inosine_m6A_2OmeA@v1"),
        directory(config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0_pseU_2OmeU@v1"),
        directory(config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0_2OmeG@v1"),

    params:
        model_dir=config["dorado_models_dir"]

    log:
        "logs/dorado_download/download.log"

    shell:
        """
        mkdir -p {params.model_dir} logs/dorado_download

        for model in \
            rna004_130bps_sup@v5.3.0 \
            rna004_130bps_sup@v5.3.0_m5C_2OmeC@v1 \
            rna004_130bps_sup@v5.3.0_inosine_m6A_2OmeA@v1 \
            rna004_130bps_sup@v5.3.0_pseU_2OmeU@v1 \
            rna004_130bps_sup@v5.3.0_2OmeG@v1; do

            echo "Downloading $model..." >> {log}
            singularity exec -B {BIND} {DORADO_SIF} \
                dorado download \
                    --model $model \
                    --directory {params.model_dir} 2>> {log}
        done

        echo "Done:" >> {log}
        ls {params.model_dir} >> {log}
        """


# =========================================================
# 1. DORADO BASECALL (1 POD5 → 1 unaligned BAM)
#
# --emit-moves:        required to preserve modification tags
# --estimate-poly-a:   polyA tail estimation
# --batchsize 64:      conservative, safe across GPU types
# --device cuda:0:     use GPU assigned by SGE
# --models-directory:  use pre-downloaded local models
# =========================================================
rule dorado_basecall:
    input:
        pod5=lambda wc: JOB_POD5[(wc.sample, wc.basename)]

    output:
        bam=temp("results/basecalled/{sample}/{basename}.unaligned.bam")

    log:
        "logs/dorado/{sample}_{basename}.log"

    params:
        model=config["dorado_model"],
        models_dir =config["dorado_models_dir"]

    threads: 2

    resources:
        gpu=1

    shell:
        """
        mkdir -p results/basecalled/{wildcards.sample} logs/dorado

        singularity exec -B {BIND} --nv {DORADO_SIF} \
        dorado basecaller {params.model} {input.pod5} \
            --emit-moves \
            --estimate-poly-a \
            --batchsize 64 \
            --device cuda:auto \
            --models-directory {params.models_dir} \
        > {output.bam} 2>> {log}
        """


# =========================================================
# 2. MERGE UNALIGNED BAM PER SAMPLE
# =========================================================
rule merge_unaligned_bam:
    input:
        lambda wc: [
            f"results/basecalled/{j['sample']}/{j['basename']}.unaligned.bam"
            for j in JOBS if j["sample"] == wc.sample
        ]

    output:
        bam="results/basecalled/{sample}/{sample}.merged.unaligned.bam"

    log:
        "logs/merge_unaligned/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 1

    shell:
        """
        mkdir -p logs/merge_unaligned
        samtools merge -f --threads {threads} {output.bam} {input} 2> {log}
        """


# =========================================================
# 3. MINIMAP2 SPLICE-AWARE ALIGNMENT
#
# junc.bed is pre-existing — path in config["gtf_bed"].
#
# -y:             transfers MM/ML/pt tags from unaligned BAM
#                 WITHOUT -y modification tags are lost
# -G 500000:      max intron length 500kb
# -L:             long CIGAR format (required by modkit)
# --secondary=no: primary alignments only
# --MD:           MD tag for mismatches (QC)
# =========================================================
rule minimap2_align:
    input:
        bam    ="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        genome =config["genome"],
        juncbed=config["gtf_bed"]

    output:
        bam="results/bams_genome/{sample}/{sample}.bam",
        bai="results/bams_genome/{sample}/{sample}.bam.bai"

    log:
        "logs/minimap2_genome/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 4

    shell:
        """
        mkdir -p results/bams_genome/{wildcards.sample} logs/minimap2_genome

        samtools bam2fq --threads {threads} -T MM,ML,pt \
            {input.bam} 2>> {log} | \
        minimap2 \
            -ax splice \
            -uf \
            -k14 \
            -G 500000 \
            -L \
            --secondary=no \
            --MD \
            -y \
            --junc-bed {input.juncbed} \
            -t {threads} \
            {input.genome} - 2>> {log} | \
        samtools sort --threads {threads} -o {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """


# =========================================================
# 4. MODKIT SAMPLE-PROBS (QC)
#
# Checks MM/ML tags are present and probability distributions
# look sensible before running the full pileup.
# Inspect counts.html and proportion.html in a browser.
# =========================================================
rule sample_probs:
    input:
        bam="results/bams_genome/{sample}/{sample}.bam",
        bai="results/bams_genome/{sample}/{sample}.bam.bai"

    output:
        probs      ="results/sample_probs_genome/{sample}/probabilities.tsv",
        thresholds ="results/sample_probs_genome/{sample}/thresholds.tsv",
        counts_html="results/sample_probs_genome/{sample}/counts.html",
        prop_html  ="results/sample_probs_genome/{sample}/proportion.html"

    log:
        "logs/sample_probs_genome/{sample}.log"

    params:
        outdir="results/sample_probs_genome/{sample}"

    threads: 1

    shell:
        """
        mkdir -p {params.outdir} logs/sample_probs_genome

        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} sample-probs \
            --hist \
            --threads {threads} \
            --out-dir {params.outdir} \
            {input.bam} \
        2> {log}
        """


# =========================================================
# 5. MODKIT PILEUP  <-- EDIT arguments here
#
# Runs on the merged aligned BAM (no strand separation).
# Strand info is preserved in col 6 of the output BED.
#
# Edit the modkit pileup arguments below to suit your analysis:
#
# --filter-threshold FLOAT
#     Global minimum base modification probability threshold.
#     Reads below this are excluded from pileup entirely.
#     Default used here: 0.7 If you want to use different threshold for different bases, use a different argument per row. Example
#     --filter-threshold C:0.8
#     --filter-threshold A:0.9
#
# --mod-thresholds mod:FLOAT
#     Per-canonical-base minimum probability for a modification
#     call to be counted. Format: <mod>:<threshold>
#     Canonical bases for RNA modifications:
#       C → m5C (m), 2OmeC (19228)
#       A → m6A (a), inosine (17596), 2OmeA (69426)
#       T → pseU (17802), 2OmeU (19227)
#       G → 2OmeG (19229)
#     Example: m:0.9,a:0.7,T:0.7,G:0.7
#
# NOTE: --min-coverage is NOT a valid modkit pileup argument.
# Coverage filtering is applied downstream in filterbed (col $5).
# The raw BED is always preserved for re-filtering.
# =========================================================
rule modkit_pileup:
    input:
        bam="results/bams_genome/{sample}/{sample}.bam",
        bai="results/bams_genome/{sample}/{sample}.bam.bai"

    output:
        bed="results/modkit_genome/{sample}.raw.bed.gz"

    log:
        "logs/modkit_genome/{sample}.log"

    params:
        outdir="results/modkit_genome",
        genome =config["genome"]

    threads: 1

    shell:
        """
        mkdir -p {params.outdir} logs/modkit_genome

        
        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} pileup \
            --mod-threshold m:0.99 \
            --modified-bases C:m \
            --reference {params.genome} \
            --bgzf \
            -t {threads} \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """

# =========================================================
# 6. MODKIT SUMMARY (QC)
# =========================================================
rule modkit_summary:
    input:
        bam="results/bams_genome/{sample}/{sample}.bam",
        bai="results/bams_genome/{sample}/{sample}.bam.bai"

    output:
        tsv="results/modkit_genome/{sample}.summary.tsv"

    log:
        "logs/modkit_genome/{sample}.summary.log"

    threads: 1

    shell:
        """
        mkdir -p logs/modkit_genome

        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} summary \
            --threads {threads} \
            --tsv \
            {input.bam} > {output.tsv} 2> {log}
        """

# =========================================================
# 7. FILTERBED  <-- thresholds in config/config.yaml
#
# modkit pileup BED columns:
#   col 5  ($5)  = valid_coverage
#   col 11 ($11) = fraction_modified (0-100)
#
# Raw BED is always preserved. To re-filter with different
# thresholds without re-running pileup:
#   snakemake --forcerun filterbed
# =========================================================
rule filterbed:
    input:
        bed="results/modkit_genome/{sample}.raw.bed.gz"

    output:
        bed="results/bedMethyl_genome/{sample}.filtered.bed.gz"

    log:
        "logs/filterbed_genome/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    params:
        min_coverage=config["min_coverage"],
        mod_pct     =config["mod_pct"]

    shell:
        """
        mkdir -p results/bedMethyl_genome logs/filterbed_genome
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | \
        awk 'NR==1 || $1~/^#/ || ($5 >= {params.min_coverage} && $11 >= {params.mod_pct})' | \
        singularity exec -B {BIND} {DORADO_SIF} bgzip -c \
        > {output.bed} 2> {log}

        echo "Raw:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | grep -vc '^#' || echo 0)" >> {log}
        echo "Filtered: $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.bed} | grep -vc '^#' || echo 0)" >> {log}
        """


# =========================================================
# 8. SPLITBED
#
# Splits filtered BED by modification type using modkit
# SAM codes (col 4 of pileup BED output).
#
# All modifications are extracted regardless of which models
# were used during basecalling — if a modification is absent
# in the BED the output file will be empty (no pipeline error).
#
# SAM codes reference (from modkit RNA modifications table):
#   m5C     → C+m      → grep code: m
#   m6A     → A+a      → grep code: a
#   inosine → A+17596  → grep code: 17596
#   pseU    → T+17802  → grep code: 17802
#   2OmeC   → C+19228  → grep code: 19228
#   2OmeA   → A+69426  → grep code: 69426
#   2OmeG   → G+19229  → grep code: 19229
#   2OmeU   → T+19227  → grep code: 19227
# =========================================================
rule splitbed:
    input:
        bed="results/bedMethyl_genome/{sample}.filtered.bed.gz"

    output:
        m5C    ="results/bedMethyl_genome/{sample}/{sample}.m5C.filtered.bed.gz",
        m6A    ="results/bedMethyl_genome/{sample}/{sample}.m6A.filtered.bed.gz",
        inosine="results/bedMethyl_genome/{sample}/{sample}.inosine.filtered.bed.gz",
        pseU   ="results/bedMethyl_genome/{sample}/{sample}.pseU.filtered.bed.gz",
        OmeC   ="results/bedMethyl_genome/{sample}/{sample}.2OmeC.filtered.bed.gz",
        OmeA   ="results/bedMethyl_genome/{sample}/{sample}.2OmeA.filtered.bed.gz",
        OmeG   ="results/bedMethyl_genome/{sample}/{sample}.2OmeG.filtered.bed.gz",
        OmeU   ="results/bedMethyl_genome/{sample}/{sample}.2OmeU.filtered.bed.gz"

    log:
        "logs/splitbed_genome/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    shell:
        """
        mkdir -p results/bedMethyl_genome/{wildcards.sample} logs/splitbed_genome

        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "m"'     | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.m5C}     || touch {output.m5C}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "a"'     | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.m6A}     || touch {output.m6A}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "17596"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.inosine} || touch {output.inosine}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "17802"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.pseU}    || touch {output.pseU}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19228"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeC}    || touch {output.OmeC}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "69426"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeA}    || touch {output.OmeA}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19229"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeG}    || touch {output.OmeG}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19227"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeU}    || touch {output.OmeU}

        echo "m5C:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.m5C}     | wc -l)" >> {log}
        echo "m6A:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.m6A}     | wc -l)" >> {log}
        echo "inosine: $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.inosine} | wc -l)" >> {log}
        echo "pseU:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.pseU}    | wc -l)" >> {log}
        echo "2OmeC:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeC}    | wc -l)" >> {log}
        echo "2OmeA:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeA}    | wc -l)" >> {log}
        echo "2OmeG:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeG}    | wc -l)" >> {log}
        echo "2OmeU:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeU}    | wc -l)" >> {log}
        """


# =========================================================
# 9. FINAL QC REPORT
# =========================================================
rule qc_report:
    input:
        bam_unaligned="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        bam          ="results/bams_genome/{sample}/{sample}.bam",
        bai          ="results/bams_genome/{sample}/{sample}.bam.bai",
        summary      ="results/modkit_genome/{sample}.summary.tsv",
        sample_probs ="results/sample_probs_genome/{sample}/probabilities.tsv",
        raw          ="results/modkit_genome/{sample}.raw.bed.gz",
        filtered     ="results/bedMethyl_genome/{sample}.filtered.bed.gz",
        m5C          ="results/bedMethyl_genome/{sample}/{sample}.m5C.filtered.bed.gz",
        m6A          ="results/bedMethyl_genome/{sample}/{sample}.m6A.filtered.bed.gz",
        inosine      ="results/bedMethyl_genome/{sample}/{sample}.inosine.filtered.bed.gz",
        pseU         ="results/bedMethyl_genome/{sample}/{sample}.pseU.filtered.bed.gz",
        OmeC        ="results/bedMethyl_genome/{sample}/{sample}.2OmeC.filtered.bed.gz",
        OmeA        ="results/bedMethyl_genome/{sample}/{sample}.2OmeA.filtered.bed.gz",
        OmeG        ="results/bedMethyl_genome/{sample}/{sample}.2OmeG.filtered.bed.gz",
        OmeU        ="results/bedMethyl_genome/{sample}/{sample}.2OmeU.filtered.bed.gz"

    output:
        report="results/qc_genome/{sample}.qc_summary.txt"

    log:
        "logs/qc_genome/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 1

    shell:
        """
        mkdir -p results/qc_genome

        {{
        echo "===== QC SUMMARY — {wildcards.sample} ====="
        echo "Date: $(date)"
        echo ""

        echo "--- Total reads (BAM non allineato) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Genomic alignment  (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- Reads per  chromosome (idxstats) ---"
        echo "chr | length | mapped | unmapped"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Modifications probabilities (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Posizioni modificate ---"
        echo "Raw:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.raw}       | grep -vc '^#'  || true)"
        echo "Filtrate: $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.filtered}  | grep -vc '^#'  || true)"
        echo "m5C:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.m5C}       | wc -l)"
        echo "m6A:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.m6A}       | wc -l)"
        echo "inosine:  $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.inosine}   | wc -l)"
        echo "pseU:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.pseU}      | wc -l)"
        echo "2OmeC:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeC}      | wc -l)"
        echo "2OmeA:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeA}      | wc -l)"  
        echo "2OmeG:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeG}      | wc -l)"
        echo "2OmeU:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeU}      | wc -l)"                      
        echo ""

        echo "--- Summary modifications (modkit summary) ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """


# =========================================================
# =========================================================
# TRASCRITTOME ALIGNED BRANCH
# Identycal to the genomuc branch but with optimizations for transcriptomic alignment.
# =========================================================
# =========================================================

# =========================================================
# 3t. MINIMAP2 ALLINEAMENTO — TRASCRITTOMA
#
# Differences from genomic:
# - target: FASTA trasnctipts (config["transcriptome"])
# - -ax map-ont:    optimized fot ONT
# - No -uf:         not relevant for transcriptome
# - No --junc-bed:  not relevant for transcriptome
# - No -G:          no introns
# =========================================================
rule minimap2_align_transcriptome:
    input:
        bam          ="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        transcriptome=config["transcriptome"]

    output:
        bam="results/bams_transcriptome/{sample}/{sample}.bam",
        bai="results/bams_transcriptome/{sample}/{sample}.bam.bai"

    log:
        "logs/minimap2_transcriptome/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 4

    shell:
        """
        mkdir -p results/bams_transcriptome/{wildcards.sample} logs/minimap2_transcriptome

        samtools bam2fq --threads {threads} -T MM,ML,pt \
            {input.bam} 2>> {log} | \
        minimap2 \
            -ax map-ont \
            -k14 \
            -L \
            --secondary=no \
            --MD \
            -y \
            -t {threads} \
            {input.transcriptome} - 2>> {log} | \
        samtools sort --threads {threads} -o {output.bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}
        """


# =========================================================
# 4t. MODKIT SAMPLE-PROBS — TRANSCRIPTOME
# =========================================================
rule sample_probs_transcriptome:
    input:
        bam="results/bams_transcriptome/{sample}/{sample}.bam",
        bai="results/bams_transcriptome/{sample}/{sample}.bam.bai"

    output:
        probs      ="results/sample_probs_transcriptome/{sample}/probabilities.tsv",
        thresholds ="results/sample_probs_transcriptome/{sample}/thresholds.tsv",
        counts_html="results/sample_probs_transcriptome/{sample}/counts.html",
        prop_html  ="results/sample_probs_transcriptome/{sample}/proportion.html"

    log:
        "logs/sample_probs_transcriptome/{sample}.log"

    params:
        outdir="results/sample_probs_transcriptome/{sample}"

    threads: 1

    shell:
        """
        mkdir -p {params.outdir} logs/sample_probs_transcriptome

        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} sample-probs \
            --hist \
            --threads {threads} \
            --out-dir {params.outdir} \
            {input.bam} \
        2> {log}
        """


# =========================================================
# 5t. MODKIT PILEUP — TRANSCRIPTOME
# =========================================================
rule modkit_pileup_transcriptome:
    input:
        bam="results/bams_transcriptome/{sample}/{sample}.bam",
        bai="results/bams_transcriptome/{sample}/{sample}.bam.bai"

    output:
        bed="results/modkit_transcriptome/{sample}.raw.bed.gz"

    log:
        "logs/modkit_transcriptome/{sample}.log"

    params:
        outdir="results/modkit_transcriptome",
        transcriptome =config["transcriptome"]

    threads: 1

    shell:
        """
        mkdir -p {params.outdir} logs/modkit_transcriptome

        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} pileup \
            --mod-threshold m:0.99 \
            --modified-bases C:m \
            --ref {params.transcriptome} \
            --preload-references \
            --bgzf \
            -t {threads} \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """


# =========================================================
# 6t. MODKIT SUMMARY — TRANSCRIPTOME
# =========================================================
rule modkit_summary_transcriptome:
    input:
        bam="results/bams_transcriptome/{sample}/{sample}.bam",
        bai="results/bams_transcriptome/{sample}/{sample}.bam.bai"

    output:
        tsv="results/modkit_transcriptome/{sample}.summary.tsv"

    log:
        "logs/modkit_transcriptome/{sample}.summary.log"

    threads: 1

    shell:
        """
        mkdir -p logs/modkit_transcriptome

        singularity exec -B {BIND} \
        {DORADO_SIF} \
        {MODKIT} summary \
            --threads {threads} \
            --tsv \
            {input.bam} > {output.tsv} 2> {log}
        """


# =========================================================
# 7t. FILTERBED — TRANSCRIPTOME
# =========================================================
rule filterbed_transcriptome:
    input:
        bed="results/modkit_transcriptome/{sample}.raw.bed.gz"

    output:
        bed="results/bedMethyl_transcriptome/{sample}.filtered.bed.gz"

    log:
        "logs/filterbed_transcriptome/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    params:
        min_coverage=config["min_coverage"],
        mod_pct     =config["mod_pct"]

    shell:
        """
        mkdir -p results/bedMethyl_transcriptome logs/filterbed_transcriptome

        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | \
        awk 'NR==1 || $1~/^#/ || ($5 >= {params.min_coverage} && $11 >= {params.mod_pct})' | \
        singularity exec -B {BIND} {DORADO_SIF} bgzip -c \
        > {output.bed} 2> {log}

        echo "Raw:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | grep -vc '^#' || echo 0)" >> {log}
        echo "Filtered: $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.bed} | grep -vc '^#' || echo 0)" >> {log}
        """


# =========================================================
# 8t. SPLITBED — TRANSCRIPTOME
# =========================================================
rule splitbed_transcriptome:
    input:
        bed="results/bedMethyl_transcriptome/{sample}.filtered.bed.gz"

    output:
        m5C    ="results/bedMethyl_transcriptome/{sample}/{sample}.m5C.filtered.bed.gz",
        m6A    ="results/bedMethyl_transcriptome/{sample}/{sample}.m6A.filtered.bed.gz",
        inosine="results/bedMethyl_transcriptome/{sample}/{sample}.inosine.filtered.bed.gz",
        pseU   ="results/bedMethyl_transcriptome/{sample}/{sample}.pseU.filtered.bed.gz",
        OmeC   ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeC.filtered.bed.gz",
        OmeA   ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeA.filtered.bed.gz",
        OmeG   ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeG.filtered.bed.gz",
        OmeU   ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeU.filtered.bed.gz"

    log:
        "logs/splitbed_transcriptome/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    shell:
        """
        mkdir -p results/bedMethyl_transcriptome/{wildcards.sample} logs/splitbed_transcriptome

        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "m"'     | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.m5C}     || touch {output.m5C}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "a"'     | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.m6A}     || touch {output.m6A}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "17596"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.inosine} || touch {output.inosine}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "17802"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.pseU}    || touch {output.pseU}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19228"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeC}    || touch {output.OmeC}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "69426"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeA}    || touch {output.OmeA}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19229"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeG}    || touch {output.OmeG}
        singularity exec -B {BIND} {DORADO_SIF} zcat {input.bed} | awk '$4 == "19227"' | singularity exec -B {BIND} {DORADO_SIF} bgzip -c > {output.OmeU}    || touch {output.OmeU}

        echo "m5C:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.m5C}     | wc -l)" >> {log}
        echo "m6A:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.m6A}     | wc -l)" >> {log}
        echo "inosine: $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.inosine} | wc -l)" >> {log}
        echo "pseU:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.pseU}    | wc -l)" >> {log}
        echo "2OmeC:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeC}    | wc -l)" >> {log}
        echo "2OmeA:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeA}    | wc -l)" >> {log}
        echo "2OmeG:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeG}    | wc -l)" >> {log}
        echo "2OmeU:   $(singularity exec -B {BIND} {DORADO_SIF} zcat {output.OmeU}    | wc -l)" >> {log}
        """


# =========================================================
# 9t. QC REPORT — TRANSCRIPTOME
# =========================================================
rule qc_report_transcriptome:
    input:
        bam_unaligned="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        bam          ="results/bams_transcriptome/{sample}/{sample}.bam",
        bai          ="results/bams_transcriptome/{sample}/{sample}.bam.bai",
        summary      ="results/modkit_transcriptome/{sample}.summary.tsv",
        sample_probs ="results/sample_probs_transcriptome/{sample}/probabilities.tsv",
        raw          ="results/modkit_transcriptome/{sample}.raw.bed.gz",
        filtered     ="results/bedMethyl_transcriptome/{sample}.filtered.bed.gz",
        m5C          ="results/bedMethyl_transcriptome/{sample}/{sample}.m5C.filtered.bed.gz",
        m6A          ="results/bedMethyl_transcriptome/{sample}/{sample}.m6A.filtered.bed.gz",
        inosine      ="results/bedMethyl_transcriptome/{sample}/{sample}.inosine.filtered.bed.gz",
        pseU         ="results/bedMethyl_transcriptome/{sample}/{sample}.pseU.filtered.bed.gz",
        OmeC        ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeC.filtered.bed.gz",
        OmeA        ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeA.filtered.bed.gz",
        OmeG        ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeG.filtered.bed.gz",
        OmeU        ="results/bedMethyl_transcriptome/{sample}/{sample}.2OmeU.filtered.bed.gz"

    output:
        report="results/qc_transcriptome/{sample}.qc_summary.txt"

    log:
        "logs/qc_transcriptome/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 1

    shell:
        """
        mkdir -p results/qc_transcriptome

        {{
        echo "===== QC SUMMARY TRASCRITTOMA — {wildcards.sample} ====="
        echo "Data: $(date)"
        echo ""

        echo "--- Total reads (BAM non allineato) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Transcriptomic alignment (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- Reads  per transcript (idxstats) ---"
        echo "transcript | length | mapped | unmapped"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Mods probabilities (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Posizioni modificate ---"
        echo "Raw:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.raw}       | grep -vc '^#'  || true)"
        echo "Filtrate: $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.filtered}  | grep -vc '^#'  || true)"
        echo "m5C:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.m5C}       | wc -l)"
        echo "m6A:      $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.m6A}       | wc -l)"
        echo "inosine:  $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.inosine}   | wc -l)"
        echo "pseU:     $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.pseU}      | wc -l)"
        echo "2OmeC:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeC}      | wc -l)"
        echo "2OmeA:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeA}      | wc -l)"  
        echo "2OmeG:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeG}      | wc -l)"
        echo "2OmeU:    $(singularity exec -B {BIND} {DORADO_SIF} zcat {input.OmeU}      | wc -l)"                      
        echo ""

        echo "--- Summary modifications (modkit summary) ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """

