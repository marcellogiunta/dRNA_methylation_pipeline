import os
import re
from pathlib import Path

# ================================================================
#
#   dRNA METHYLATION PIPELINE — Snakemake / SGE cluster
#   Based on DOGME (nanoporeModule.nf v1.2.3)
#   + Laurens Lambrechts approach (dorado v1.4.0 + junc.bed)
#
#   WORKFLOW:
#   0.  dorado_models_download — download dorado models once
#   1.  dorado_basecall        — basecall each pod5 → unaligned BAM
#   2.  merge_unaligned_bam    — merge BAMs per sample
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
MODKIT_SIF = config["modkit_sif"]
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
        expand("results/bedMethyl/{sample}.m5C.filtered.bed",      sample=SAMPLES),
        expand("results/bedMethyl/{sample}.m6A.filtered.bed",      sample=SAMPLES),
        expand("results/bedMethyl/{sample}.inosine.filtered.bed",  sample=SAMPLES),
        expand("results/bedMethyl/{sample}.pseU.filtered.bed",     sample=SAMPLES),
        expand("results/bedMethyl/{sample}.2OmeC.filtered.bed",    sample=SAMPLES),
        expand("results/bedMethyl/{sample}.2OmeA.filtered.bed",    sample=SAMPLES),
        expand("results/bedMethyl/{sample}.2OmeG.filtered.bed",    sample=SAMPLES),
        expand("results/bedMethyl/{sample}.2OmeU.filtered.bed",    sample=SAMPLES),
        expand("results/sample_probs/{sample}/probabilities.tsv",  sample=SAMPLES),
        expand("results/qc/{sample}.qc_summary.txt",               sample=SAMPLES)


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
        pod5=lambda wc: JOB_POD5[(wc.sample, wc.basename)],
        # ensures models are downloaded before basecalling starts
        models=expand(
            config["dorado_models_dir"] + "/rna004_130bps_sup@v5.3.0_{mod}",
            mod=["m5C_2OmeC@v1",
                 "inosine_m6A_2OmeA@v1",
                 "pseU_2OmeU@v1",
                 "2OmeG@v1"]
        )

    output:
        bam=temp("results/basecalled/{sample}/{basename}.unaligned.bam")

    log:
        "logs/dorado/{sample}_{basename}.log"

    params:
        model     =config["dorado_model"],
        models_dir=config["dorado_models_dir"]

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
            --device cuda:0 \
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

    threads: 8

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
        bam="results/bams/{sample}/{sample}.bam",
        bai="results/bams/{sample}/{sample}.bam.bai"

    log:
        "logs/minimap2/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 4

    shell:
        """
        mkdir -p results/bams/{wildcards.sample} logs/minimap2

        samtools bam2fq --threads {threads} -T MM,ML,pt \
            {input.bam} 2>> {log} | \
        minimap2 \
            -ax splice \
            -uf \
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
        bam="results/bams/{sample}/{sample}.bam",
        bai="results/bams/{sample}/{sample}.bam.bai"

    output:
        probs      ="results/sample_probs/{sample}/probabilities.tsv",
        thresholds ="results/sample_probs/{sample}/thresholds.tsv",
        counts_html="results/sample_probs/{sample}/counts.html",
        prop_html  ="results/sample_probs/{sample}/proportion.html"

    log:
        "logs/sample_probs/{sample}.log"

    params:
        outdir="results/sample_probs/{sample}"

    threads: 4

    shell:
        """
        mkdir -p {params.outdir} logs/sample_probs

        singularity exec -B {BIND} {MODKIT_SIF} \
        modkit sample-probs \
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
        bam="results/bams/{sample}/{sample}.bam",
        bai="results/bams/{sample}/{sample}.bam.bai"

    output:
        bed="results/modkit/{sample}.raw.bed"

    log:
        "logs/modkit/{sample}.log"

    threads: 8

    shell:
        """
        mkdir -p results/modkit logs/modkit

        singularity exec -B {BIND} {MODKIT_SIF} \
        modkit pileup \
            -t {threads} \
            --filter-threshold 0.7 \
            --mod-thresholds m:0.99 \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """


# =========================================================
# 6. MODKIT SUMMARY (QC)
# =========================================================
rule modkit_summary:
    input:
        bam="results/bams/{sample}/{sample}.bam",
        bai="results/bams/{sample}/{sample}.bam.bai"

    output:
        tsv="results/modkit/{sample}.summary.tsv"

    log:
        "logs/modkit/{sample}.summary.log"

    threads: 4

    shell:
        """
        mkdir -p logs/modkit

        singularity exec -B {BIND} {MODKIT_SIF} \
        modkit summary \
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
        bed="results/modkit/{sample}.raw.bed"

    output:
        bed="results/bedMethyl/{sample}.filtered.bed"

    log:
        "logs/filterbed/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    params:
        min_coverage=config["min_coverage"],
        mod_pct     =config["mod_pct"]

    shell:
        """
        mkdir -p results/bedMethyl logs/filterbed

        awk 'NR==1 || $1~/^#/ || ($5 >= {params.min_coverage} && $11 >= {params.mod_pct})' \
            {input.bed} > {output.bed} 2> {log}

        echo "Raw:      $(grep -vc '^#' {input.bed}  || echo 0)" >> {log}
        echo "Filtered: $(grep -vc '^#' {output.bed} || echo 0)" >> {log}
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
        bed="results/bedMethyl/{sample}.filtered.bed"

    output:
        m5C   ="results/bedMethyl/{sample}.m5C.filtered.bed",
        m6A   ="results/bedMethyl/{sample}.m6A.filtered.bed",
        inosine="results/bedMethyl/{sample}.inosine.filtered.bed",
        pseU  ="results/bedMethyl/{sample}.pseU.filtered.bed",
        ome2C ="results/bedMethyl/{sample}.2OmeC.filtered.bed",
        ome2A ="results/bedMethyl/{sample}.2OmeA.filtered.bed",
        ome2G ="results/bedMethyl/{sample}.2OmeG.filtered.bed",
        ome2U ="results/bedMethyl/{sample}.2OmeU.filtered.bed"

    log:
        "logs/splitbed/{sample}.log"

    wildcard_constraints:
        sample="|".join(re.escape(s) for s in SAMPLES)

    shell:
        """
        mkdir -p logs/splitbed

        grep -w 'm'     {input.bed} > {output.m5C}    2>> {log} || touch {output.m5C}
        grep -w 'a'     {input.bed} > {output.m6A}    2>> {log} || touch {output.m6A}
        grep -w '17596' {input.bed} > {output.inosine} 2>> {log} || touch {output.inosine}
        grep -w '17802' {input.bed} > {output.pseU}   2>> {log} || touch {output.pseU}
        grep -w '19228' {input.bed} > {output.ome2C}  2>> {log} || touch {output.ome2C}
        grep -w '69426' {input.bed} > {output.ome2A}  2>> {log} || touch {output.ome2A}
        grep -w '19229' {input.bed} > {output.ome2G}  2>> {log} || touch {output.ome2G}
        grep -w '19227' {input.bed} > {output.ome2U}  2>> {log} || touch {output.ome2U}

        echo "m5C:     $(wc -l < {output.m5C})"    >> {log}
        echo "m6A:     $(wc -l < {output.m6A})"    >> {log}
        echo "inosine: $(wc -l < {output.inosine})" >> {log}
        echo "pseU:    $(wc -l < {output.pseU})"   >> {log}
        echo "2OmeC:   $(wc -l < {output.ome2C})"  >> {log}
        echo "2OmeA:   $(wc -l < {output.ome2A})"  >> {log}
        echo "2OmeG:   $(wc -l < {output.ome2G})"  >> {log}
        echo "2OmeU:   $(wc -l < {output.ome2U})"  >> {log}
        """


# =========================================================
# 9. FINAL QC REPORT
# =========================================================
rule qc_report:
    input:
        bam_unaligned="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        bam          ="results/bams/{sample}/{sample}.bam",
        bai          ="results/bams/{sample}/{sample}.bam.bai",
        summary      ="results/modkit/{sample}.summary.tsv",
        sample_probs ="results/sample_probs/{sample}/probabilities.tsv",
        raw          ="results/modkit/{sample}.raw.bed",
        filtered     ="results/bedMethyl/{sample}.filtered.bed",
        m5C          ="results/bedMethyl/{sample}.m5C.filtered.bed",
        m6A          ="results/bedMethyl/{sample}.m6A.filtered.bed",
        inosine      ="results/bedMethyl/{sample}.inosine.filtered.bed",
        pseU         ="results/bedMethyl/{sample}.pseU.filtered.bed",
        ome2C        ="results/bedMethyl/{sample}.2OmeC.filtered.bed",
        ome2A        ="results/bedMethyl/{sample}.2OmeA.filtered.bed",
        ome2G        ="results/bedMethyl/{sample}.2OmeG.filtered.bed",
        ome2U        ="results/bedMethyl/{sample}.2OmeU.filtered.bed"

    output:
        report="results/qc/{sample}.qc_summary.txt"

    log:
        "logs/qc/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 2

    shell:
        """
        mkdir -p results/qc

        {{
        echo "===== QC SUMMARY — {wildcards.sample} ====="
        echo "Date: $(date)"
        echo ""

        echo "--- Total reads (unaligned BAM) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Alignment (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- Reads per chromosome (idxstats, top 25) ---"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Modification probabilities (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Modified sites ---"
        echo "Raw:      $(grep -vc '^#' {input.raw}      || true)"
        echo "Filtered: $(grep -vc '^#' {input.filtered} || true)"
        echo ""

        echo "--- Sites per modification type ---"
        echo "m5C:     $(wc -l < {input.m5C})"
        echo "m6A:     $(wc -l < {input.m6A})"
        echo "inosine: $(wc -l < {input.inosine})"
        echo "pseU:    $(wc -l < {input.pseU})"
        echo "2OmeC:   $(wc -l < {input.ome2C})"
        echo "2OmeA:   $(wc -l < {input.ome2A})"
        echo "2OmeG:   $(wc -l < {input.ome2G})"
        echo "2OmeU:   $(wc -l < {input.ome2U})"
        echo ""

        echo "--- modkit summary ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """
