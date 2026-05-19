import os
from pathlib import Path

# ================================================================
#
#   WORKFLOW:
#   1.  dorado_basecall    — basecall each pod5 file separately
#   2.  merge_unaligned    — merging all the BAM files per sample
#   3.  minimap2_align     — alignement
#   4.  sample_probs       — modkit sample-probs to inspect base modification probabilities (QC)
#   5.  separate_strands   — separate strands +/- (RNA strand-specific)
#   6.  modkit_pileup      — calling modification
#   7.  modkit_summary     — modkit summary (QC)
#   8.  filterbed          — filtering for coverage and minimum percentage of modification
#   9.  splitbed           — split BED files (according to the modification)
#   10. qc_report          — final report 
#
# ================================================================

configfile: "config/config.yaml"

DORADO_SIF  = config["dorado_sif"]
MODKIT_SIF  = config["modkit_sif"]


# =========================================================
# Files detecting (the pipeline need the pod files organized in separated folders per sample)
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
        expand("results/bedMethyl/{sample}.plus.m5C.filtered.bed",      sample=SAMPLES),
        expand("results/bedMethyl/{sample}.minus.m5C.filtered.bed",     sample=SAMPLES),
        expand("results/bedMethyl/{sample}.plus.m6A.filtered.bed",      sample=SAMPLES),
        expand("results/bedMethyl/{sample}.minus.m6A.filtered.bed",     sample=SAMPLES),
        expand("results/bedMethyl/{sample}.plus.pseU.filtered.bed",     sample=SAMPLES),
        expand("results/bedMethyl/{sample}.minus.pseU.filtered.bed",    sample=SAMPLES),
        expand("results/bedMethyl/{sample}.plus.inosine.filtered.bed",  sample=SAMPLES),
        expand("results/bedMethyl/{sample}.minus.inosine.filtered.bed", sample=SAMPLES),
        # Sample probs (QC modificazioni)
        expand("results/sample_probs/{sample}/probabilities.tsv",       sample=SAMPLES),
        # QC report finale
        expand("results/qc/{sample}.qc_summary.txt",                    sample=SAMPLES)


# =========================================================
# 1. DORADO BASECALL (1 POD5 file→ 1 not aligned BAM )
#
# - NO --reference: minimap2 alignment in the next step
# - --emit-moves: required to keep modification tags

# Use a conteiner ans specify the path in config.yam
# =========================================================
rule dorado_basecall:
    input:
        pod5=lambda wc: JOB_POD5[(wc.sample, wc.basename)]

    output:
        bam=temp("results/basecalled/{sample}/{basename}.unaligned.bam")

    log:
        "logs/dorado/{sample}_{basename}.log"

    params:
        model=config["dorado_model"]

    threads: 2

    resources:
        gpu=1

    shell:
        """
        mkdir -p results/basecalled/{wildcards.sample} logs/dorado

        singularity exec -B /SAN/vyplab:/SAN/vyplab --nv {DORADO_SIF} \
        dorado basecaller {params.model} {input.pod5} \
            --emit-moves \
            --estimate-poly-a \
            --device cuda:auto \
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
# Requires junc.bed files — path to be provided in config["gtf_bed"].
#
# Notes about arguments:
# - bam2fq -T MM,ML,pt → exports modifications tags  from the unaligned BAM
# - minimap2 -y         → transfers the tags to the output BAM
# - -G 500000           → max Intron lenght 
# - -L                  → long CIGAR (required by modkit)
# - --secondary=no      → only primary alignment
# - --MD                → tag for MD  mismatch (for QC)
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
# 4. MODKIT SAMPLE-PROBS (QC for modifications probabilities)
#
# Works on the merged-aligned BAM.
# It plots the probability of modifications according to the thresholds values-
#
# Outputs per sample:
#   - probabilities.tsv  — probabilities distribution per mod
#   - thresholds.tsv     — estimate threshold per mod
#   - counts.html        — histogram of counts (--hist)
#   - proportion.html    — histogram of proportions (--hist)
#
# Use a container and specify the path in (config["modkit_sif"])
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

        singularity exec -B /SAN/vyplab:/SAN/vyplab {MODKIT_SIF} \
        modkit sample-probs \
            --hist \
            --threads {threads} \
            --out-dir {params.outdir} \
            {input.bam} \
        2> {log}
        """


# =========================================================
# 5. SEPARAZIONE STRAND +/- SEPARATION 
#
# In RNA direct ONT mods are strand-specific.
# samtools view:
#   -F 16 = strand +
#   -f 16 = strand -
# =========================================================
rule separate_strands:
    input:
        bam="results/bams/{sample}/{sample}.bam",
        bai="results/bams/{sample}/{sample}.bam.bai"

    output:
        plus_bam ="results/bams/{sample}/{sample}.plus.bam",
        plus_bai ="results/bams/{sample}/{sample}.plus.bam.bai",
        minus_bam="results/bams/{sample}/{sample}.minus.bam",
        minus_bai="results/bams/{sample}/{sample}.minus.bam.bai"

    log:
        "logs/strands/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 8

    shell:
        """
        mkdir -p logs/strands

        samtools view -b -F 16 --threads {threads} \
            {input.bam} -o {output.plus_bam} 2>> {log}
        samtools index -@ {threads} {output.plus_bam} 2>> {log}

        samtools view -b -f 16 --threads {threads} \
            {input.bam} -o {output.minus_bam} 2>> {log}
        samtools index -@ {threads} {output.minus_bam} 2>> {log}
        """


# =========================================================
# 6. MODKIT PILEUP per strand
#
# In the example the defaults thresholds are used (automatically calculated by modkit)
# For specific values use --filter-threshold {params.filter_threshold} (passing the requested values via the config file).
# Same for --mod-threshold (e.g. --mod-threshold C: 0.9, for m5c specific threshold)
# =========================================================
rule modkit_pileup:
    input:
        bam="results/bams/{sample}/{sample}.{strand}.bam",
        bai="results/bams/{sample}/{sample}.{strand}.bam.bai"

    output:
        bed="results/modkit/{sample}.{strand}.raw.bed"

    log:
        "logs/modkit/{sample}.{strand}.log"

    params:
        outdir          ="results/modkit",
        filter_threshold=config["filter_threshold"],
        mod_threshold   =config["mod_threshold"]

    wildcard_constraints:
        strand="plus|minus"

    threads: 8

    shell:
        """
        mkdir -p {params.outdir} logs/modkit

        singularity exec -B /SAN/vyplab:/SAN/vyplab {MODKIT_SIF} \
        modkit pileup \
            -t {threads} \
            {input.bam} {output.bed} \
            --log-filepath {log}
        """
## removed      --filter-threshold {params.filter_threshold} \


# =========================================================
# 7. MODKIT SUMMARY (global QC , on the merged BAM )
#
# --tsv         → output 
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

        singularity exec -B /SAN/vyplab:/SAN/vyplab {MODKIT_SIF} \
        modkit summary \
            --threads {threads} \
            --tsv \
            {input.bam} > {output.tsv} 2> {log}
        """


# =========================================================
# 8. FILTERBED
#
# Filtering on the columns (BED modkit pileup output):
#   col 5  ($5)  = valid_coverage
#   col 11 ($11) = fraction_modified (0-100)
# =========================================================
rule filterbed:
    input:
        bed="results/modkit/{sample}.{strand}.raw.bed"

    output:
        bed="results/bedMethyl/{sample}.{strand}.filtered.bed"

    log:
        "logs/filterbed/{sample}.{strand}.log"

    wildcard_constraints:
        strand="plus|minus"

    params:
        min_coverage=config["min_coverage"],
        mod_pct=config["mod_pct"]

    shell:
        """
        mkdir -p results/bedMethyl logs/filterbed

        awk 'NR==1 || $1~/^#/ || ($5 >= {params.min_coverage} && $11 >= {params.mod_pct})' \
            {input.bed} > {output.bed} 2> {log}

        echo "Raw (post modkit):           $(grep -vc '^#' {input.bed}   || echo 0)" >> {log}
        echo "Filtrate (cov>={params.min_coverage}, mod>={params.mod_pct}%): $(grep -vc '^#' {output.bed} || echo 0)" >> {log}
        """


# =========================================================
# 9. SPLITBED 
#
# Codici modkit per RNA:
#   'm'                        = m5C
#   'a'                        = m6A / inosine_m6A (stesso codice)
#   '17596'                    = inosine (A→I editing)
#   '17802'                    = pseudouridine (pseU)
#   '19228|19229|19227|69426'  = Nm (2'-O-methylation)
# =========================================================
rule splitbed:
    input:
        bed="results/bedMethyl/{sample}.{strand}.filtered.bed"

    output:
        m5C    ="results/bedMethyl/{sample}.{strand}.m5C.filtered.bed",
        m6A    ="results/bedMethyl/{sample}.{strand}.m6A.filtered.bed",
        inosine="results/bedMethyl/{sample}.{strand}.inosine.filtered.bed",
        pseU   ="results/bedMethyl/{sample}.{strand}.pseU.filtered.bed",
        Nm     ="results/bedMethyl/{sample}.{strand}.Nm.filtered.bed"

    log:
        "logs/splitbed/{sample}.{strand}.log"

    wildcard_constraints:
        strand="plus|minus"

    shell:
        """
        mkdir -p logs/splitbed

        grep -w 'm'                        {input.bed} > {output.m5C}     2>> {log} || touch {output.m5C}
        grep -w 'a'                        {input.bed} > {output.m6A}     2>> {log} || touch {output.m6A}
        grep -w '17596'                    {input.bed} > {output.inosine} 2>> {log} || touch {output.inosine}
        grep -w '17802'                    {input.bed} > {output.pseU}    2>> {log} || touch {output.pseU}
        grep -Ew '19228|19229|19227|69426' {input.bed} > {output.Nm}      2>> {log} || touch {output.Nm}

        echo "m5C:     $(wc -l < {output.m5C})"     >> {log}
        echo "m6A:     $(wc -l < {output.m6A})"     >> {log}
        echo "inosine: $(wc -l < {output.inosine})" >> {log}
        echo "pseU:    $(wc -l < {output.pseU})"    >> {log}
        echo "Nm:      $(wc -l < {output.Nm})"      >> {log}
        """


# =========================================================
# 10. QC FINAL REPORT 
# =========================================================
rule qc_report:
    input:
        bam_unaligned ="results/basecalled/{sample}/{sample}.merged.unaligned.bam",
        bam           ="results/bams/{sample}/{sample}.bam",
        bai           ="results/bams/{sample}/{sample}.bam.bai",
        summary       ="results/modkit/{sample}.summary.tsv",
        sample_probs  ="results/sample_probs/{sample}/probabilities.tsv",
        plus_raw      ="results/modkit/{sample}.plus.raw.bed",
        minus_raw     ="results/modkit/{sample}.minus.raw.bed",
        plus_filtered ="results/bedMethyl/{sample}.plus.filtered.bed",
        minus_filtered="results/bedMethyl/{sample}.minus.filtered.bed",
        m5C_plus      ="results/bedMethyl/{sample}.plus.m5C.filtered.bed",
        m5C_minus     ="results/bedMethyl/{sample}.minus.m5C.filtered.bed",
        m6A_plus      ="results/bedMethyl/{sample}.plus.m6A.filtered.bed",
        m6A_minus     ="results/bedMethyl/{sample}.minus.m6A.filtered.bed"

    output:
        report="results/qc/{sample}.qc_summary.txt"

    log:
        "logs/qc/{sample}.log"

    conda:
        "env/minimap.yaml"

    threads: 2

    shell:
        """
        set -x
        mkdir -p results/qc

        {{
        echo "===== QC SUMMARY — {wildcards.sample} ====="
        echo "Data: $(date)"
        echo ""

        echo "--- Total reads (Unaligned BAM) ---"
        echo "Reads: $(samtools view -c -@ {threads} {input.bam_unaligned})"
        echo ""

        echo "--- Alignment  (flagstat) ---"
        samtools flagstat --threads {threads} {input.bam}
        echo ""

        echo "--- Reads per chromosome (idxstats) ---"
        echo "chr | length | mapped | unmapped"
        samtools idxstats {input.bam} | sort -k3 -rn | head -25 || true
        echo ""

        echo "--- Modification prob (sample-probs) ---"
        cat {input.sample_probs}
        echo ""

        echo "--- Modified sites STRAND + ---"
        echo "Raw:      $(grep -vc '^#' {input.plus_raw}      || true)"
        echo "Filtrate: $(grep -vc '^#' {input.plus_filtered} || true)"
        echo "m5C:      $(wc -l < {input.m5C_plus})"
        echo "m6A:      $(wc -l < {input.m6A_plus})"
        echo ""

        echo "--- PModified sites STRAND - ---"
        echo "Raw:      $(grep -vc '^#' {input.minus_raw}      || true)"
        echo "Filtrate: $(grep -vc '^#' {input.minus_filtered} || true)"
        echo "m5C:      $(wc -l < {input.m5C_minus})"
        echo "m6A:      $(wc -l < {input.m6A_minus})"
        echo ""

        echo "--- Summary modifications (modkit summary) ---"
        cat {input.summary}
        }} > {output.report} 2> {log}
        """    
