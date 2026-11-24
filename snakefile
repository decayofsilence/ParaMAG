"""
ParaMAG: assembly and binning of intracellular parasite MAGs
from single-cell (MDA) sequencing data.

Input:  samples.tsv - sample name and paths to R1/R2 reads
        config.yaml - parameters
Output: results/<sample>/...

Run from the repository root:
    snakemake --cores 32

"""
import csv

configfile: "config.yaml"


# ---------- Samples ----------
SAMPLES = {}
with open(config["samples"]) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        SAMPLES[row["sample"]] = row


def reads_r1(wildcards):
    """Path to the R1 reads of a sample."""
    return SAMPLES[wildcards.sample]["r1"]


def reads_r2(wildcards):
    """Path to the R2 reads of a sample."""
    return SAMPLES[wildcards.sample]["r2"]


# ---------- Bins ----------
BIN_SAMPLES, BIN_IDS = glob_wildcards("results/{sample}/maxbin/bin.{binid}.fasta")


def bin_ids(sample):
    """Bin ids found on disk for one sample, e.g. ["001", "002"]."""
    ids = []
    for s, b in zip(BIN_SAMPLES, BIN_IDS):
        if s == sample:
            ids.append(b)
    return ids


def bin_fastas(wildcards):
    """All bin FASTA files of a sample."""
    return expand("results/{sample}/maxbin/bin.{binid}.fasta",
                  sample=wildcards.sample, binid=bin_ids(wildcards.sample))


def busco_summaries(wildcards):
    """BUSCO summary files for all bins of a sample."""
    return expand("results/{sample}/busco/bin.{binid}/short_summary.specific.fungi_odb10.txt",
                  sample=wildcards.sample, binid=bin_ids(wildcards.sample))


# ---------- Final targets ----------
rule all:
    input:
        expand("results/{sample}/qc/fastqc_done", sample=SAMPLES),
        expand("results/{sample}/assembly/contigs.fasta", sample=SAMPLES),
        expand("results/{sample}/mapping/abundance_data.tsv", sample=SAMPLES),
        # the report needs BUSCO of every bin and the tree, so they are built too
        expand("results/{sample}/report/summary_report.md", sample=SAMPLES),


# ---------- Steps ----------
rule fastqc:
    input:
        r1=reads_r1,
        r2=reads_r2,
    output:
        "results/{sample}/qc/fastqc_done"
    shell:
        """
        fastqc {input.r1} {input.r2} -o results/{wildcards.sample}/qc/
        touch {output}
        """


rule trimmomatic:
    input:
        r1=reads_r1,
        r2=reads_r2,
    output:
        r1_paired="results/{sample}/trimmed/R1_paired.fastq",
        r1_unpaired="results/{sample}/trimmed/R1_unpaired.fastq",
        r2_paired="results/{sample}/trimmed/R2_paired.fastq",
        r2_unpaired="results/{sample}/trimmed/R2_unpaired.fastq",
    params:
        steps=config["trimmomatic_steps"],
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            {params.steps}
        """


rule spades:
    input:
        r1="results/{sample}/trimmed/R1_paired.fastq",
        r2="results/{sample}/trimmed/R2_paired.fastq",
    output:
        "results/{sample}/assembly/contigs.fasta"
    params:
        mode=config["spades_mode"],
    shell:
        "spades.py {params.mode} -1 {input.r1} -2 {input.r2} -o results/{wildcards.sample}/assembly/"


rule abundance:
    input:
        contigs="results/{sample}/assembly/contigs.fasta",
        r1="results/{sample}/trimmed/R1_paired.fastq",
        r2="results/{sample}/trimmed/R2_paired.fastq",
    output:
        bam="results/{sample}/mapping/aln.bam",
        depth="results/{sample}/mapping/abundance_data.tsv",
    params:
        index="results/{sample}/mapping/contigs_index",
    shell:
        """
        bowtie2-build {input.contigs} {params.index}
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output.bam}
        samtools index {output.bam}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}
        """


rule coverage:
    input:
        "results/{sample}/mapping/aln.bam"
    output:
        "results/{sample}/mapping/coverage.tsv"
    shell:
        "samtools idxstats {input} > {output}"


rule maxbin:
    input:
        contigs="results/{sample}/assembly/contigs.fasta",
        abund="results/{sample}/mapping/abundance_data.tsv",
    output:
        "results/{sample}/maxbin/bin.001.fasta",
        "results/{sample}/maxbin/bin.002.fasta",
        "results/{sample}/maxbin/bin.003.fasta",
        "results/{sample}/maxbin/bin.004.fasta",
        "results/{sample}/maxbin/bin.005.fasta",
    shell:
        "run_MaxBin.pl -contig {input.contigs} -out results/{wildcards.sample}/maxbin/bin -abund {input.abund}"


rule busco:
    input:
        "results/{sample}/maxbin/bin.{binid}.fasta"
    output:
        "results/{sample}/busco/bin.{binid}/short_summary.specific.fungi_odb10.txt"
    shell:
        "busco -i {input} -l fungi_odb10 -o results/{wildcards.sample}/busco/bin.{wildcards.binid} -m genome"


rule align_bins:
    input:
        bin_fastas
    output:
        "results/{sample}/phylogeny/aligned_bins.fasta"
    shell:
        """
        cat {input} > results/{wildcards.sample}/phylogeny/concatenated_bins.fasta
        mafft --auto results/{wildcards.sample}/phylogeny/concatenated_bins.fasta > {output}
        """


rule phylogenomics:
    input:
        "results/{sample}/phylogeny/aligned_bins.fasta"
    output:
        "results/{sample}/phylogeny/super_tree_output.tre"
    shell:
        """
        iqtree2 -s {input} -m MFP -bb 1000 -nt AUTO -pre results/{wildcards.sample}/phylogeny/iqtree_output
        cp results/{wildcards.sample}/phylogeny/iqtree_output.treefile {output}
        """


rule report:
    input:
        summaries=busco_summaries,
        tree="results/{sample}/phylogeny/super_tree_output.tre",
    output:
        "results/{sample}/report/summary_report.md"
    shell:
        """
        echo "# Pipeline Summary Report" > {output}
        echo "\n## BUSCO Summaries" >> {output}
        for f in {input.summaries}; do echo "\n### $$f" >> {output}; tail -n 10 $$f >> {output}; done
        echo "\n## Phylogenetic Tree File" >> {output}
        echo "Tree saved as: {input.tree}" >> {output}
        """
