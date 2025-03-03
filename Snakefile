## Snakefile for processing Micro-C data

## conda activate microc_env

configfile: "config/config.yml"

import re
import glob
import os

# Auto-detect samples from raw data directory
raw_fastq_files = glob.glob("raw-data/*.fastq*")
print("Raw fastq files:", raw_fastq_files)

reads = [re.sub(r"\.fastq(\.gz)?$", "", os.path.basename(f)) for f in raw_fastq_files]
print("Detected reads:", reads)

samples = list(set(re.sub(r"_[Rr][12].*", "", os.path.basename(f)) for f in raw_fastq_files))
print("Detected samples:", samples)

#####################################
#### RULES
#####################################

rule all:
    input:
        expand("results/fastqc/{read}.html", read=reads),
        expand("results/{sample}.mapped.PT.bam", sample=samples),
        expand("results/{sample}.duplication_stats_summary.txt", sample=samples),
        expand("results/{sample}.complexity.txt", sample=samples),
        expand("results/{sample}.contact_map.hic", sample=samples),



rule fastqc:
    input:
        "raw-data/{read}.fastq.gz"
    output:
        html="results/fastqc/{read}.html",
        zip="results/fastqc/{read}_fastqc.zip"
    log:
        "results/logs/fastqc_{read}.log"
    params:
        "--threads 4"
    wrapper:
        "v1.5.0/bio/fastqc"





# Single-Step Pipeline: FASTQ → Aligned BAM → Filtered Pairs → Final BAM
rule fastq_to_valid_pairs_bam:
    input:
        ref=config["reference_genome"],
        r1="raw-data/{sample}_R1.fastq.gz",
        r2="raw-data/{sample}_R2.fastq.gz",
        genome_file=config["genome_file"],
    output:
        bam="results/{sample}.mapped.PT.bam",
        pairs="results/{sample}.mapped.pairs",
        stats="results/{sample}.duplication_stats.txt"
    conda:
        "envs/microc_env.yaml"
    log:
        "results/logs/valid_pairs_{sample}.log"
    threads: config["threads"]
    shell:
        """
        bwa mem -5SP -T0 -t {threads} {input.ref} {input.r1} {input.r2} | \
        pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {input.genome_file} | \
        pairtools sort --nproc {threads} | \
        pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {output.stats} | \
        pairtools split --nproc-in {threads} --nproc-out {threads} --output-pairs {output.pairs} --output-sam - | \
        samtools view -bS -@ {threads} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index the final BAM file
        samtools index {output.bam}
        """

# Summarize stats
rule summarize_stats:
    input:
        stats="results/{sample}.duplication_stats.txt"
    output:
        summary="results/{sample}.duplication_stats_summary.txt"
    shell:
        """
        python ./scripts/get_qc.py -p {input.stats} > {output.summary}
        """
        
# Analyze Library Complexity - in this example the output file out.preseq will detail the extrapolated complexity curve of your library, with the number of reads in the first column and the expected distinct read value in the second column. For a typical experiment (human sample) check the expected complexity at 300M reads (to show the content of the file, type cat out.preseq). Expected unique pairs at 300M sequencing is at least ~ 120 million.
rule lib_complexity:
    input:
        bam="results/{sample}.mapped.PT.bam",
    output:
        complexity="results/{sample}.complexity.txt"
    conda:
        "envs/microc_env.yaml"
    shell:
        """
        preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output {output.complexity} {input.bam}
        """


rule contact_matrix:
    input:
        pairs="results/{sample}.mapped.pairs",
        genome_file=config["genome_file"],
    output:
        cmap="results/{sample}.contact_map.hic",
    conda:
        "envs/microc_env.yaml"
    threads: config["threads"]
    shell:
        """
        java -Xmx48000m  -Djava.awt.headless=true -jar ./scripts/juicer_tools_1.22.01.jar pre --threads {threads} {input.pairs} {output.cmap} {input.genome_file}
        """















# # run hicup pipeline 
# rule hicup:
#     input:
#        "raw-data/{sample}_1.fastq.gz",
#        "raw-data/{sample}_2.fastq.gz"
#     output:
#        "data/hicup/{sample}/{sample}_1_2.hicup.bam"
#     params:
#         index = config["GENOME"],
#         digest = config["DIGEST"]
#     threads: 16 
#     conda:
#         "envs/hic.yml"
#     log:
#         "data/logs/hicup_{sample}.log"
#     shell:
# #         "hicup --bowtie2 $(which bowtie2) "
#         "hicup --bowtie2 ~/miniconda3/envs/snakemake/bin/bowtie2 "
#         "--digest {params.digest} "
#         "--format Sanger "
#         "--index {params.index} "
#         "--longest 800 "
#         "--zip "
#         "--outdir data/hicup/{wildcards.sample} "
#         "--shortest 50 "
#         "--threads {threads} "
#         "{input} >{log} 2>&1"

# rule bam2pairs:
#     input:
#         "data/hicup/{sample}/{sample}_1_2.hicup.bam"
#     output:
#         "data/pairix/{sample}.bsorted.pairs.gz"
#     conda:
#         "envs/hic.yml"
#     log:
#         "data/logs/bam2pairs.{sample}.log"
#     shell:
#     # -c -p to uniqify @SQ and @PG
#         "bam2pairs -c {config[CHRSIZES]} {input} data/pairix/{wildcards.sample} >{log} 2>&1"

# # merge samples for best resolution = best quality boundary calls
# rule merge_pairs:
#     input:
#         expand("data/pairix/{sample}.bsorted.pairs.gz", sample=samples)
#     output:
#         "data/pairix/{merged_prefix}.pairs.gz"
#     conda:
#         "envs/hic.yml"
#     params:
#         prefix=merged_prefix
#     shell:
#         "merge-pairs.sh data/pairix/{params.prefix} {input}"

# # make coolers
# rule cooler_cload:
#     input:
# #         expand("data/pairix/{sample}.bsorted.pairs.gz", sample=samples)
#         "data/pairix/{sample}.bsorted.pairs.gz"
#     output:
#         "data/cool/{sample}.cool"
#     conda:
#         "envs/hic.yml"
#     log:
#         "data/logs/cooler_cload.{sample}.log"  # Remove expand and use the same wildcard
#     threads: 16
#     shell:
#         "cooler cload pairix -p {threads} --assembly {config[ASMBLY]} {config[CHRSIZES]}:10000 {input} {output} >{log} 2>&1"
   

# # create multiple resolution coolers for visualization in higlass
# rule zoomify:
#     input:
#         rules.cooler_cload.output
#     output:
#         "data/cool/{sample}.mcool"
#     log:
#         "data/logs/zoomify.{sample}.log"
#     conda:
#         "envs/hic.yml"
#     threads: 16 
#     shell:
#         "cooler zoomify -p {threads} --balance -o {output} {input} >{log} 2>&1"

# # find TADs using hicFindTADs from hicExplorer at 10kb resolution
# rule hicFindTADs:
#     input:
#         "data/cool/{sample}.mcool"
#     output:
#         "data/TADs/{sample}_min10_max60_fdr01_d01_boundaries.bed"
#     conda:
#         "envs/hicexplorer.yml"
#     log:
#         "data/logs/hicFindTADs_narrow.{sample}.log"
#     threads: 16
#     shell:
#         "hicFindTADs -m {input}::resolutions/10000 --minDepth 100000 --maxDepth 600000 --outPrefix data/TADs/{wildcards.sample}_min10_max60_fdr01_d01 --correctForMultipleTesting fdr -p {threads} >{log} 2>&1"

# rule multiqc:
#     input:
#        expand("data/hicup/{sample}/{sample}_1_2.hicup.bam", sample=samples), "data"
#     output:
#         "data/multiqc/multiqc_report.html"
#     log:
#         "data/logs/multiqc.log"
#     wrapper:
#         "v1.5.0/bio/multiqc"

