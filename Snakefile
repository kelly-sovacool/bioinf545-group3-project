configfile: "config.yml"
num_threads = config['threads']

with open('data/SRR_Acc_List.txt', 'r') as infile:
    all_samples = [line.strip() for line in infile]
with open('data/metagenome/SRR_Acc_List_metagen.txt', 'r') as infile:
    metag_samples = [line.strip() for line in infile]
with open('data/virome/SRR_Acc_List_virome.txt', 'r') as infile:
    virome_samples = [line.strip() for line in infile]

include: "code/metagenome/workflow.smk"
include: "code/virome/workflow.smk"

rule targets:
    input:
        "docs/proposal.pdf",
        expand("data/metagenome/gene_abundance_results/{sample}_keggCount.txt", sample=metag_samples),
        "data/metagenome/metaphlan2_results/merged.txt",
        expand("data/qc/bwa_GRCh38_results/{sample}_unmapped_{pair}.fastq.gz", sample=all_samples, pair=[1,2])

rule render_pdf:
    input:
        code="code/render.R",
        rmd="submission/{doc}.Rmd",
        preamble="submission/preamble.tex"
    output:
        file="docs/{doc}.pdf"
    params:
        format="pdf_document"
    script:
        "{input.code}"

rule download:
    input:
        expand("data/raw/{sample}_{pair}.fastq.gz", sample=all_samples, pair=[1,2])

rule fastq_dump:
    input:
        sra_list="data/SRR_Acc_List.txt"
    output:
        "data/raw/{sample}_1.fastq.gz",
        "data/raw/{sample}_2.fastq.gz"
    params:
        outdir="data/raw/"
    shell:
        """
        prefetch {wildcards.sample}
        fasterq-dump --split-files {wildcards.sample} -O {params.outdir}
        gzip data/raw/{wildcards.sample}_1.fastq data/raw/{wildcards.sample}_2.fastq
        """

rule trim:
    input:
        R1="data/raw/{sample}_1.fastq.gz",
        R2="data/raw/{sample}_2.fastq.gz"
    params:
        adapters="data/qc/adapters.fna"
    output:
        R1_P="data/qc/trimm_results/{sample}_paired_1.fastq.gz",
        R2_P="data/qc/trimm_results/{sample}_paired_2.fastq.gz",
        R1_U="data/qc/trimm_results/{sample}_unpaired_1.fastq.gz",
        R2_U="data/qc/trimm_results/{sample}_unpaired_2.fastq.gz"
    threads: num_threads
    log:
        "log/qc/trimmomatic_{sample}.log"
    benchmark:
        "benchmarks/qc/trimmomatic_{sample}.txt"
    shell:
        """
        trimmomatic PE -phred33 -threads {threads} \
                       {input.R1} {input.R2} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
                       ILLUMINACLIP:{params.adapters}:2:40:15 \
                       LEADING:3 TRAILING:3 MINLEN:24 SLIDINGWINDOW:4:15 \
        &> {log}
        """

rule re_pair:
    input:
        R1=rules.trim.output.R1_P,
        R2=rules.trim.output.R2_P
    output:
        R1="data/qc/trimm_results/{sample}_repaired_1.fastq.gz",
        R2="data/qc/trimm_results/{sample}_repaired_2.fastq.gz",
        single="data/qc/trimm_results/{sample}_singleton.fastq.gz"
    conda:
        "../../environment_bwa.yml"
    log:
        "log/qc/repair_GRCh38_{sample}.log"
    benchmark:
        "benchmarks/qc/repair_{sample}.txt"
    shell:
        """
        repair.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} outs={output.single} 2> {log}
        """

rule bwa_mem_GRCh38:
    input:
        R1=rules.re_pair.output.R1,
        R2=rules.re_pair.output.R1
    params:
        index="data/qc/bwa_DB/GRCh38/GRCh38"
    output:
        mapped="data/qc/bwa_GRCh38_results/{sample}_GRCh38_mapped.bam",
        unmapped="data/qc/bwa_GRCh38_results/{sample}_GRCh38_unmapped.bam",
        flagstat="data/qc/bwa_GRCh38_results/{sample}_flagstat.txt"
    conda:
       "environment_bwa.yml"
    threads: num_threads
    log:
        "log/qc/bwa-mem_GRCh38_{sample}.log"
    benchmark:
        "benchmarks/qc/bwa-mem_GRCh38_{sample}.txt"
    shell:
        """
        bwa mem -t {threads} {params.index} {input.R1} {input.R2} |
        samtools view -bh - > {output.unmapped} 2> {log}
        samtools view -bh -f 8 {output.unmapped} > {output.unmapped}
        samtools flagstat {output.unmapped} > {output.flagstat}
        """

rule bam_to_fastq:
    input:
        rules.bwa_mem_GRCh38.output.unmapped
    output:
        R1="data/qc/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        R2="data/qc/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    conda:
       "environment_bwa.yml"
    log:
        "log/qc/bamtofastq_{sample}.log"
    benchmark:
        "benchmarks/qc/bamtofastq_{sample}.txt"
    shell:
        """
        samtools sort -n {input} |
        samtools fastq -1 {output.R1} -2 {output.R2} - 2> {log}
        """
