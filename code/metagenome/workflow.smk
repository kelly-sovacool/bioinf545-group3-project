with open('data/metagenomes/SRR_Acc_List_metagen.txt', 'r') as infile:
    samples = [line.strip() for line in infile]

rule targets:
    input:
        expand("data/metagenome/bwa_IGC_results/{sample}_IGC.bam", sample=samples)

rule bwa_mem:
    input:
        R1="data/metagenome/bwa_GRC_sbatch/{sample}_1_unmapped.fastq.gz",
        R2="data/metagenome/bwa_GRC_sbatch/{sample}_2_unmapped.fastq.gz",
        ref='data/metagenome/bwa_DB/IGC/IGC'
    output:
        bam="data/metagenome/bwa_IGC_results/{sample}_IGC.bam",
        flagstat="data/metagenome/bwa_IGC_results/{sample}_flagstat.txt"
    threads: 16
    log:
        "log/metagenome/bwa-mem_{sample}.log"
    benchmark:
        "benchmarks/metagenome/bwa-mem_{sample}.txt"
    shell:
        """
        bwa mem -t {threads} {reference_file} {input.R1} {input.R2} |
        samtools view -Sb - > {output.bam} 2> {log}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """
