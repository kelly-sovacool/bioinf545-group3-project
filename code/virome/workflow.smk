with open('data/virome/SRR_Acc_List_virome.txt', 'r') as infile:
    samples = [line.strip() for line in infile]

rule trim:
    input:
        R1="data/raw/{sample}_1.fastq.gz",
        R2="data/raw/{sample}_2.fastq.gz"
    params:
        adapters="data/virome/virome_adapters.fasta"
    output:
        R1_P="data/virome/trimm_results/{sample}_paired_1.fastq.gz",
        R2_P="data/virome/trimm_results/{sample}_paired_2.fastq.gz",
        R1_U="data/virome/trimm_results/{sample}_unpaired_1.fastq.gz",
        R2_U="data/virome/trimm_results/{sample}_unpaired_2.fastq.gz"
    threads: 16
    log:
        "log/virome/trimmomatic_{sample}.log"
    benchmark:
        "benchmarks/virome/trimmomatic_{sample}.txt"
    shell:
        """
        trimmomatic PE -phred33 -threads {threads} \
                       {input.R1} {input.R2} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
                       ILLUMINACLIP:{params}:2:40:15 \
                       LEADING:3 TRAILING:3 MINLEN:24 SLIDINGWINDOW:4:15 \
        &> {log}
        """