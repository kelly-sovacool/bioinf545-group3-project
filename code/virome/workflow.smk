
rule virome:
    input:
        expand("data/qc/trimm_results/{sample}_repaired_{pair}.fastq.gz", sample=virome_samples, pair=[1,2])

rule remove_human:
    input:
        fq="data/qc/trimm_results/{sample}_repaired_{pair}.fastq.gz"