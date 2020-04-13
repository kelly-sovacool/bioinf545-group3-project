with open('data/metagenome/SRR_Acc_List_metagen.txt', 'r') as infile:
    samples = [line.strip() for line in infile]

rule trim:
    input:
        R1="data/fastq/{sample}_1.fastq.gz",
        R2="data/fastq/{sample}_2.fastq.gz"
    params:
        adapters="data/metagenome/metagenome_adapters.fasta"
    output:
        R1_P="data/metagenome/trimm_results/{sample}_paired_1.fastq.gz",
        R2_P="data/metagenome/trimm_results/{sample}_paired_2.fastq.gz",
        R1_U="data/metagenome/trimm_results/{sample}_unpaired_1.fastq.gz",
        R2_U="data/metagenome/trimm_results/{sample}_unpaired_2.fastq.gz"
    threads: 16
    log:
        "log/metagenome/trimmomatic_{sample}.log"
    benchmark:
        "benchmarks/metagenome/trimmomatic_{sample}.txt"
    shell:
        """
        trimmomatic PE -phred33 -threads {threads} \
                       {input.R1} {input.R2} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
                       ILLUMINACLIP:{params}:2:40:15 \
                       LEADING:3 TRAILING:3 MINLEN:24 SLIDINGWINDOW:4:15 \
        &> {log}
        """

rule re_pair:
    input:
        R1="data/metagenome/trimm_results/{sample}_paired_1.fastq.gz",
        R2="data/metagenome/trimm_results/{sample}_paired_2.fastq.gz"
    output:
        R1="data/metagenome/trimm_results/{sample}_repaired_1.fastq.gz",
        R2="data/metagenome/trimm_results/{sample}_repaired_2.fastq.gz",
        single="data/metagenome/trimm_results/{sample}_singleton.fastq.gz"
    conda:
        "../../environment_bwa.yml"
    log:
        "log/metagenome/repair_GRCh38_{sample}.log"
    benchmark:
        "benchmarks/metagenome/repair_{sample}.txt"
    shell:
        """
        repair.sh in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} outs={output.single} 2> {log}
        """

rule bwa_mem_GRCh38:
    input:
        R1="data/metagenome/trimm_results/{sample}_repaired_1.fastq.gz",
        R2="data/metagenome/trimm_results/{sample}_repaired_2.fastq.gz"
    params:
        index="data/metagenome/bwa_DB/GRCh38/GRCh38"
    output:
        bam="data/metagenome/bwa_GRCh38_results/{sample}_GRCh38.bam",
        flagstat="data/metagenome/bwa_GRCh38_results/{sample}_flagstat.txt"
    conda: 
       "../../environment_bwa.yml"
    threads: 16
    log:
        "log/metagenome/bwa-mem_GRCh38_{sample}.log"
    benchmark:
        "benchmarks/metagenome/bwa-mem_GRCh38_{sample}.txt"
    shell:
        """
        bwa mem -t {threads} {params.index} {input.R1} {input.R2} |
        samtools view -bh - > {output.bam} 2> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule metaphlan2:
    input:
        "data/metagenome/bwa_GRCh38_results/{sample}_GRCh38.bam"
    output:
        mtphln2="data/metagenome/metaphlan_test/{sample}_mtphln2.txt",
        bowtie2="data/metagenome/metaphlan_test/{sample}_bowtie2.out.bz2"
    threads: 4
    log:
        "log/metagenome/metaphlan2_{sample}.log"
    benchmark:
        "benchmarks/metagenome/metaphlan2_{sample}.txt"
    shell:
        """ 
        samtools fasta {input} | cat |
                 metaphlan2.py --input_type fasta --nproc {threads} --bowtie2out {output.bowtie2} > {output.mtphln2}
        2> {log}
        """

rule bam_to_fastq:
    input:
        "data/metagenome/bwa_GRCh38_results/{sample}_GRCh38.bam"
    output:
        R1="data/metagenome/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        R2="data/metagenome/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    conda: 
       "../../environment_bwa.yml"
    log:
        "log/metagenome/bamtofastq_{sample}.log"
    benchmark:
        "benchmarks/metagenome/bamtofastq_{sample}.txt"
    shell:
        """
        samtools sort -n {input} |
        samtools fastq -1 {output.R1} -2 {output.R2} - 2> {log}
        """                

rule bwa_mem_IGC:
    input:
        R1="data/metagenome/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        R2="data/metagenome/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    params:
        index="data/metagenome/bwa_DB/IGC/IGC"
    output:
        bam="data/metagenome/bwa_IGC_results/{sample}_IGC.bam",
        flagstat="data/metagenome/bwa_IGC_results/{sample}_flagstat.txt"
    conda: 
       "../../environment_bwa.yml"
    threads: 16
    log:
        "log/metagenome/bwa-mem_IGC_{sample}.log"
    benchmark:
        "benchmarks/metagenome/bwa-mem_IGC_{sample}.txt"
    shell:
        """
        bwa mem -t {threads} {params.index} {input.R1} {input.R2} |
        samtools view -bh - > {output.bam} 2> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule extract_geneList:
    input:
        "data/metagenome/bwa_IGC_results/{sample}_IGC.bam"
    params:
        "data/metagenome/bwa_DB/IGC.annotation_OF.summary"
    output:
        gene="data/metagenome/gene_abundance_results/{sample}.gene",
        list="data/metagenome/gene_abundance_results/{sample}.list",
        anno="data/metagenome/gene_abundance_results/{sample}_anno.txt",
        kegg="data/metagenome/gene_abundance_results/{sample}.kegg"
    conda: 
       "../../environment_bwa.yml"
    shell:
        """
        samtools view {input} |
        cut -f 3 - | sort - | uniq -c - | sort -b -nr -k 1,1 - | grep -v ":" - > {output.gene}
        sed -i 's/^ *//' {output.gene}
        cut -f 2 -d " " {output.gene} > {output.list}
        grep -Fw -f {output.list} {params} > {output.anno}
        cut -f 8 {output.anno} | grep -v "unknown" - | sort - | uniq -c - | sort -b -nr -k 1,1 - > {output.kegg}
        """

rule countKegg:
    input:
        "data/metagenome/gene_abundance_results/{sample}.gene"
    params:
        "data/metagenome/bwa_DB/IGC.kegg"
    output:
        "data/metagenome/gene_abundance_results/{sample}_keggCount.txt"
    conda:
        "../../environment.yml"
    log:
        "log/metagenome/countKegg_{sample}.log"
    shell:
        "python code/countKegg.py {input} {output} 2> {log}"
