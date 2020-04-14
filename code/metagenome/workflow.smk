

rule metaphlan2_samples:
    input:
        "data/qc/bwa_GRCh38_results/{sample}_GRCh38_unmapped.bam"
    output:
        mtphln2="data/metagenome/metaphlan2_samples/{sample}_mtphln2.txt",
        bowtie2="data/metagenome/metaphlan_samples/{sample}_bowtie2.out.bz2"
    threads: num_threads
    log:
        "log/metagenome/metaphlan2_{sample}.log"
    benchmark:
        "benchmarks/metagenome/metaphlan2_{sample}.txt"
    shell:
        """
        samtools view {input} | samtools fasta - | cat |
                 metaphlan2.py --input_type fasta --nproc {threads} --bowtie2out {output.bowtie2} > {output.mtphln2}
        2> {log}
        """

rule metaphlan2_results:
    input:
        expand("data/metagenome/metaphlan2_samples/{sample}_mtphln2.txt", sample=metag_samples)
    output:
        merged="data/metagenome/metaphlan2_results/merged.txt",
        phylum="data/metagenome/metaphlan2_results/merged_phylum.txt",
        family="data/metagenome/metaphlan2_results/merged_family.txt",
        genus="data/metagenome/metaphlan2_results/merged_genus.txt",
        species="data/metagenome/metaphlan2_results/merged_species.txt"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output.merged}
        grep -E "(s__)|(^clade_name)" {output.merged} | grep -v "t__" | sed 's/^.*s__//g' > {output.species}
        grep -E "(g__)|(^clade_name)" {output.merged} | grep -v "t__" | sed 's/^.*g__//g' | sed '/|s_/d' > {output.genus}
        grep -E "(f__)|(^clade_name)" {output.merged} | grep -v "t__" | sed 's/^.*f__//g' | sed '/|g_/d' > {output.family}
        grep -E "(p__)|(^clade_name)" {output.merged} | grep -v "t__" | sed 's/^.*p__//g' | sed '/|c_/d' > {output.phylum}
        """

rule bwa_mem_IGC:
    input:
        R1="data/qc/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        R2="data/qc/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    params:
        index="data/metagenome/bwa_DB/IGC/IGC"
    output:
        bam="data/metagenome/bwa_IGC_results/{sample}_IGC.bam",
        flagstat="data/metagenome/bwa_IGC_results/{sample}_flagstat.txt"
    conda:
       "../../environment_bwa.yml"
    threads: num_threads
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
        "data/qc/bwa_GRCh38_results/{sample}_GRCh38_unmapped.bam"
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
        samtools view -f 2 {input} |
        cut -f 3 - | sort - | uniq -c - | sort -b -nr -k 1,1 - | grep -v ":" - > {output.gene}
        sed -i 's/^ *//' {output.gene}
        cut -f 2 -d " " {output.gene} > {output.list}
        grep -Fw -f {output.list} {params} > {output.anno}
        cut -f 8 {output.anno} | grep -v "unknown" - | sort - | uniq -c - | sort -b -nr -k 1,1 - > {output.kegg}
        """

rule countKegg:
    input:
        rules.extract_geneList.output.gene
    params:
        "data/metagenome/bwa_DB/IGC.kegg"
    output:
        "data/metagenome/gene_abundance_results/{sample}_keggCount.txt"
    log:
        "log/metagenome/countKegg_{sample}.log"
    shell:
        """
        python code/metagenome/countKegg.py {input} {output} 2> {log}
        code/metagenome/clean_countKegg.sh {output}
        """
