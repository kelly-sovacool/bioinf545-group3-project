rule virome_assembly:
    input:
        r1="data/qc/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        r2="data/qc/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    output:
        fna="data/virome/assembly/{sample}.megahit_asm/final.contigs.fa"
    params:
        dir="data/virome/assembly/{sample}.megahit_asm"
    threads: num_threads
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -o {params.dir} -t {threads} \
            --min-contig-len 1000 --k-min 21 --k-max 101 --k-step 20
        """

rule concat_contigs:
    input:
        expand("data/virome/assembly/{sample}.megahit_asm/final.contigs.fa", sample=virome_samples)
    output:
        fna="data/virome/contigs/contigs.fna"
    run:
        import Bio
        contigs = set()  # only keep unique sequences
        for infile in input:
            contigs.update({record.seq
                            for record in Bio.SeqIO.parse(infile, 'fasta')})
        Bio.SeqIO.write(contigs, output, 'fasta')

rule index_contigs:
    input:
        fna=rules.concat_contigs.output.fna
    output:
        bwt="data/virome/contigs/contigs.bwt"
    threads: num_threads
    shell:
        """
        bwa index -t {threads} {input.fna}
        """

rule align:
    input:
        fna=rules.virome_assembly.output.fna,
        bwt=rules.index_contigs.output.bwt
    output:
        bam="data/virome/align/{sample}_mapped.bam",
        tsv="data/virome/align/{sample}_flagstat.tsv",
        txt="data/virome/align/{sample}_idxstats.txt"
    params:
        index="data/virome/contigs/contigs"
    threads: num_threads
    shell:
        """
        bwa mem -t {threads} {params.index} {input.fna} |
        samtools sort - - |
        samtools view -bh - > {output.bam} 2> {log}
        samtools flagstat -O tsv {output.bam} > {output.tsv}
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.txt}
        """

rule contig_counts:
    input:
        code="code/virome/get_count_table.py",
        files=expand("data/virome/align/{sample}_idxstats.txt", sample=virome_samples)
    output:
        txt="data/virome/contig_counts.tsv"
    shell:
        """
        python {input.code} {input.files} > {output.txt}
        """