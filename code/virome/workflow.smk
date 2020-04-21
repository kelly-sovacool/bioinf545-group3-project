rule virome_assembly:
    input:
        r1="data/qc/bwa_GRCh38_results/{sample}_unmapped_1.fastq.gz",
        r2="data/qc/bwa_GRCh38_results/{sample}_unmapped_2.fastq.gz"
    output:
        fna="data/virome/assembly/{sample}.megahit_asm/final.contigs.fa"
    log:
        "log/virome/assembly/{sample}.txt"
    params:
        dir="data/virome/assembly/{sample}.megahit_asm"
    threads: num_threads
    shell:
        """
        rm -rf {params.dir}
        megahit -1 {input.r1} -2 {input.r2} -o {params.dir} -t {threads} \
            --min-contig-len 1000 --k-min 21 --k-max 101 --k-step 20 \
            2> {log}
        """

rule concat_contigs:
    input:
        fastas=expand("data/virome/assembly/{sample}.megahit_asm/final.contigs.fa", sample=virome_samples)
    output:
        fna="data/virome/contigs/contigs.fna"
    run:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        contig_seqs = set()  # only keep unique sequences
        contig_records = list()
        for infile in input.fastas:
            for record in SeqIO.parse(infile, 'fasta'):
                if record.seq not in contig_seqs:
                    contig_records.append(record)
                    contig_seqs.add(record.seq)
        SeqIO.write(contig_records, output.fna, 'fasta')

rule index_contigs:
    input:
        fna=rules.concat_contigs.output.fna
    output:
        bwt="data/virome/contigs/contigs.bwt"
    conda:
        "../../environment_bwa.yml"
    params:
        index="data/virome/contigs/contigs"
    shell:
        """
        bwa index -p {params.index} {input.fna}
        """

rule map:
    input:
        fna=rules.virome_assembly.output.fna,
        bwt=rules.index_contigs.output.bwt
    output:
        map="data/virome/mapping/{sample}_mapped.bam",
        bam="data/virome/mapping/{sample}_mapped.sorted.bam",
        tsv="data/virome/mapping/{sample}_flagstat.tsv",
        txt="data/virome/mapping/{sample}_idxstats.txt"
    log:
        "log/virome/mapping/{sample}.txt"
    params:
        index="data/virome/contigs/contigs"
    conda:
        "../../environment_bwa.yml"
    threads: num_threads
    shell:
        """
        bwa mem -t {threads} {params.index} {input.fna} |
        samtools view -bh - > {output.map} 2> {log}
        samtools sort {output.map} -o {output.bam}
        samtools flagstat -O tsv {output.bam} > {output.tsv}
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.txt}
        """

rule concoct_prep:
    input:
        contigs=rules.concat_contigs.output.fna,
        bams=expand("data/virome/mapping/{sample}_mapped.sorted.bam", sample=virome_samples)
    output:
        bed="rules/virome/contigs/contigs_10K.bed",
        fna="rules/virome/contigs/contigs_10K.fna",
        tsv="rules/virome/contigs/coverage_table.tsv"
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 > {output.fna}
        concoct_coverage_table.py {output.bed} {input.bams} > {output.tsv}
        """

rule concoct_cluster:
    input:
        tsv=rules.concoct_prep.output.tsv,
        fna=rules.concoct_prep.output.fna
    output:
        csv1="data/virome/concoct/clustering_gt1000.csv",
        csv2="data/virome/concoct/clustering_merged.csv"
    params:
        dir="data/virome/concoct/"
    threads: num_threads
    shell:
        """
        concoct --coverage_file {input.tsv} --composition_file = {input.fna} \
            --clusters 500 --kmer_length 4 --length_threshold 1000 \
            --read_length 150 --basename {params.dir} --no_total_coverage \
            --iterations 50 -t {threads} --seed 545
        merge_cutup_clustering.py {output.csv1} > {output.csv2}
        """

rule extract_viromes:
    input:
        fna=rules.concat_contigs.output.fna,
        csv=rules.concoct_cluster.output.csv2
    output:
        directory("data/virome/concoct/fasta_bins")
    shell:
        """
        extract_fasta_bins.py {input.fna} {input.csv} --output_path
        """

# TODO: Get OVU abundance
