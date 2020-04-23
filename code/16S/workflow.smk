include: "references.smk"

rule make_names_file:
    input:
        code="code/16S/make_names_file.R",
        metadata="data/SraRunTable.txt",
        fasta=expand("data/raw/{sample}_{pair}", sample=metag_samples, pair=[1,2])
    output:
        file="data/16S/crc.files"
    shell:
        "Rscript {input.code}"

rule process_seqs:
    input:
        names=rules.make_names_file.output.file,
        silva="data/16S/refs/silva.v4.fasta",
        rdp_fna="data/16S/refs/rdp.bacteria.fasta",
        rdp_tax="data/16S/refs/rdp.bacteria.tax",
        mock="data/16S/refs/HMP_MOCK.v35.fasta"
    output:
        fasta="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta",
        count="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table",
        tax="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy"
    log:
        "log/16S/process_seqs.log"
    threads: num_threads
    params:
        indir="data/raw/",
        outdir="data/16S/processed/"
    shell:
        """
        mothur '#
        set.logfile(name={log});
        set.dir(input={params.indir}, output={params.outdir});
        make.contigs(file={input.names}, inputdir={params.indir}, processors={threads});
        summary.seqs(fasta=current);
        screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275);
        summary.seqs(fasta=current);
        unique.seqs(fasta=current);
        summary.seqs(fasta=current, name=current);
        count.seqs(name=current, group=current);
        summary.seqs(fasta=current, count=current);
        align.seqs(fasta=current, reference={input.silva});
        summary.seqs(fasta=current, count=current);
        screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8);
        summary.seqs(fasta=current, count=current);
        filter.seqs(fasta=current, vertical=T, trump=.);
        unique.seqs(fasta=current, count=current);
        pre.cluster(fasta=current, count=current, diffs=2);
        chimera.uchime(fasta=current, count=current, dereplicate=t);
        remove.seqs(fasta=current, accnos=current);
        summary.seqs(fasta=current, count=current);
        classify.seqs(fasta=current, count=current, reference={input.rdp_fna}, taxonomy={input.rdp_tax}, cutoff=80);
        remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
        get.groups(fasta=current, count=current, groups=mock1-mock2);
        seq.error(fasta=current, count=current, reference={input.mock}, aligned=F)
        '
        """

rule cluster_OTUs:
    input:
        fasta=rules.process_seqs.output.fasta,
        count=rules.process_seqs.output.count,
        tax=rules.process_seqs.output.tax
    output:
        dist="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist",
        list="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list",
        steps="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps",
        sensspec="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sensspec",
        shared="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
        tax="data/16S/processed/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
    log:
        "log/16S/cluster.log"
    params:
        workdir="data/16S/processed/",
        seed=545
    threads: num_threads
    shell:
        """
        mothur '#
        set.logfile(name={log});
        set.dir(input={params.workdir}, output={params.workdir});
        set.seed(seed={params.seed});
        dist.seqs(fasta={input.fasta}, cutoff=0.03, processors={threads});
        cluster(column=current, count={params.count}, cutoff=0.03);
        make.shared(list=current, count={input.count}, label=0.03);
        classify.otu(list=current, count=current, taxonomy={input.tax}, label=0.03);
        get.oturep(fasta={input.fasta}, count=current, list=current, label=0.03, method=abundance);
        remove.groups(shared=current, groups=mock1-mock2)
        '
        """

rule mothur_16S:
    input:
        expand("data/raw/{sample}_{r}.fastq.gz", sample=all_samples, r=[1,2])
    output:
        "data/16S/processed/tmp.txt"
    shell:
        """
        # real mothur pipeline was run interactively, not incorporated in Snakemake
        touch {output}
        """

rule cluster_OTU:
    input:
        rules.mothur_16S.output
    output:
        "data/16S/mothur_output/stability.opti_mcc.shared"
    shell:
        """
        # real mothur pipeline was run interactively, not incorporated in Snakemake
        touch {output}
        """