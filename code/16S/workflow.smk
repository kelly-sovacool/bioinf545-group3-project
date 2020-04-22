include: "references.smk"

rule make_names_file:
    input:
        code="code/16S/make_names_file.R",
        metadata="data/SraRunTable.txt"
    output:
        file="data/16S/crc.files"
    script:
        "{input.code}"

rule process_seqs:
    input:
        names=rules.make_names_file.output.file,
        silva="data/16S/refs/silva.v4.fasta",
        rdp_fna="data/16S/refs/rdp.bacteria.fasta",
        rdp_tax="data/16S/refs/rdp.bacteria.tax",
        mock="data/16S/refs/HMP_MOCK.v35.fasta"
    log: "log/16S/process_seqs.log"
    threads: num_threads
    params:
        indir="data/raw/",
        outdir="data/16S/processed/"
    shell:
        """
        mothur '#
        set.logfile(name={log});
        set.dir(input={params.indir}, output={params.outdir});
        make.contigs(file={input.names}, inputDir={params.indir}, processors={threads});
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
        summary.seqs(fasta=current,count=current);
        filter.seqs(fasta=current, vertical=T, trump=.);
        unique.seqs(fasta=current, count=current);
        pre.cluster(fasta=current, count=current, diffs=2);
        chimera.uchime(fasta=current, count=current, dereplicate=t);
        remove.seqs(fasta=current, accnos=current);
        summary.seqs(fasta=current,count=current);
        classify.seqs(fasta=current, count=current, reference={input.rdp_fna}, taxonomy={input.rdp_tax}, cutoff=80);
        remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
        get.groups(fasta=current, count=current, groups=mock1-mock2);
        seq.error(fasta=current, count=current, reference={input.mock}, aligned=F)
        '
        """
