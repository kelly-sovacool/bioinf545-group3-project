""" Download the reference databases and process with mothur """

rule download_silva:
    output:
        tar=temp("data/16S/refs/silva.nr_v138.tgz"),
        fasta="data/16S/refs/silva.fasta",
        tax="data/16S/refs/silva.tax"
    params:
        fasta="data/16S/refs/silva.nr_v138.align",
        tax="data/16S/refs/silva.nr_v138.tax",
        workdir="data/16S/refs/"
    benchmark:
        "benchmarks/silva/download.txt"
    shell:
        """
        wget -N -P {params.workdir} https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138.tgz
        tar xzvf {output.tar} -C {params.workdir}
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        """

rule download_rdp:
    output:
        tar=temp("data/16S/refs/trainset16_022016.pds.tgz"),
        fasta="data/16S/refs/rdp.fasta",
        tax="data/16S/refs/rdp.tax"
    params:
        fasta="data/16S/refs/trainset16_022016.pds/trainset16_022016.pds.fasta",
        tax="data/16S/refs/trainset16_022016.pds/trainset16_022016.pds.tax"
    benchmark:
        "benchmarks/rdp/download.txt"
    shell:
        """
        wget -N -P data/16S/refs/ https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz
        tar xvzf {output.tar} -C data/16S/refs/
        mv {params.fasta} {output.fasta}
        mv {params.tax} {output.tax}
        rm -rf data/16S/refs/trainset16_022016.pds/
        """

rule download_mock:
    output:
        "data/16S/refs/HMP_MOCK.v35.fasta"
    params:
        workdir="data/16S/refs/"
    shell:
        """
        wget -N -P {params.workdir} https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
        unzip -d {params.workdir} miseqsopdata.zip MiSeq_SOP/HMP_MOCK.v35.fasta
        """

rule get_bacteria:
    input:
        fasta="data/16S/refs/{ref}.fasta",
        tax="data/16S/refs/{ref}.tax"
    output:
        fasta="data/16S/refs/{ref}.bacteria.fasta",
        tax="data/16S/refs/{ref}.bacteria.tax"
    params:
        fasta="data/16S/refs/{ref}.pick.fasta",
        tax="data/16S/refs/{ref}.pick.tax"
    log:
        "log/16S/{ref}_get_bacteria.log"
    benchmark:
        "benchmarks/16S/{ref}_get_bacteria.txt"
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(input=data/{wildcards.ref}/);

        get.lineage(fasta={input.fasta}, taxonomy={input.tax}, taxon=Bacteria);
        summary.seqs(fasta=current);

        rename.file(input={params.fasta}, new={output.fasta});
        rename.file(input={params.tax}, new={output.tax}) '
        """

rule subset_v4_region:
    input:
        align="data/16S/refs/silva.bacteria.fasta"
    output:
        sum="data/16S/refs/silva.bacteria.summary",
        fasta_v4="data/16S/refs/silva.v4.fasta"
    params:
        screen="data/16S/refs/silva.bacteria.good.pcr.align",
        workdir="data/16S/refs/"
    threads: num_threads
    log:
        "log/16S/silva_subset_v4_region.log"
    benchmark:
        "benchmarks/16S/silva_subset_v4_region.txt"
    shell:
        """
        mothur '#set.logfile(name={log}); set.dir(input={params.workdir});
        summary.seqs(fasta={input.align}, processors={threads});

        screen.seqs(fasta={input.align}, maxambig=0, maxhomop=8, start=11894, end=25319);
        summary.seqs();
        pcr.seqs(fasta=current, start=11894, end=25319, keepdots=F);
        summary.seqs();
        rename.file(input={params.screen}, new={output.fasta_v4})'
        """
