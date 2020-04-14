
with open('data/SRR_Acc_List.txt', 'r') as infile:
    all_samples = [line.strip() for line in infile]

include: "code/metagenome/workflow.smk"

rule targets:
    input:
        "docs/proposal.pdf",
        expand("data/metagenome/gene_abundance_results/{sample}_keggCount.txt", sample=samples)

rule render_pdf:
    input:
        code="code/render.R",
        rmd="submission/{doc}.Rmd",
        preamble="submission/preamble.tex"
    output:
        file="docs/{doc}.pdf"
    params:
        format="pdf_document"
    script:
        "{input.code}"

rule download:
    input:
        expand("data/raw/{sample}_{pair}.fastq.gz", sample=all_samples, pair=[1,2])

rule fastq_dump:
    input:
        sra_list="data/SRR_Acc_List.txt"
    output:
        "data/raw/{sample}_1.fastq.gz",
        "data/raw/{sample}_2.fastq.gz"
    params:
        outdir="data/raw/"
    shell:
        """
        prefetch {wildcards.sample}
        fasterq-dump --split-files {wildcards.sample} -O {params.outdir}
        gzip data/raw/{wildcards.sample}_1.fastq data/raw/{wildcards.sample}_2.fastq
        """
