include: "code/metagenome/workflow.smk"

rule targets:
    input:
        "docs/proposal.pdf",
        expand("data/metagenome/bwa_IGC_results/{sample}_IGC.bam", sample=samples)

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
        code="code/download.sh",
        sra_list="data/SRR_Acc_List.txt"
    shell:
        "bash {input.code} {input.sra_list} data/raw/"
