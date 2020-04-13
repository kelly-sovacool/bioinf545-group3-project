rule targets:
    input:
        "docs/proposal.pdf",
        "docs/report.pdf"

rule render_pdf:
    input:
        code="code/render.R",
        rmd="submission/{doc}.Rmd",
        preamble="submission/preamble_{doc}.tex"
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
