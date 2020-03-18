rule targets:
    input:
        "docs/proposal.pdf"

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

rule render_slides:
    input:
        code="code/render.R",
        rmd="submission/{doc}.Rmd",
        preamble="submission/preamble.tex"
    output:
        file="docs/{doc}.html"
    params:
        format="ioslides_presentation"
    script:
        "{input.code}"
