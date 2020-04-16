rule make_names_file:
    input:
        code="code/16S/make_names_file.R",
        metadata="data/SraRunTable.txt"
    output:
        file="data/16S/crc.files"
    script:
        "{input.code}"

rule get_silva:

rule get_rdp:

rule get_hmp:

rule preprocess:
    input:
        script="code/16S/1_preprocess.sh"
    shell:
        "bash {input.script}"