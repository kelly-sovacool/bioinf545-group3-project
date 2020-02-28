format <- snakemake@params[["format"]]
input_file <- here::here(snakemake@input[["rmd"]])
output_file <- here::here(snakemake@output[["file"]])

# TODO: dictionary/list to map output file extensions to RMarkdown formats
rmarkdown::render(input_file,
  output_format = format,
  output_file = output_file
)
