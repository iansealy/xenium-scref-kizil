library(tidyverse)
library(Seurat)

seurat_obj <- readRDS("wtPBSofabca7.rds")

write_tsv(
    as.data.frame(rownames(seurat_obj)),
    "gene_names.tsv",
    col_names = FALSE
)
