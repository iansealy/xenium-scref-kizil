# Single Cell Reference for Xenium Custom Panel Design

Obtained wtPBSofabca7.rds file (a Seurat object) from Elanur Yilmaz & Caghan Kizil

Seurat object doesn't contain Ensembl IDs (only gene names), but https://www.10xgenomics.com/analysis-guides/creating-single-cell-references-for-xenium-custom-panel-design-from-seurat-or-anndata requires them

Get gene names from Seurat object:

```
Rscript get_gene_names.R
```

Get Ensembl IDs and gene names for all Ensembl releases:

```
script/dump_ensembl_genes_wrapper.sh
```

Find each gene name in the Ensembl files:

```
split -d -a 1 -n 8 gene_names.tsv chunk_
ls chunk_* | parallel --jobs 8 "while read -r name; do grep \"\t\$name\$\" data/*/ensembl.txt; done < {} > {}.out"
cat chunk_?.out | sed -e 's/:/\t/' > gene_names_ens.tsv
rm chunk_*
```

Check which Ensembl release matches best:

```
cut -f1 gene_names_ens.tsv | sort | uniq -c | sort -nr | head
23458 data/104/ensembl.txt
23436 data/105/ensembl.txt
23434 data/102/ensembl.txt
23430 data/107/ensembl.txt
23427 data/106/ensembl.txt
23423 data/108/ensembl.txt
23418 data/103/ensembl.txt
23417 data/111/ensembl.txt
23415 data/100/ensembl.txt
23413 data/101/ensembl.txt
```

Map gene names to Ensembl 104 IDs:

```
cat /dev/null > gene_names_ens_104_meta.tsv
cat /dev/null > gene_names_not_found_ens_104.tsv
cat /dev/null > gene_names_multiple_ens_104.tsv
for name in `cat gene_names.tsv`; do
  count=`grep -c "\t$name$" data/104/ensembl.txt`
  if [ "$count" -eq 1 ]; then
    grep "\t$name$" data/104/ensembl.txt
    grep "\t$name$" data/104/ensembl.txt | awk -F"\t" '{ print $2 "\t" $1 "\tTRUE" }' >> gene_names_ens_104_meta.tsv
  elif echo "$name" | grep -q "^ENSDARG"; then
    echo -e "$name\t$name"
    echo -e "$name\t$name\tTRUE" >> gene_names_ens_104_meta.tsv
  elif [ "$count" -eq 0 ]; then
    echo "$name" >> gene_names_not_found_ens_104.tsv
    echo -e "$name\t\tFALSE" >> gene_names_ens_104_meta.tsv
  else
    echo "$name" >> gene_names_multiple_ens_104.tsv
    echo -e "$name\t\tFALSE" >> gene_names_ens_104_meta.tsv
  fi
done > gene_names_ens_104.tsv
```

Count genes:

```
wc -l gene_names.tsv gene_names_ens_104_meta.tsv gene_names_ens_104.tsv gene_names_multiple_ens_104.tsv gene_names_not_found_ens_104.tsv
   24122 gene_names.tsv
   24122 gene_names_ens_104_meta.tsv
   23268 gene_names_ens_104.tsv
     387 gene_names_multiple_ens_104.tsv
     467 gene_names_not_found_ens_104.tsv
```

So of the 24122 names, 23268 can be uniquely mapped to an Ensembl ID, with the remaining 854 mapping to multiple Ensembl IDs or being simply missing

Most of the missing gene names have a numerical suffix attached so actually probably map to multiple Ensembl IDs

Check the list of genes supplied by Amanda:

```
sed -e 's/.*,//' gene_list.csv | grep ENSDARG | sort -u | wc -l
336
```

```
sed -e 's/.*,//' gene_list.csv | grep ENSDARG | grep -c -f - gene_names_ens_104.tsv
323
```

So 13 genes missing:

```
for gene in `sed -e 's/.*,//' gene_list.csv | grep ENSDARG`; do echo -ne "$gene\t"; grep -c $gene gene_names_ens_104.tsv; done | grep 0$ | cut -f1 | grep -f - gene_list.csv | sort
crhr2,ENSDARG00000103704
ebf3a,ENSDARG00000099849
foxg1c,ENSDARG00000068380
her4.2,ENSDARG00000094426
meis2a,ENSDARG00000098240
meis2b,ENSDARG00000077840
npvf,ENSDARG00000036227
pmch,ENSDARG00000073959
pth1b,ENSDARG00000091961
pth2,ENSDARG00000022951
vav3b,ENSDARG00000075962
vipb,ENSDARG00000079443
wnt8a,ENSDARG00000078507
```

6 can't be mapped uniquely to a single Ensembl ID (crhr2, ebf3a, her4.2, meis2a, meis2b, vav3b) and 7 are missing entirely from the annotation used for the single cell data (foxg1c, npvf, pmch, pth1b, pth2, vipb)

```
for gene in `sed -e 's/.*,//' gene_list.csv | grep ENSDARG`; do echo -ne "$gene\t"; grep -c $gene gene_names_ens_104.tsv; done | grep 0$ | cut -f1 \
  | grep -f - gene_list.csv | sed -e 's/,ENSDARG.*//' | sed -e 's/.*,//' | grep -wif - gene_names_ens_104_meta.tsv | sort
crhr2		FALSE
crhr2.1		FALSE
ebf3a		FALSE
ebf3a.1		FALSE
her4.2		FALSE
her4.2.1		FALSE
meis2a		FALSE
meis2a.1		FALSE
vav3b		FALSE
vav3b.1		FALSE
```

In R:

```
library(tidyverse)
library(Seurat)

seurat_obj <- readRDS("wtPBSofabca7.rds")

# Subset to genes with Ensembl IDs
ens.meta <- read_tsv("gene_names_ens_104_meta.tsv", col_names=c("gene_name", "ensembl_id", "ensembl_mapped"))
genes.to.keep <- ens.meta$gene_name[ens.meta$ensembl_mapped]
seurat_subset_obj <- subset(seurat_obj, features = genes.to.keep)

# Function for creating MEX files
writeCounts <- function(out_dir, data, barcodes = colnames(data), gene.id = rownames(data), gene.symbol = rownames(data), feature.type = "Gene Expression", subset = 1:length(barcodes)) {
  require("R.utils")
  require("Matrix")

  if (file.exists(out_dir) || (dir.exists(out_dir) && length(list.files(out_dir)) > 0)) {
    stop("The specified output directory already exists! Not overwriting")
  }
  dir.create(out_dir, recursive = TRUE)

  write.table(
    data.frame(barcode = barcodes[subset]),
    gzfile(file.path(out_dir, "barcodes.tsv.gz")),
    sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )

  write.table(
    data.frame(
      gene_id = gene.id,
      gene_symbol = gene.symbol,
      feature_type = feature.type
    ),
    gzfile(file.path(out_dir, "features.tsv.gz")),
    sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )

  Matrix::writeMM(data[, subset], file.path(out_dir, "matrix.mtx"))
  R.utils::gzip(file.path(out_dir, "matrix.mtx"), remove = TRUE)
}

# Create MEX files
writeCounts(
  "reference_data",
  GetAssayData(seurat_subset_obj, assay="RNA", slot="counts"),
  gene.id = ens.meta$ensembl_id[ens.meta$ensembl_mapped],
  gene.symbol = rownames(seurat_subset_obj)
)

# Function for saving cell type annotation and bundling reference data
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  data.table::fwrite(
    data.table::data.table(
      barcode = barcodes,
      annotation = unlist(data[[cell_type]])
    )[subset, ],
    file.path(out_dir, "annotations.csv")
  )

  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))

  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )

  if (file.info(bundle)$size / 1e6 > 2000) {
    warning("The output file is more than 2G and will need to be subset further.")
  }
}

# Create cell type annotation and bundle reference data
bundleOutputs(out_dir = "reference_data", data = seurat_subset_obj, cell_type="main.cells")
```
