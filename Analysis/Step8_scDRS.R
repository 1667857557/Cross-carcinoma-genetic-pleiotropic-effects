#####
data <- Read10X_h5("TabulaTIME_integrated_expression.h5")
meta<-readRDS("TabulaTIME_metacell_meta_information.rds")
colnames(meta)
A <- CreateSeuratObject(counts = data,meta.data = meta)

umap_matrix <- as.matrix(A@meta.data[, c("UMAP1", "UMAP2")])

new_umap_reduction <- CreateDimReducObject(
  embeddings = umap_matrix, 
  key = "UMAP_",           
  assay = "RNA"            
)
A[["umap"]] <- new_umap_reduction

non_epi <- c("GIST","PBMC","PPB","SKCM","SS","OS","UVM")
cat("cells before:", ncol(A), "\n")
A <- subset(A, subset = !(Cancer_type %in% non_epi) & !is.na(Cancer_type))
cat("cells after: ", ncol(A), "\n")
table(A$Cancer_type)
A@assays$RNA@data<-A@assays$RNA@counts
sceasy::convertFormat(A, from="seurat", to="anndata",
                      outFile='metacell_pancarcinoma.h5ad')

######
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils","Seurat")
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
setwd("D:/")
A<-fread("F3.genes.out")
m <- useMart("ensembl","hsapiens_gene_ensembl")
res <- getBM(c("ensembl_gene_id","hgnc_symbol"),"ensembl_gene_id",A$GENE,m)
A <- merge(A,res,by.x="GENE",by.y="ensembl_gene_id",all.x=TRUE)
A$hgnc_symbol <- ifelse(is.na(A$hgnc_symbol), A$GENE, A$hgnc_symbol)
A <- A[!duplicated(A$hgnc_symbol), ]
A<-dplyr::select(A,c(hgnc_symbol,P))
colnames(A)<-c("GENE","P")
write_tsv(A,file = "F3.genes.txt")

scdrs munge-gs \
--out-file F3.gs \
--pval_file F3.genes.txt \
--weight zscore \
--n-max 1000

scdrs compute-score \
--h5ad-file metacell_pancarcinoma.h5ad \
--h5ad-species human \
--gs-file ALL.gs \
--gs-species human \
--out-folder /mnt/d/scDRS/ \
--flag-filter-data True \
--flag-raw-count True \
--n-ctrl 1000 \
--flag-return-ctrl-raw-score False \
--flag-return-ctrl-norm-score True


scdrs perform-downstream \
--h5ad-file metacell_pancarcinoma.h5ad \
--score-file /mnt/d/scDRS/EF.full_score.gz \
--out-folder  /mnt/d/scDRS/ \
--group-analysis curated_anno \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

scdrs perform-downstream \
--h5ad-file metacell_pancarcinoma.h5ad \
--score-file /mnt/d/scDRS/CF.full_score.gz \
--out-folder  /mnt/d/scDRS/ \
--group-analysis curated_anno \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

scdrs perform-downstream \
--h5ad-file metacell_pancarcinoma.h5ad \
--score-file /mnt/d/scDRS/F1.full_score.gz \
--out-folder  /mnt/d/scDRS/ \
--group-analysis curated_anno \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

scdrs perform-downstream \
--h5ad-file metacell_pancarcinoma.h5ad \
--score-file /mnt/d/scDRS/F2.full_score.gz \
--out-folder  /mnt/d/scDRS/ \
--group-analysis curated_anno \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

scdrs perform-downstream \
--h5ad-file metacell_pancarcinoma.h5ad \
--score-file /mnt/d/scDRS/F3.full_score.gz \
--out-folder  /mnt/d/scDRS/ \
--group-analysis curated_anno \
--gene-analysis \
--flag-filter-data True \
--flag-raw-count True

