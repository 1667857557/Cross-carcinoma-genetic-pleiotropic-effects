setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")
tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"# 
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
set.seed(1000)
setwd("G:/Database/LAVA")

A<-fread("LAVA_s2500_m25_f1_w200.blocks")
A$LOC<-row.names(A)
colnames(A)<-c("CHR","START","STOP","LOC")

A[, LOC := sprintf("locus_%04d", .I)]

library(data.table)
library(LAVA)
load('D:/pan_cancer/ldsc.covstruct_ALL.Rdata')
phenos <- c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA")
R = round(cov2cor(ldsc.covstruct$I),5)
rownames(R) <- colnames(R) <- phenos     
write.table(R,
            file = "sample.overlap.txt",
            sep  = "\t",
            quote= FALSE,
            col.names = NA,   
            row.names = TRUE)
input <- process.input(
  input.info.file      = "G:/Database/LAVA/input.info_all.txt",
  sample.overlap.file  = "./sample.overlap_all.txt",
  ref.prefix           = "G:/Database/LAVA/lava-ukb-v1.1",
  phenos               = c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA")
)

library(parallel);library(data.table)
cl<-makeCluster(detectCores()-5)
clusterEvalQ(cl,pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","LAVA"))
clusterExport(cl,c("A","input","process.locus","run.univ.bivar"))

bivar_res<-parLapply(cl,seq_len(nrow(A)),function(i){
  loc.info<-A[i]
  tryCatch({
    locus<-process.locus(loc.info,input)
    b<-run.univ.bivar(locus,univ.thresh=0.05/2495)
    if(!is.null(b$bivar)&&nrow(b$bivar)>0){
      tmp<-b$bivar
      tmp$LOC<-loc.info$LOC
      tmp
    }
  },error=function(e)NULL)
})

stopCluster(cl)
bivar_all<-rbindlist(Filter(Negate(is.null),bivar_res),fill=TRUE)
gc()
write_tsv(bivar_all,file = "all_lava.txt")
