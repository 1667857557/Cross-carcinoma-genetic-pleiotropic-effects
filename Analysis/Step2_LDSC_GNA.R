######LDSC-----
setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")

tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"# 修改为d盘路径
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
files<-c("LC_GCST004748.txt.gz","BRCA_GCST010098.txt.gz","ESCA_GCST003739.txt.gz","EC_GCST006464.txt.gz","OC_GCST004462.txt.gz","RCC_GCST90320057.txt.gz","CRC_GCST90255675.txt.gz",
         "THCA_GCST90399736.txt.gz","PRCA_GCST90274714.txt.gz")
hm3<-"eur_w_ld_chr/w_hm3.snplist"
trait.names<-c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA")
N=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
info.filter=0.9
maf.filter=0.01
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

ldsc.covstruct <- ldsc(traits =c("LC.sumstats.gz","BRCA.sumstats.gz","ESCA.sumstats.gz", "EC.sumstats.gz","OC.sumstats.gz","RCC.sumstats.gz","CRC.sumstats.gz","THCA.sumstats.gz","PRCA.sumstats.gz"),
                       sample.prev = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                       population.prev = c(0.059,0.137,0.006,0.03,0.011,0.018,0.039,0.012,0.127),
                       ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/",
                       trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA"))
save(ldsc.covstruct,file = "ldsc.covstruct_ALL.Rdata")

########GNA-----
setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")
require(GNA)
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"# 修改为d盘路径
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
load("ldsc.covstruct_ALL.Rdata")
covstruc<-ldsc.covstruct
fix_omega<-"full"
prune<-TRUE
p.adjust<-"fdr"
alpha<-0.05
reestimate<-TRUE
recursive<-TRUE
cancernetwork <- traitNET(covstruc=covstruc,fix_omega=fix_omega,
                          prune=prune,p.adjust=p.adjust,alpha=alpha,reestimate=reestimate,
                          recursive=recursive)
cancernetwork$model_results$sparse$parameters
