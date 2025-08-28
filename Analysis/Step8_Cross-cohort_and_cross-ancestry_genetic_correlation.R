#####genetic correlation-----
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")
tempdir<-function() "D:\\rtemp"
unlockBinding("tempdir",baseenv());assignInNamespace("tempdir",tempdir,ns="base",envir=baseenv());lockBinding("tempdir",baseenv())
ldsc.covstruct <- ldsc(traits =c("EF.sumstats.gz","F1.sumstats.gz","F2.sumstats.gz","F3.sumstats.gz","CF.sumstats.gz","LC.sumstats.gz","BRCA.sumstats.gz","ESCA.sumstats.gz", "EC.sumstats.gz","OC.sumstats.gz","RCC.sumstats.gz","CRC.sumstats.gz","THCA.sumstats.gz","PRCA.sumstats.gz"),
                       sample.prev = c(NA,NA,NA,NA,NA,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                       population.prev = c(NA,NA,NA,NA,NA,0.059,0.137,0.006,0.03,0.011,0.018,0.039,0.012,0.127),
                       ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/",
                       trait.names=c("EF","F1","F2","F3","CF","LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA"))

files<-c("meta_analysis_ukbb_summary_stats_finngen_R12_C3_CANCER_EXALLC_meta_out.txt.gz","meta_analysis_ukbb_summary_stats_finngen_R12_C3_STOMACH_EXALLC_meta_out.txt.gz", "meta_analysis_ukbb_summary_stats_finngen_R12_C3_PANCREAS_EXALLC_meta_out.txt.gz", 
         "meta_analysis_ukbb_summary_stats_finngen_R12_C3_HEPATOCELLU_CARC_EXALLC_meta_out.txt.gz","meta_analysis_ukbb_summary_stats_finngen_R12_C3_HEAD_AND_NECK_EXALLC_meta_out.txt.gz",
         "meta_analysis_ukbb_summary_stats_finngen_R12_C3_CERVIX_UTERI_EXALLC_meta_out.txt.gz","meta_analysis_ukbb_summary_stats_finngen_R12_C3_BLADDER_EXALLC_meta_out.txt.gz",
         "meta_analysis_ukbb_summary_stats_finngen_R12_C3_CORPUS_UTERI_EXALLC_meta_out.txt.gz","meta_analysis_ukbb_summary_stats_finngen_R12_C3_THYMUS_EXALLC_meta_out.txt.gz")

hm3<-"eur_w_ld_chr/w_hm3.snplist"

trait.names<-c("All_Cancer_ukbb_finngen","STAD_ukbb_finngen","PAAD_ukbb_finngen","LIHC_ukbb_finngen","HNSC_ukbb_finngen","CESC_ukbb_finngen","BLCA_ukbb_finngen","UCEC_ukbb_finngen","THYM_ukbb_finngen")

N=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)

info.filter=0.9
maf.filter=0.01
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

ldsc.covstruct <- ldsc(traits =c("EF.sumstats.gz","F1.sumstats.gz","F2.sumstats.gz","F3.sumstats.gz","CF.sumstats.gz","All_Cancer_ukbb_finngen.sumstats.gz","CESC_ukbb_finngen.sumstats.gz", "BLCA_ukbb_finngen.sumstats.gz","HNSC_ukbb_finngen.sumstats.gz","PAAD_ukbb_finngen.sumstats.gz","STAD_ukbb_finngen.sumstats.gz","LIHC_ukbb_finngen.sumstats.gz"),
                       sample.prev = c(NA,NA,NA,NA,NA,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                       population.prev = c(NA,NA,NA,NA,NA,0.403,0.006,0.025,0.014,0.016,0.006,0.009),
                       ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/",
                       trait.names = c("EF","F1","F2","F3","CF","ALL_CANCER","CESC","BLCA","HNSC","PAAD","STAD","LIHC"))

######EAS_process-----
setwd("G:/")
A<-vroom("allbc_eas_hg19.txt.gz")
colnames(A)
A<-dplyr::select(A,c(SNPID,CHR,POS,EA,NEA,BETA,SE,P,EAF))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF")
gc()
Ncases<-21319 
Ncontrol<-106116
Ntotal<-Ncases+Ncontrol
v<-Ncases/(Ntotal)
TotalNeff<-4*v*(1-v)*(Ntotal)
A$Neff <- 4 / ((2 * A$EAF * (1 - A$EAF)) * (A$SE^2))
A$Neff<-ifelse(A$Neff > 1.1*TotalNeff, 1.1*TotalNeff, A$Neff)
A$Neff<-ifelse(A$Neff < 0.5*TotalNeff, 0.5*TotalNeff, A$Neff)
A$Neff <- round(A$Neff, 2)
A$Neff<-as.numeric(A$Neff)
A<-subset(A,SE>0.00001)
A <- subset(A, abs(BETA) <= 3)
A$MAF<-eaf2maf(A$EAF)
A<-subset(A,MAF>0.01)
A<-na.omit(A)
A<-dplyr::select(A,c(SNP,CHR,POS,effect_allele,other_allele,BETA,SE,P,EAF,Neff))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N")
write_tsv(A,file = "BRCA_GCST90429846.txt.gz")


setwd("G:/")
A<-vroom("PRCA_GCST90274716.txt.gz")
colnames(A)
A<-dplyr::select(A,c(SNP,CHR,BP,A2,A1,BETA,SE,P,FRQ))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF")
gc()
Ncases<-10809 
Ncontrol<-95790
Ntotal<-Ncases+Ncontrol
v<-Ncases/(Ntotal)
TotalNeff<-4*v*(1-v)*(Ntotal)
A$Neff <- 4 / ((2 * A$EAF * (1 - A$EAF)) * (A$SE^2))
A$Neff<-ifelse(A$Neff > 1.1*TotalNeff, 1.1*TotalNeff, A$Neff)
A$Neff<-ifelse(A$Neff < 0.5*TotalNeff, 0.5*TotalNeff, A$Neff)
A$Neff <- round(A$Neff, 2)
A$Neff<-as.numeric(A$Neff)
A<-subset(A,SE>0.00001)
A <- subset(A, abs(BETA) <= 3)
A$MAF<-eaf2maf(A$EAF)
A<-subset(A,MAF>0.01)
A<-na.omit(A)
A<-dplyr::select(A,c(SNP,CHR,POS,effect_allele,other_allele,BETA,SE,P,EAF,Neff))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N")
write_tsv(A,file = "PRCA_GCST90274716.txt.gz")

setwd("G:/")
A<-vroom("STAD_GCST90018629.txt.gz")
colnames(A)
A<-dplyr::select(A,c(SNP,CHR,BP,A2,A1,BETA,SE,P,FRQ))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF")
gc()
Ncases<-7921 
Ncontrol<-159201
Ntotal<-Ncases+Ncontrol
v<-Ncases/(Ntotal)
TotalNeff<-4*v*(1-v)*(Ntotal)
A$Neff <- 4 / ((2 * A$EAF * (1 - A$EAF)) * (A$SE^2))
A$Neff<-ifelse(A$Neff > 1.1*TotalNeff, 1.1*TotalNeff, A$Neff)
A$Neff<-ifelse(A$Neff < 0.5*TotalNeff, 0.5*TotalNeff, A$Neff)
A$Neff <- round(A$Neff, 2)
A$Neff<-as.numeric(A$Neff)
A<-subset(A,SE>0.00001)
A <- subset(A, abs(BETA) <= 3)
A$MAF<-eaf2maf(A$EAF)
A<-subset(A,MAF>0.01)
A<-na.omit(A)
A<-dplyr::select(A,c(SNP,CHR,POS,effect_allele,other_allele,BETA,SE,P,EAF,Neff))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N")
write_tsv(A,file = "STAD_GCST90018629.txt.gz")

setwd("G:/")
A<-vroom("CRC_GCST005591.txt.gz")
colnames(A)
A<-dplyr::select(A,c(SNP,CHR,BP,A2,A1,BETA,SE,P,FRQ))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF")
gc()
Ncases<-23572 
Ncontrol<-48700
Ntotal<-Ncases+Ncontrol
v<-Ncases/(Ntotal)
TotalNeff<-4*v*(1-v)*(Ntotal)
A$Neff <- 4 / ((2 * A$EAF * (1 - A$EAF)) * (A$SE^2))
A$Neff<-ifelse(A$Neff > 1.1*TotalNeff, 1.1*TotalNeff, A$Neff)
A$Neff<-ifelse(A$Neff < 0.5*TotalNeff, 0.5*TotalNeff, A$Neff)
A$Neff <- round(A$Neff, 2)
A$Neff<-as.numeric(A$Neff)
A<-subset(A,SE>0.00001)
A <- subset(A, abs(BETA) <= 3)
A$MAF<-eaf2maf(A$EAF)
A<-subset(A,MAF>0.01)
A<-na.omit(A)
A<-dplyr::select(A,c(SNP,CHR,POS,effect_allele,other_allele,BETA,SE,P,EAF,Neff))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","EAF","N")
write_tsv(A,file = "CRC_GCST005591.txt.gz")


setwd("D:/pan_cancer")
files<-c("STAD_GCST90018629.txt.gz","BRCA_GCST90429846.txt.gz","CRC_GCST005591.txt.gz","PRCA_GCST90274716.txt.gz")
hm3<-"eur_w_ld_chr/w_hm3.snplist"
trait.names<-c("STAD_EAS","BRCA_EAS","CRC_EAS","PRCA_EAS")
N=c(NA,NA,NA,NA)
info.filter=0.9
maf.filter=0.01
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

ldsc.covstruct <- ldsc(traits =c("BRCA_EAS.sumstats.gz","STAD_EAS.sumstats.gz","CRC_EAS.sumstats.gz","PRCA_EAS.sumstats.gz"),
                       sample.prev = c(0.5,0.5,0.5,0.5),
                       population.prev = c(0.125,0.006,0.038,0.092),
                       ld = "./eas_w_ld_chr/",wld = "./eas_w_ld_chr/",
                       trait.names=c("BRCA_EAS","STAD_EAS","CRC_EAS","PRCA_EAS"))

set.seed(1000)

S  <- ldsc.covstruct$S
corS <- cov2cor(S)
ev_cor <- eigen(corS)$values
library(nFactors)
ns <- nScree(ev_cor, cor = FALSE)
print(ns$Components)
save(ldsc.covstruct,file = "eas_ldsc.covstruct.Rdata")
require(Matrix)
Ssmooth<-as.matrix((nearPD(ldsc.covstruct$S, corr = FALSE))$mat)
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 1, rotation = "promax")
EFA$loadings
gc()
CommonFactor_DWLS<- commonfactor(covstruc = ldsc.covstruct, estimation="DWLS")
CommonFactor_DWLS$results
trait.names=c("BRCA","STAD","CRC","PRCA"),)
A<-vroom("reference.1000G.maf_eas.0.005.txt.gz")
colnames(A)<-c("SNP","CHR","BP","MAF","A1","A2")
write_tsv(A,file = "reference.1000G.maf_eas.0.005.txt.gz")
head(A)
files<-c("BRCA_GCST90429846.txt.gz","STAD_GCST90018629.txt.gz","CRC_GCST005591.txt.gz","PRCA_GCST90274716.txt.gz")
ref= "D:/pan_cancer/reference.1000G.maf_eas.0.005.txt.gz"

trait.names<-c("BRCA","STAD","CRC","PRCA")
gc()
se.logit=c(T,T,T,T)
linprob=c(F,F,F,F)
info.filter=.6
maf.filter=0
OLS=NULL
N=NULL
betas=NULL

CANCER_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=NULL,linprob=linprob,N=NULL,betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE,parallel=T,cores=5)
write_tsv(CANCER_sumstats,file = "EAS_CANCER_sumstats.txt.gz")
CANCER_sumstats<-vroom("EAS_CANCER_sumstats.txt.gz")
load("eas_ldsc.covstruct.Rdata")
gc()
ncores <- detectCores() - 5
EAS_factor <- commonfactorGWAS(covstruc = ldsc.covstruct, SNPs = CANCER_sumstats, estimation = "DWLS", cores = ncores ,toler = 1e-60,SNPSE = FALSE, parallel = TRUE, GC="standard",MPI=FALSE)
saveRDS(EAS_factor,file = "commonfactor_eas.rds.gz")

A<-readRDS("commonfactor_eas.rds.gz")
A<-subset(A, A$MAF >= .1)
A<-subset(A, A$Q_pval >= 5e-8)
N_hat_F2<-mean(1/((2*A$MAF*(1-A$MAF))*A$se_c^2))
A$N<-1/((2*A$MAF*(1-A$MAF))*A$se_c^2)
A <- dplyr::select(A, c("SNP","CHR","BP","A1","A2","est","se_c","Pval_Estimate","MAF","N"))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")
write_tsv(A,file = "common_factor_EAS_GWAS.txt.gz")
A<-vroom("./common_factor.txt.gz")
A<-dplyr::select(A,c(SNP,effect_allele,other_allele,BETA,SE,FRQ,N))
colnames(A)<-c("rsid","a2","a1","beta","SE","af","N")
gc()
write_tsv(A,file = "popcorn_common_factor_EUR.txt")     

#####POPCORN-----
cd /mnt/g/linux/Popcorn
source activate popcorn3
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --use_mle --sfile1 /mnt/d/pan_cancer/F1.sumstats.gz --sfile2 /mnt/d/pan_cancer/common_factor_EAS.sumstats.gz F1_EUR_EAS.txt
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --use_mle --sfile1 /mnt/d/pan_cancer/F2.sumstats.gz --sfile2 /mnt/d/pan_cancer/common_factor_EAS.sumstats.gz F2_EUR_EAS.txt
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --use_mle --sfile1 /mnt/d/pan_cancer/F3.sumstats.gz --sfile2 /mnt/d/pan_cancer/common_factor_EAS.sumstats.gz F3_EUR_EAS.txt
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --use_mle --sfile1 /mnt/d/pan_cancer/CF.sumstats.gz --sfile2 /mnt/d/pan_cancer/common_factor_EAS.sumstats.gz common_factor_EUR_EAS.txt
popcorn fit -v 1 --cfile EUR_EAS_all_gen_imp.cscore --use_mle --sfile1 /mnt/d/pan_cancer/EF.sumstats.gz --sfile2 /mnt/d/pan_cancer/common_factor_EAS.sumstats.gz e_factor_EUR_EAS.txt
