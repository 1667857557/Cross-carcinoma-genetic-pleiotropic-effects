#####EFA-------
setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")
tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
load("ldsc.covstruct_ALL.Rdata")
set.seed(100)
gc()
S_Stand <- ldsc.covstruct$S
V_Stand <- ldsc.covstruct$V
paLDSC(S=S_Stand,V=V_Stand,r=100)
ggsave('paLDSC.pdf', width= 6 , height= 6 , units='in')
S  <- ldsc.covstruct$S
corS <- cov2cor(S)
ev_cor <- eigen(corS)$values
library(nFactors)
ns <- nScree(ev_cor, cor = FALSE)
print(ns$Components)
set.seed(100)
require(Matrix)
Ssmooth<-as.matrix((nearPD(ldsc.covstruct$S, corr = FALSE))$mat)
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")
EFA$loadings
gc()
ldsc.covstruct <- ldsc(traits =c("LC.sumstats.gz","BRCA.sumstats.gz","ESCA.sumstats.gz", "EC.sumstats.gz","OC.sumstats.gz","RCC.sumstats.gz","CRC.sumstats.gz","THCA.sumstats.gz"),
                       sample.prev = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),
                       population.prev = c(0.059,0.137,0.006,0.03,0.011,0.018,0.039,0.012),
                       ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/",
                       trait.names=-c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA"))
S_Stand <- ldsc.covstruct$S
V_Stand <- ldsc.covstruct$V
paLDSC(S=S_Stand,V=V_Stand,r=100)
ggsave('paLDSC.pdf', width= 6 , height= 6 , units='in')
S  <- ldsc.covstruct$S
corS <- cov2cor(S)
ev_cor <- eigen(corS)$values
library(nFactors)
ns <- nScree(ev_cor, cor = FALSE)
print(ns$Components)
set.seed(100)
require(Matrix)
Ssmooth<-as.matrix((nearPD(ldsc.covstruct$S, corr = FALSE))$mat)
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 3, rotation = "promax")
EFA$loadings
######CF------
CommonFactor_DWLS<- commonfactor(covstruc = ldsc.covstruct, estimation="DWLS")
save(ldsc.covstruct,file = "ldsc.covstruct.Rdata")
load("ldsc.covstruct.Rdata")
CommonFactor_DWLS
setwd("D:/pan_cancer/CLEAN_ORIGINAL")
files<-c("LC_GCST004748.txt.gz","BRCA_GCST010098.txt.gz","ESCA_GCST003739.txt.gz","EC_GCST006464.txt.gz","OC_GCST004462.txt.gz","RCC_GCST90320057.txt.gz","CRC_GCST90255675.txt.gz",
         "THCA_GCST90399736.txt.gz")
ref= "D:/pan_cancer/reference.1000G.maf.0.005.txt.gz"

trait.names<-c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
gc()
se.logit=c(T,T,T,T,T,T,T,T)
linprob=c(F,F,F,F,F,F,F,F)
info.filter=.6
maf.filter=0
OLS=NULL
N=NULL
betas=NULL
CANCER_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=NULL,linprob=linprob,N=NULL,betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE,parallel=T,cores=5)
write_tsv(CANCER_sumstats,file = "CANCER_sumstats.txt.gz")
head(CANCER_sumstats)
gc()
CANCER_factor  <- commonfactorGWAS(covstruc = ldsc.covstruct,  parallel = TRUE,SNPs = CANCER_sumstats,cores = 15, smooth_check = TRUE, toler= 1e-60)

write_tsv(CANCER_factor,file = "common_factor_GWAS_LUNG_BRCA_ESC_ENDO_OVCA_RCC_CRC_THCA.txt.gz")
rm(CANCER_factor)
A<-vroom("D:/pan_cancer/CLEAN_ORIGINAL/common_factor_GWAS_LUNG_BRCA_ESC_ENDO_OVCA_RCC_CRC_THCA.txt.gz")
A <- subset(A, fail == 0 & warning == "0")
A<-subset(A, A$MAF >= .1)
A<-subset(A, A$Q_pval >= 5e-8)
N_hat_F2<-mean(1/((2*A$MAF*(1-A$MAF))*A$se_c^2))
A$N<-1/((2*A$MAF*(1-A$MAF))*A$se_c^2)
A <- dplyr::select(A, c("SNP","CHR","BP","A1","A2","est","se_c","Pval_Estimate","MAF","N"))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")
write_tsv(A,file = "common_factor.txt.gz")
######EF------
##cross relationship
model <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  
  F1 ~~ F2
  F2 ~~ F3
  F1 ~~ F3
'
library(GenomicSEM)
Anthro<-usermodel(ldsc.covstruct, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = FALSE, imp_cov = FALSE)
Anthro

####E-factor
model <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
'
library(GenomicSEM)
Anthro<-usermodel(ldsc.covstruct, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = FALSE, imp_cov = FALSE)
Anthro

model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e  ~ SNP
'
fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = CANCER_sumstats, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("e ~ SNP"),
               toler = 1e-60,  parallel = TRUE, GC="standard",MPI=FALSE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE,Q_SNP=FALSE)
saveRDS(fit1,file = "e_model.rds")
model2 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  F1 ~~ F2
  F2 ~~ F3
  F1 ~~ F3
  F1 ~ SNP
  F2 ~ SNP
  F3 ~ SNP
'
fit2<-userGWAS(covstruc = ldsc.covstruct, SNPs = CANCER_sumstats, estimation = "DWLS",
               model = model2, printwarn = TRUE, cores = 24,sub = c("F1 ~ SNP","F2 ~ SNP","F3 ~ SNP"),
               toler = 1e-60,  parallel = TRUE, GC="standard",MPI=FALSE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE,Q_SNP=TRUE)
saveRDS(fit2,file = "Multivariate_GWAS_model.rds.gz")

fit2<-readRDS("Multivariate_GWAS_model.rds.gz")
F1<-fit2[[1]]
F1<-subset(F1, F1$Q_SNP_pval >= 5e-8)
F2<-fit2[[2]]
F2<-subset(F2, F2$Q_SNP_pval >= 5e-8)
F3<-fit2[[3]]
F3<-subset(F3, F3$Q_SNP_pval >= 5e-8)
F1 <- subset(F1, error == 0 & warning == "0")
F2 <- subset(F2, error == 0 & warning == "0")
F3 <- subset(F3, error == 0 & warning == "0")


A<-vroom("e_factor_clean.txt.gz")
res1 <- A[, c("SNP","chisq","chisq_df")]       
res2 <- F1[, c("SNP","chisq","chisq_df")]       
colnames(res2) <- c("SNP","chisq_m2","df_m2")
het <- merge(res1, res2, by="SNP")
het$delta_chisq_F1 <- het$chisq   - het$chisq_m2
het$delta_df_F1    <- het$df_m2    - het$chisq_df
het$Q_SNP_pval_F1  <- pchisq(het$delta_chisq_F1,het$delta_df_F1,lower.tail = FALSE)
het<-select(het,c(SNP,delta_chisq_F1,delta_df_F1,Q_SNP_pval_F1))
A<-left_join(A,het,by = "SNP")
res2 <- F2[, c("SNP","chisq","chisq_df")]       
colnames(res2) <- c("SNP","chisq_m2","df_m2")
het <- merge(res1, res2, by="SNP")
het$delta_chisq_F2 <- het$chisq   - het$chisq_m2
het$delta_df_F2    <- het$df_m2    - het$chisq_df
het$Q_SNP_pval_F2  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(SNP,delta_chisq_F2,delta_df_F2,Q_SNP_pval_F2))
A<-left_join(A,het,by = "SNP")
res2 <- F3[, c("SNP","chisq","chisq_df")]       
colnames(res2) <- c("SNP","chisq_m2","df_m2")
het <- merge(res1, res2, by="SNP")
het$delta_chisq_F3 <- het$chisq   - het$chisq_m2
het$delta_df_F3    <- het$df_m2    - het$chisq_df
het$Q_SNP_pval_F3  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(SNP,delta_chisq_F3,delta_df_F3,Q_SNP_pval_F3))
A<-left_join(A,het,by = "SNP")
write_tsv(A,file= "origincal_e_factor.txt.gz")
gc()
A<-vroom("origincal_e_factor.txt.gz")
A <- subset(A, error == 0 & warning == "0")

A<-subset(A, A$MAF >= .1)
A<-subset(A, A$Q_SNP_pval_F1 >= 5e-8 & A$Q_SNP_pval_F2 >= 5e-8 & A$Q_SNP_pval_F3 >= 5e-8)
A<-as.data.frame(A)
library(dplyr)
A <- A %>% filter(!is.na(Q_SNP_pval_F1))
A <- A %>% filter(!is.na(Q_SNP_pval_F2))
A <- A %>% filter(!is.na(Q_SNP_pval_F3))
rm(fit2)
gc()

A$MAF <- as.numeric(as.character(A$MAF))
B<-subset(A, A$MAF <= 0.4 & A$MAF >= 0.1)
N_hat<-mean(1/((2*B$MAF*(1-B$MAF))*B$SE^2))
N_hat
A$N<-1/((2*A$MAF*(1-A$MAF))*A$SE^2)
A <- dplyr::select(A, c("SNP","CHR","BP","A1","A2","est","SE","Pval_Estimate","MAF","N"))
colnames(A)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")

write_tsv(A,file = "e_factor.txt.gz")

F1<-subset(F1, F1$MAF >= .1)
A<-subset(F1, F1$MAF <= 0.4 & F1$MAF >= 0.1)
N_hat_F1<-mean(1/((2*A$MAF*(1-A$MAF))*A$SE^2))
N_hat_F1
F1$N<-1/((2*F1$MAF*(1-F1$MAF))*F1$SE^2)
F1 <- dplyr::select(F1, c("SNP","CHR","BP","A1","A2","est","SE","Pval_Estimate","MAF","N"))
colnames(F1)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")
write_tsv(F1,file = "F1.txt.gz")
F2<-subset(F2, F2$MAF >= .1)
A<-subset(F2, F2$MAF <= 0.4 & F2$MAF >= 0.1)
N_hat_F2<-mean(1/((2*A$MAF*(1-A$MAF))*A$SE^2))
N_hat_F2

F2$N<-1/((2*F2$MAF*(1-F2$MAF))*F2$SE^2)
F2 <- dplyr::select(F2, c("SNP","CHR","BP","A1","A2","est","SE","Pval_Estimate","MAF","N"))
colnames(F2)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")
write_tsv(F2,file = "F2.txt.gz")
F3<-subset(F3, F3$MAF >= .1)
A<-subset(F3, F3$MAF <= 0.4 & F3$MAF >= 0.1)
N_hat_F3<-mean(1/((2*A$MAF*(1-A$MAF))*A$SE^2))
N_hat_F3
F3$N<-1/((2*F3$MAF*(1-F3$MAF))*F3$SE^2)
F3 <- dplyr::select(F3, c("SNP","CHR","BP","A1","A2","est","SE","Pval_Estimate","MAF","N"))
colnames(F3)<-c("SNP","CHR","POS","effect_allele","other_allele","BETA","SE","P","FRQ","N")
write_tsv(F3,file = "F3.txt.gz")
rm (list=ls ())
gc()

files<-c("common_factor.txt.gz","e_factor.txt.gz","F1.txt.gz","F2.txt.gz","F3.txt.gz")
hm3<-"eur_w_ld_chr/w_hm3.snplist"
trait.names<-c("CF","EF","F1","F2","F3")
N=c(NA,NA,NA,NA)
info.filter=0.9
maf.filter=0.01
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)
ldsc.covstruct <- ldsc(traits =c("CF.sumstats.gz","EF.sumstats.gz","F1.sumstats.gz", "F2.sumstats.gz","F3.sumstats.gz"),
                       sample.prev = c(NA,NA,NA,NA,NA),
                       population.prev = c(NA,NA,NA,NA,NA),
                       ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/",
                       trait.names=c("CF","EF","F1","F2","F3"))
