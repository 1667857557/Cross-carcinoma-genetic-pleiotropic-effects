###UKBPPP
setwd("G:/BLISS/results")
files<-list("UKBPPP_EUR_LC.txt.finished","UKBPPP_EUR_BRCA.txt.finished","UKBPPP_EUR_ESCA.txt.finished","UKBPPP_EUR_EC.txt.finished", "UKBPPP_EUR_OC.txt.finished", "UKBPPP_EUR_RCC.txt.finished", "UKBPPP_EUR_CRC.txt.finished","UKBPPP_EUR_THCA.txt.finished")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_PSEM_UKBPPP<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)
gc()
write_tsv(gfactor_PSEM_UKBPPP,file = "PSEM_UKBPPP.txt")
###deCODE
setwd("G:/BLISS/results")
files<-list("deCODE_LC.txt.finished","deCODE_BRCA.txt.finished","deCODE_ESCA.txt.finished","deCODE_EC.txt.finished", "deCODE_OC.txt.finished", "deCODE_RCC.txt.finished", "deCODE_CRC.txt.finished","deCODE_THCA.txt.finished")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_PSEM_deCODE<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)

write_tsv(gfactor_PSEM_deCODE,file = "PSEM_deCODE.txt")
gc()
###ARIC_EA
setwd("G:/BLISS/results")
files<-list("ARIC_EA_LC.txt.finished","ARIC_EA_BRCA.txt.finished","ARIC_EA_ESCA.txt.finished","ARIC_EA_EC.txt.finished", "ARIC_EA_OC.txt.finished", "ARIC_EA_RCC.txt.finished", "ARIC_EA_CRC.txt.finished","ARIC_EA_THCA.txt.finished")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_PSEM_ARIC_EA<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)
write_tsv(gfactor_PSEM_ARIC_EA,file = "PSEM_ARIC_EA.txt")

#######
library(GenomicSEM)
setwd("G:/BLISS/results")
files<-list("UKBPPP_EUR_LC.txt.finished","UKBPPP_EUR_BRCA.txt.finished","UKBPPP_EUR_ESCA.txt.finished","UKBPPP_EUR_EC.txt.finished", "UKBPPP_EUR_OC.txt.finished", "UKBPPP_EUR_RCC.txt.finished", "UKBPPP_EUR_CRC.txt.finished","UKBPPP_EUR_THCA.txt.finished")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  F1 ~ Gene
  F2 ~ Gene
  F3 ~ Gene
'
fit2<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("e ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
gc()
A<-fit1[[1]]
F1<-fit2[[1]]
F2<-fit2[[2]]
F3<-fit2[[3]]
res1 <- A[, c("Gene","chisq","chisq_df")]       
res2 <- F1[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F1 <- het$chisq   - het$chisq_m2
het$delta_df_F1    <- het$chisq_df - het$df_m2 
het$Q_SNP_pval_F1  <- pchisq(het$delta_chisq_F1,het$delta_df_F1,lower.tail = FALSE)
het<-select(het,c(Gene,delta_chisq_F1,delta_df_F1,Q_SNP_pval_F1))
A<-left_join(A,het,by = "Gene")
res2 <- F2[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F2 <- het$chisq   - het$chisq_m2
het$delta_df_F2    <- het$chisq_df - het$df_m2
het$Q_SNP_pval_F2  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F2,delta_df_F2,Q_SNP_pval_F2))
A<-left_join(A,het,by = "Gene")
res2 <- F3[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F3 <- het$chisq   - het$chisq_m2
het$delta_df_F3    <- het$chisq_df - het$df_m2   
het$Q_SNP_pval_F3  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F3,delta_df_F3,Q_SNP_pval_F3))
A<-left_join(A,het,by = "Gene")
gc()
write_tsv(A,file = "UKBPPP_e_factor.txt")
write_tsv(F1,file = "UKBPPP_F1.txt")
write_tsv(F2,file = "UKBPPP_F2.txt")
write_tsv(F3,file = "UKBPPP_F3.txt")
######
library(GenomicSEM)
setwd("G:/BLISS/results")
files<-list("deCODE_LC.txt.finished","deCODE_BRCA.txt.finished","deCODE_ESCA.txt.finished","deCODE_EC.txt.finished", "deCODE_OC.txt.finished", "deCODE_RCC.txt.finished", "deCODE_CRC.txt.finished","deCODE_THCA.txt.finished")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  F1 ~ Gene
  F2 ~ Gene
  F3 ~ Gene
'
fit2<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("e ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
gc()
A<-fit1[[1]]
F1<-fit2[[1]]
F2<-fit2[[2]]
F3<-fit2[[3]]
res1 <- A[, c("Gene","chisq","chisq_df")]       
res2 <- F1[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F1 <- het$chisq   - het$chisq_m2
het$delta_df_F1    <- het$chisq_df - het$df_m2 
het$Q_SNP_pval_F1  <- pchisq(het$delta_chisq_F1,het$delta_df_F1,lower.tail = FALSE)
het<-select(het,c(Gene,delta_chisq_F1,delta_df_F1,Q_SNP_pval_F1))
A<-left_join(A,het,by = "Gene")
res2 <- F2[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F2 <- het$chisq   - het$chisq_m2
het$delta_df_F2    <- het$chisq_df - het$df_m2
het$Q_SNP_pval_F2  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F2,delta_df_F2,Q_SNP_pval_F2))
A<-left_join(A,het,by = "Gene")
res2 <- F3[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F3 <- het$chisq   - het$chisq_m2
het$delta_df_F3    <- het$chisq_df - het$df_m2   
het$Q_SNP_pval_F3  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F3,delta_df_F3,Q_SNP_pval_F3))
A<-left_join(A,het,by = "Gene")

write_tsv(A,file = "deCODE_e_factor.txt")
write_tsv(F1,file = "deCODE_F1.txt")
write_tsv(F2,file = "deCODE_F2.txt")
write_tsv(F3,file = "deCODE_F3.txt")
rm (list=ls ())
gc()

###ARIC_EA
library(GenomicSEM)
setwd("G:/BLISS/results")
files<-list("ARIC_EA_LC.txt.finished","ARIC_EA_BRCA.txt.finished","ARIC_EA_ESCA.txt.finished","ARIC_EA_EC.txt.finished", "ARIC_EA_OC.txt.finished", "ARIC_EA_RCC.txt.finished", "ARIC_EA_CRC.txt.finished","ARIC_EA_THCA.txt.finished")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
setwd("G:/BLISS/results")
load("D:/pan_cancer/ldsc.covstruct.Rdata")
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  F1 ~ Gene
  F2 ~ Gene
  F3 ~ Gene
'
fit2<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 8,sub = c("e ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
gc()
A<-fit1[[1]]
F1<-fit2[[1]]
F2<-fit2[[2]]
F3<-fit2[[3]]
res1 <- A[, c("Gene","chisq","chisq_df")]       
res2 <- F1[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F1 <- het$chisq   - het$chisq_m2
het$delta_df_F1    <- het$chisq_df - het$df_m2 
het$Q_SNP_pval_F1  <- pchisq(het$delta_chisq_F1,het$delta_df_F1,lower.tail = FALSE)
het<-select(het,c(Gene,delta_chisq_F1,delta_df_F1,Q_SNP_pval_F1))
A<-left_join(A,het,by = "Gene")
res2 <- F2[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F2 <- het$chisq   - het$chisq_m2
het$delta_df_F2    <- het$chisq_df - het$df_m2
het$Q_SNP_pval_F2  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F2,delta_df_F2,Q_SNP_pval_F2))
A<-left_join(A,het,by = "Gene")
res2 <- F3[, c("Gene","chisq","chisq_df")]       
colnames(res2) <- c("Gene","chisq_m2","df_m2")
het <- merge(res1, res2, by="Gene")
het$delta_chisq_F3 <- het$chisq   - het$chisq_m2
het$delta_df_F3    <- het$chisq_df - het$df_m2   
het$Q_SNP_pval_F3  <- pchisq(het$delta_chisq, het$delta_df, lower.tail=FALSE)
het<-select(het,c(Gene,delta_chisq_F3,delta_df_F3,Q_SNP_pval_F3))
A<-left_join(A,het,by = "Gene")

write_tsv(A,file = "ARIC_EA_e_factor.txt")
write_tsv(F1,file = "ARIC_EA_F1.txt")
write_tsv(F2,file = "ARIC_EA_F2.txt")
write_tsv(F3,file = "ARIC_EA_F3.txt")
#######
A<-vroom("deCODE_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
B<-subset(A, A$Q_SNP_pval_F1 >= (0.05/2129/5) & A$Q_SNP_pval_F2 >= (0.05/2129/5) & A$Q_SNP_pval_F3 >= (0.05/2129/5))
A<-vroom("ARIC_EA_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval_F1 >= (0.05/2129/5) & A$Q_SNP_pval_F2 >= (0.05/2129/5) & A$Q_SNP_pval_F3 >= (0.05/2129/5))
B<-rbind(A,B)
B<-select(B,-c(P_adj))

A<-vroom("UKBPPP_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval_F1 >= (0.05/2129/5) & A$Q_SNP_pval_F2 >= (0.05/2129/5) & A$Q_SNP_pval_F3 >= (0.05/2129/5))
B<-rbind(A,B)
B<-subset(B, B$HSQ >= 0.01)

write_tsv(B,file = "e_factor_P_SEM.txt")

A<-vroom("deCODE_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
B<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
A<-vroom("ARIC_EA_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-select(B,-c(P_adj))

A<-vroom("UKBPPP_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-subset(B, B$HSQ >= 0.01)

write_tsv(B,file = "F1_P_SEM.txt")

A<-vroom("deCODE_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
B<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
A<-vroom("ARIC_EA_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-select(B,-c(P_adj))

A<-vroom("UKBPPP_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-subset(B, B$HSQ >= 0.01)

write_tsv(B,file = "F2_P_SEM.txt")

A<-vroom("deCODE_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
B<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
A<-vroom("ARIC_EA_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-select(B,-c(P_adj))

A<-vroom("UKBPPP_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-subset(B, B$HSQ >= 0.01)

write_tsv(B,file = "F3_P_SEM.txt")

A<-vroom("PSEM_deCODE.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
B<-subset(A, A$Q_pval >= (0.05/2129/5))
A<-vroom("PSEM_ARIC_EA.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-select(B,-c(P_adj))

A<-vroom("PSEM_UKBPPP.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/2129/5))
A<-subset(A, A$Q_pval >= (0.05/2129/5))
B<-rbind(A,B)
B<-subset(B, B$HSQ >= 0.01)

write_tsv(B,file = "CF_P_SEM.txt")
gc()
