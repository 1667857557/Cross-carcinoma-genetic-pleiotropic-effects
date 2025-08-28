pacman::p_load(vroom, data.table, dplyr, readr, tidyr, ggplot2, tidyverse,
               GWAS.utils, Seurat)

setwd("D:/pan_cancer")
dt <- vroom("RCC_GCST90320057.txt.gz")
FRQ_ref <- vroom("G:/Database/1000G/UK10K_1KG_freq.afreq",
                 col_select = c("ID", "ALT", "ALT_FREQS")) |>
  transmute(SNP = ID, ALT, ALT_FREQS)
dt <- dt |>
  left_join(FRQ_ref, by = "SNP") |>
  mutate(FRQ = if_else(effect_allele == ALT, ALT_FREQS, 1 - ALT_FREQS)) |>
  drop_na(FRQ) |>
  rename(EA   = effect_allele,
         NEA  = other_allele,
         Neff = N,
         BP   = POS) |>
  mutate(OR = exp(BETA)) |>
  select(SNP, CHR, BP, EA, NEA, FRQ, OR, SE, P, Neff)
write_tsv(dt, file = "RCC_CC_GWAS.txt.gz")

rm(dt)
gc()

pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils","Seurat")
setwd("D:/pan_cancer")

library(R.utils)
library(CCGWAS)
library(ldscr)
m = (6096+4530)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_PRCA.out" , A_name = "LUNG" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()


m = (7599+4530)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_PRCA.out" , A_name = "BRCA" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()

m = (1484+4530)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_PRCA.out" , A_name = "CRC" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()

m = (1052+4530)/2
A<-vroom("ENDO.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ENDO_PRCA.out" , A_name = "ENDO" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./ENDO_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.03 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 12906 , N_B1 = 108979 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()
m = (1015+4530)/2
A<-vroom("OVCA.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "OVCA_PRCA.out" , A_name = "OVCA" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./OVCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.011 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 16924 , N_B1 = 68502 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()

m = (2220+4530)/2
A<-vroom("RCC.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "RCC_PRCA.out" , A_name = "RCC" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./RCC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.018 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 25890 , N_B1 = 743585 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()
m = (3641+4530)/2
A<-vroom("ESC.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ESC_PRCA.out" , A_name = "ESC" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./ESC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.006 , K_A1A0_high = 0.01 , K_A1A0_low = 0.001 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 4112 , N_B1 = 17159 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()

m = (1500+4530)/2
A<-vroom("THCA.sumstats.gz")
B<-vroom("PRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "THCA_PRCA.out" , A_name = "THCA" , B_name = "PRCA" , 
        sumstats_fileA1A0 = "./THCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./PRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.012 , K_A1A0_high = 0.05 , K_A1A0_low = 0.001 ,  
        K_B1B0 = 0.127 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 6015 , N_B1 = 1333754 , N_A0 = 56450 , N_B0 = 604640 , N_overlap_A0B0 = 0)

gc()

m = (6096+7599)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("BRCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_BRCA.out" , A_name = "LUNG" , B_name = "BRCA" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./BRCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.137 , K_B1B0_high = 0.2 , K_B1B0_low = 0.1 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 133384 , N_B0 = 113789 , N_overlap_A0B0 = 0)

gc()

m = (6096+1484)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("CRC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_CRC.out" , A_name = "LUNG" , B_name = "CRC" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./CRC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.039 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 78473 , N_B0 = 28670 , N_overlap_A0B0 = 0)

gc()
m = (6096+1052)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("ENDO.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_ENDO.out" , A_name = "LUNG" , B_name = "ENDO" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ENDO_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.03 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 12906 , N_B0 = 108979 , N_overlap_A0B0 = 0)

gc()

m = (6096+1015)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("OVCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_OVCA.out" , A_name = "LUNG" , B_name = "OVCA" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./OVCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.011 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 16924 , N_B0 = 68502 , N_overlap_A0B0 = 0)

gc()
m = (6096+2220)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("RCC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_RCC.out" , A_name = "LUNG" , B_name = "RCC" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./RCC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.018 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 25890 , N_B0 = 743585 , N_overlap_A0B0 = 0)

gc()
m = (6096+3641)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_ESC.out" , A_name = "LUNG" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.05 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (6096+1500)/2
A<-vroom("LUNG.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "LUNG_THCA.out" , A_name = "LUNG" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./LUNG_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.059 , K_A1A0_high = 0.1 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 29266 , N_B1 = 122188 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()

m = (7599+1484)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("CRC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_CRC.out" , A_name = "BRCA" , B_name = "CRC" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./CRC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.039 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 78473 , N_B0 = 28670 , N_overlap_A0B0 = 0)

gc()
m = (7599+1052)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("ENDO.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_ENDO.out" , A_name = "BRCA" , B_name = "ENDO" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ENDO_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.03 , K_B1B0_high = 0.05 , K_B1B0_low = 0.011 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 12906 , N_B0 = 108979 , N_overlap_A0B0 = 0)

gc()

m = (7599+1015)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("OVCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_OVCA.out" , A_name = "BRCA" , B_name = "OVCA" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./OVCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.011 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 16924 , N_B0 = 68502 , N_overlap_A0B0 = 0)

gc()
m = (7599+2220)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("RCC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_RCC.out" , A_name = "BRCA" , B_name = "RCC" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./RCC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.018 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 25890 , N_B0 = 743585 , N_overlap_A0B0 = 0)

gc()
m = (7599+3641)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_ESC.out" , A_name = "BRCA" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.01 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (7599+1500)/2
A<-vroom("BRCA.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "BRCA_THCA.out" , A_name = "BRCA" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./BRCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.137 , K_A1A0_high = 0.2 , K_A1A0_low = 0.1 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 133384 , N_B1 = 113789 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()




m = (1484+1052)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("ENDO.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_ENDO.out" , A_name = "CRC" , B_name = "ENDO" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ENDO_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,   
        K_B1B0 = 0.03 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 12906 , N_B0 = 108979 , N_overlap_A0B0 = 0)

gc()

m = (1484+1015)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("OVCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_OVCA.out" , A_name = "CRC" , B_name = "OVCA" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./OVCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 , 
        K_B1B0 = 0.011 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 16924 , N_B0 = 68502 , N_overlap_A0B0 = 0)

gc()
m = (1484+2220)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("RCC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_RCC.out" , A_name = "CRC" , B_name = "RCC" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./RCC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 , 
        K_B1B0 = 0.018 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 25890 , N_B0 = 743585 , N_overlap_A0B0 = 0)

gc()
m = (1484+3641)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_ESC.out" , A_name = "CRC" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.01 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (1484+1500)/2
A<-vroom("CRC.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "CRC_THCA.out" , A_name = "CRC" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./CRC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.039 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 78473 , N_B1 = 28670 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()

m = (1052+1015)/2
A<-vroom("ENDO.sumstats.gz")
B<-vroom("OVCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ENDO_OVCA.out" , A_name = "ENDO" , B_name = "OVCA" , 
        sumstats_fileA1A0 = "./ENDO_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./OVCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.03 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 , 
        K_B1B0 = 0.011 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 12906 , N_B1 = 108979 , N_A0 = 16924 , N_B0 = 68502 , N_overlap_A0B0 = 0)

gc()
m = (1052+2220)/2
A<-vroom("ENDO.sumstats.gz")
B<-vroom("RCC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ENDO_RCC.out" , A_name = "ENDO" , B_name = "RCC" , 
        sumstats_fileA1A0 = "./ENDO_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./RCC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.03 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 , 
        K_B1B0 = 0.018 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 12906 , N_B1 = 108979 , N_A0 = 25890 , N_B0 = 743585 , N_overlap_A0B0 = 0)

gc()
m = (1052+3641)/2
A<-vroom("ENDO.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ENDO_ESC.out" , A_name = "ENDO" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./ENDO_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.03 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.01 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 12906 , N_B1 = 108979 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (1052+1500)/2
A<-vroom("ENDO.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ENDO_THCA.out" , A_name = "ENDO" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./ENDO_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.03 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 12906 , N_B1 = 108979 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()
gc()
m = (1015+2220)/2
A<-vroom("OVCA.sumstats.gz")
B<-vroom("RCC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "OVCA_RCC.out" , A_name = "OVCA" , B_name = "RCC" , 
        sumstats_fileA1A0 = "./OVCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./RCC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.011 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 , 
        K_B1B0 = 0.018 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 16924 , N_B1 = 68502 , N_A0 = 25890 , N_B0 = 743585 , N_overlap_A0B0 = 0)

gc()
m = (1015+3641)/2
A<-vroom("OVCA.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "OVCA_ESC.out" , A_name = "OVCA" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./OVCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.011 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.01 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 16924 , N_B1 = 68502 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (1015+1500)/2
A<-vroom("OVCA.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "OVCA_THCA.out" , A_name = "OVCA" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./OVCA_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.011 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 16924 , N_B1 = 68502 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()
gc()
m = (2220+3641)/2
A<-vroom("RCC.sumstats.gz")
B<-vroom("ESC.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "RCC_ESC.out" , A_name = "RCC" , B_name = "ESC" , 
        sumstats_fileA1A0 = "./RCC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./ESC_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.018 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.006 , K_B1B0_high = 0.01 , K_B1B0_low = 0.001 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 25890 , N_B1 = 743585 , N_A0 = 4112 , N_B0 = 17159 , N_overlap_A0B0 = 0)

gc()
m = (2220+1500)/2
A<-vroom("RCC.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "RCC_THCA.out" , A_name = "RCC" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./RCC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.018 , K_A1A0_high = 0.05 , K_A1A0_low = 0.01 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 25890 , N_B1 = 743585 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()
m = (3641+1500)/2
A<-vroom("ESC.sumstats.gz")
B<-vroom("THCA.sumstats.gz")

rg_res <- ldsc_rg(munged_sumstats = list(
  "A" =A,"B" = B),ld = "./eur_w_ld_chr/",wld = "./eur_w_ld_chr/")
rm(A,B)
gc()
CCGWAS( outcome_file = "ESC_THCA.out" , A_name = "ESC" , B_name = "THCA" , 
        sumstats_fileA1A0 = "./ESC_CC_GWAS.txt.gz" ,
        sumstats_fileB1B0 = "./THCA_CC_GWAS.txt.gz" ,
        K_A1A0 = 0.006 , K_A1A0_high = 0.01 , K_A1A0_low = 0.001 ,  
        K_B1B0 = 0.012 , K_B1B0_high = 0.05 , K_B1B0_low = 0.01 ,intercept_A1A0_A1B1 =rg_res$h2$intercept[1],intercept_B1B0_A1B1=rg_res$h2$intercept[2],
        h2l_A1A0 = rg_res$h2$h2_observed[1] , h2l_B1B0 = rg_res$h2$h2_observed[2] , rg_A1A0_B1B0 =rg_res$rg$rg,
        intercept_A1A0_B1B0 = rg_res[["raw"]][["I"]][1,2] , 
        m = m ,  subtype_data = FALSE,
        N_A1 = 4112 , N_B1 = 17159 , N_A0 = 6015 , N_B0 = 1333754 , N_overlap_A0B0 = 0)

gc()
