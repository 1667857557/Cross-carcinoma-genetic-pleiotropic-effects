setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","GenomicSEM","ggplot2","tidyverse","GWAS.utils","Seurat")

tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
######EF------
load("ldsc.covstruct.Rdata")
model <- '
  F1 =~ LUNG + ESC + CRC
  F2 =~ OVCA + ENDO + BRCA
  F3 =~ RCC + THCA 
  e  =~ F1 + F2 + F3
'
library(GenomicSEM)
Anthro<-usermodel(ldsc.covstruct, estimation = "DWLS", model = model, 
                  CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Anthro
res <- Anthro$results        

λ <- subset(res, op == "=~",
            select = c(lhs, rhs, STD_Genotype, STD_Genotype_SE))

θ <- subset(res, op == "~~" & lhs == rhs,
            select = c(lhs, STD_Genotype, STD_Genotype_SE))
names(θ)[2:3] <- c("theta", "theta_se")
library(dplyr)
pie_F <- λ %>%
  filter(lhs == "e") %>%                  
  transmute(factor = rhs,
            var_e   = STD_Genotype^2) %>%   
  left_join(θ, by = c("factor"="lhs")) %>%  
  mutate(total = var_e + theta,
         prop_e   = var_e / total,         
         prop_uni = theta / total)        
library(dplyr)
lam <- res %>%
  filter(op == "=~", lhs %in% c("F1","F2","F3")) %>%
  transmute(factor = lhs,    
            trait  = rhs,
            lambda = Unstand_Est)

theta <- res %>%
  filter(op == "~~", lhs == rhs, lhs %in% lam$trait) %>%
  transmute(trait = lhs,
            theta = Unstand_Est)

pie_trait <- lam %>%
  left_join(theta, by = "trait") %>%
  mutate(var_fac  = lambda^2,
         total    = var_fac + theta,
         prop_fac = var_fac / total,
         prop_res = theta   / total) %>%
  arrange(factor, trait)

pie_trait
#####CommonFactor------
CommonFactor_DWLS<- commonfactor(covstruc = ldsc.covstruct, estimation="DWLS")
CommonFactor_DWLS
res <- CommonFactor_DWLS$results        
lam <- res %>%
  filter(op == "=~", lhs %in% c("F1")) %>%
  transmute(factor = lhs,     
            trait  = rhs,
            lambda = Unstandardized_Estimate)
theta <- res %>%
  filter(op == "~~", lhs == rhs, lhs %in% lam$trait) %>%
  transmute(trait = lhs,
            theta = Unstandardized_Estimate)
pie_trait <- lam %>%
  left_join(theta, by = "trait") %>%
  mutate(var_fac  = lambda^2,
         total    = var_fac + theta,
         prop_fac = var_fac / total,
         prop_res = theta   / total) %>%
  arrange(factor, trait)

pie_trait
