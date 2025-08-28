pacman::p_load(RadialMR,ieugwasr,vroom,TwoSampleMR,data.table,readr,tidyr,dplyr,plyr,mr.raps,MRMix,devtools,ggplot2,MendelianRandomization,future.apply)

Sys.setenv(TMPDIR="G:/rtemp",OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1)
results_df <- data.frame()
setwd("G:/MR/EUR_MR")
file_list <- list.files(pattern="*.txt.gz")

process_one <- function(file){
  exp_dat <- tryCatch(fread(file,data.table=FALSE),error=function(e) NULL); if(is.null(exp_dat) || nrow(exp_dat)==0) return(NULL)
  exp_dat$SNP <- as.character(exp_dat$SNP); exp_dat <- exp_dat[!is.na(exp_dat$SNP) & exp_dat$SNP!="",,drop=FALSE] %>% distinct(SNP,.keep_all=TRUE); exp_dat$id.exposure <- file
  outcome_raw <- tryCatch(vroom("D:/pan_cancer/F1.txt.gz",show_col_types=FALSE,progress=FALSE,num_threads=1),error=function(e) NULL); if(is.null(outcome_raw)) return(NULL)
  need <- c("SNP","BETA","SE","CHR","POS","effect_allele","other_allele","P"); if(!all(need %in% names(outcome_raw))) return(NULL)
  outcome_raw$SNP <- as.character(outcome_raw$SNP)
  outcome_raw <- outcome_raw %>% filter(!is.na(SNP) & SNP!="") %>% distinct(SNP,.keep_all=TRUE)
  outcome_raw <- dplyr::semi_join(outcome_raw,exp_dat %>% dplyr::select(SNP),by="SNP"); if(nrow(outcome_raw)==0) return(NULL)
  outcome_dat <- TwoSampleMR::format_data(outcome_raw,type="outcome",snp_col="SNP",beta_col="BETA",se_col="SE",chr_col="CHR",pos_col="POS",effect_allele_col="effect_allele",other_allele_col="other_allele",pval_col="P") %>% mutate(outcome="PAN")
  if(nrow(outcome_dat)==0) return(NULL)
  dat3 <- harmonise_data(exp_dat,outcome_dat,action=2) %>% filter(mr_keep); if(nrow(dat3)<1) return(NULL)
  
  iter <- 0
  repeat{
    dat4 <- tsmr_to_rmr_format(dat3); if(is.null(dat4) || nrow(dat4)<=3) break
    eggrad <- tryCatch(egger_radial(r_input=dat4,alpha=0.05,weights=1,summary=TRUE),error=function(e) NULL)
    if(is.null(eggrad) || is.null(eggrad$outliers) || identical(eggrad$outliers,"No significant outliers") || nrow(eggrad$outliers)==0) break
    dat3 <- dat3[!dat3$SNP %in% unique(eggrad$outliers$SNP), ]; if(nrow(dat3)<=3) break
    iter <- iter+1; if(iter>=20) break
  }
  
  if(nrow(dat3)==0) return(NULL)
  if(nrow(dat3)==1) res <- mr(dat3,method_list="mr_wald_ratio") else if(nrow(dat3)>=2) res <- mr(dat3,method_list="mr_ivw") else return(NULL)
  res[1,]
}

N <- max(1, parallel::detectCores()-5)
future::plan(future::multisession, workers=N)
res_list <- future_lapply(file_list, process_one, future.seed=TRUE)
results_df <- dplyr::bind_rows(Filter(function(x)!is.null(x) && nrow(x)>0, res_list))
write.csv(results_df,"PhWas_phenotypes_F1.csv",row.names=FALSE)
gc()
setwd("G:/MR/EUR_MR")

results_df<-fread("PhWas_phenotypes_F1.csv")
A<-subset(results_df,pval<6.06e-6)
file_list<-A$id.exposure
results_df <- data.frame()

for (i in seq_along(file_list)) {
  setwd("G:/MR/EUR_MR")
  exp_dat <- fread(file_list[i], head = TRUE, data.table = FALSE)
  exp_dat$id.exposure <- file_list[i]
  if (nrow(exp_dat) > 3) {
    setwd("D:/pan_cancer/")
    outcome_dat <- vroom("F1.txt.gz") %>%
      format_data(., type = "outcome", snps = exp_dat$SNP,  snp_col = "SNP",
                  beta_col = "BETA",  se_col = "SE",chr_col = "CHR",pos_col = "POS",
                  effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "P",samplesize_col="N", eaf_col="FRQ") %>%
      mutate(outcome = "pan")
    if (nrow(outcome_dat) > 3) {
      set.seed(1000)
      dat3 <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat, action = 2)
      dat3<-subset(dat3, mr_keep != "FALSE")
      if (nrow(dat3) > 3) {
        dat4 <- tsmr_to_rmr_format(dat3)
        if (nrow(dat4) > 3) {
          
          eggrad <- egger_radial(r_input = dat4, alpha = 0.05,weights = 1, summary = TRUE)
          if(!is.null(eggrad[["outliers"]]) && !identical(eggrad[["outliers"]], "No significant outliers")) {dat3 <- dat3[!dat3$SNP %in% eggrad$outliers$SNP, ]
          }
          
        }
        dat4 <- tsmr_to_rmr_format(dat3)
        if (nrow(dat4) > 3) {
          eggrad <- egger_radial(r_input = dat4, alpha = 0.05,weights = 1, summary = TRUE)
          if(!is.null(eggrad[["outliers"]]) && !identical(eggrad[["outliers"]], "No significant outliers")) {dat3 <- dat3[!dat3$SNP %in% eggrad$outliers$SNP, ]
          }
        }
        dat4 <- tsmr_to_rmr_format(dat3)
        if (nrow(dat4) > 3) {
          eggrad <- egger_radial(r_input = dat4, alpha = 0.05,weights = 1, summary = TRUE)
          if(!is.null(eggrad[["outliers"]]) && !identical(eggrad[["outliers"]], "No significant outliers")) {dat3 <- dat3[!dat3$SNP %in% eggrad$outliers$SNP, ]
          }
        }
        
        set.seed(1000)
        results <- mr(dat3, method_list = c("mr_ivw","mr_weighted_median","mr_simple_median"))
        select=results[1,]
        select1=results[2,]
        select2=results[3,]
        set.seed(1000)
        df <- data.frame(select,select1,select2)
        results_df <- rbind(results_df, df)
      } else {
        next
      }
    }
    else {
      next
    }
  }
  else {
    next
  }
  rm(list = setdiff(ls(), c("B","results_df","file_list")))
  gc()
  gc()
}
keep_rows <- apply(results_df[, c("b", "b.1", "b.2")], 1, function(x) {
  all(sign(x) == sign(x[1]))
})

filtered_df <- results_df[keep_rows, ]
