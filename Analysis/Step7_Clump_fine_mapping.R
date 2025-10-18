
######clump-------
cd G:/Database/1000G

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Database\1000G\CF.txt.gz `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Database\1000G\CF

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Database\1000G\EF.txt.gz `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Database\1000G\EF


./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Database\1000G\F1.txt.gz `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Database\1000G\F1

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Database\1000G\F2.txt.gz `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Database\1000G\F2

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Database\1000G\F3.txt.gz `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Database\1000G\F3

cd G:/Database/1000G

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_BRCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_BRCA

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_ENDO.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_ENDO


./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_OVCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_OVCA

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_CRC.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_CRC

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_ESC.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_ESC

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_LUNG.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_LUNG

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_RCC.txt `
--clump-snp-field SNP `
--clump-field P `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_RCC


./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_maxFDR_THCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_THCA

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_MVP_maxFDR_CRC.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_MVP_CRC

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_MVP_maxFDR_ESC.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_MVP_ESC

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_MVP_maxFDR_LUNG.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_MVP_LUNG

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_MVP_maxFDR_RCC.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_MVP_RCC

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_MVP_maxFDR_THCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_MVP_THCA

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_UKB_FIN_maxFDR_BRCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_UKB_FIN_BRCA

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_UKB_FIN_maxFDR_ENDO.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_UKB_FIN_ENDO

./plink2.exe --bfile G:\Database\1000G\UK10K_1KG_qc_rsid `
--clump G:\Linux\mtag\A\mtag_UKB_FIN_maxFDR_OVCA.txt `
--clump-snp-field SNP `
--clump-field mtag_pval `
--clump-p1 5e-8 `
--clump-r2 0.1 `
--clump-kb 500 `
--out G:\Linux\mtag\A\mtag_UKB_FIN_OVCA




#####perpare annotation
setwd("D:/CARMA/baselineLDv2.2")
data_files_txt <- list.files(pattern = ".annot.gz")
results_df <- data.frame()

for (file_rds in data_files_txt) {
  A<-vroom(file_rds)
  results_df<-rbind(results_df,A)
  rm(A)
  gc()
}
results_df<-dplyr::select(results_df,-c(CHR,BP))
delete_annots <- c(
  # all .flanking.500 
  grep("\\.flanking\\.500$", colnames(results_df), value = TRUE),
  # MAF related
  paste0("MAFbin", 1:10),
  c("MAF_Adj_Predicted_Allele_Age", "MAF_Adj_LLD_AFR", "MAF_Adj_ASMC"),
  # LD /Sequence composition correction
  "Recomb_Rate_10kb", "Nucleotide_Diversity_10kb",
  "Backgrd_Selection_Stat", "CpG_Content_50kb",
  c("base", "CM")
)
results_df<-dplyr::select(results_df,-delete_annots)
write_tsv(results_df,file = "baselineLDv2.2_annot_60.txt.gz")

######fine-mapping----
cd /mnt/d/CARMA
R

suppressPackageStartupMessages(pacman::p_load(dplyr, mapgen, susieR, GWAS.utils, CARMA, data.table, vroom, readr, Matrix))
gc(); tempdir <- function() "/mnt/d/rtemp1"; unlockBinding("tempdir", baseenv()); utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv()); assign("tempdir", tempdir, baseenv()); lockBinding("tempdir", baseenv())
setwd("/mnt/d/CARMA")
position1 <- "/mnt/d/CARMA/common_factor/"; plink_exec <- "/mnt/d/CARMA/plink"; dir.create(position1, showWarnings=FALSE, recursive=TRUE)

bim <- data.table::fread("/mnt/d/LD/UK10K_1KG_qc_rsid.bim")
pos_idx <- bim %>% transmute(SNP=V2, CHR=V1, POS=V4, V5=V5, V6=V6)

Aall_try <- vroom::vroom("common_factor.txt.gz") %>%
  mutate(BETA=as.numeric(BETA), SE=as.numeric(SE),
         effect_allele=toupper(effect_allele),
         other_allele=toupper(other_allele),
         allele_pair = paste0(sort(c(effect_allele, other_allele)), collapse="")) %>%
  filter(nchar(effect_allele)==1 & nchar(other_allele)==1) %>%
  filter(!allele_pair %in% c("AT","TA","CG","GC"))
abf_only <- nrow(Aall_try)==0
Aall <- if(abf_only) vroom::vroom("common_factor.txt.gz") %>% mutate(P=2*pnorm(-abs(BETA/SE)), Z=BETA/SE) else Aall_try
if(!all(c("CHR","POS") %in% names(Aall))) Aall <- Aall %>% left_join(pos_idx %>% select(SNP,CHR,POS), by="SNP")
n <- mean(Aall$N, na.rm=TRUE)
if(!abf_only) B <- vroom::vroom("./baselineLDv2.2_annot_60.txt.gz")

sig.loci <- data.table::fread("common_factor.clumps") %>% mutate(SP2=gsub("\\(\\d+\\)","",SP2), GenomicLocus=dplyr::row_number())
existing_carma <- list.files(position1, pattern="^CARMA_\\d+\\.rds\\.gz$", full.names=FALSE)
existing_abf   <- list.files(position1, pattern="^ABF_\\d+\\.txt$", full.names=FALSE)
done_loci <- c(as.integer(sub("^CARMA_(\\d+)\\.rds\\.gz$","\\1",existing_carma)), as.integer(sub("^ABF_(\\d+)\\.txt$","\\1",existing_abf)))
loci <- setdiff(sig.loci$GenomicLocus, done_loci); message("剩余要处理的 locus: ", paste(loci, collapse=", "))

safe_cs_idx <- function(res) tryCatch({cs <- res[[1]][["Credible set"]]; if(is.null(cs)) integer(0) else unique(unlist(cs, use.names=FALSE))}, error=function(e) integer(0))
make_psd <- function(R){R <- (R+t(R))/2; ev <- eigen(R, only.values=TRUE)$values; if(min(ev)<1e-8) diag(R) <- diag(R)+(1e-6-min(ev)); R}
abf_wakefield <- function(beta,se,V=0.04){ve <- se^2; z <- beta/se; sqrt(ve/(ve+V))*exp((z^2*V)/(2*(ve+V)))}
abf_write <- function(df,locus,outdir){
  df <- df %>% select(CHR,POS,Ref=other_allele,Alt=effect_allele,SNP,BETA,SE,P) %>% filter(is.finite(BETA),is.finite(SE))
  df$Z<-df$BETA/df$SE
  if(nrow(df)==0) return(FALSE)
  abf <- abf_wakefield(df$BETA,df$SE,0.04); pip <- {w <- abf; w/sum(w)}
  ord <- order(pip,decreasing=TRUE); cs <- integer(length(pip))
  if(length(pip)){k <- which(cumsum(pip[ord])>=0.95)[1]; if(!is.na(k)) cs[ord[seq_len(k)]] <- 1L}
  out <- df %>% transmute(CHR,POS,Ref,Alt,SNP,Pval=P,Z=Z,PIP=pip,CS=cs)
  outfile <- file.path(outdir,paste0("ABF_",locus,".txt")); readr::write_tsv(out,outfile)
  message("ABF written: ",outfile," (p=",nrow(df),")"); TRUE
}

for(value in loci){
  tryCatch({
    locus <- value
    row <- sig.loci %>% filter(GenomicLocus==locus) %>% slice(1); if(nrow(row)==0) next
    clump_snps <- unique(c(row$ID[1], if(is.na(row$SP2[1])||row$SP2[1]=="") character(0) else strsplit(row$SP2[1], ",")[[1]]))
    span_tbl <- pos_idx %>% filter(SNP %in% clump_snps) %>% distinct(SNP,.keep_all=TRUE)
    if(nrow(span_tbl)==0 || dplyr::n_distinct(span_tbl$CHR)!=1){message("Skip ",locus,": clump missing or multi-chr"); next}
    
    chr <- span_tbl$CHR[1]
    start <- max(1L, min(span_tbl$POS,na.rm=TRUE)-1e5)
    end_cap <- max(pos_idx$POS[pos_idx$CHR==chr], na.rm=TRUE)
    end <- min(end_cap, max(span_tbl$POS,na.rm=TRUE)+1e5)
    
    Aloc0 <- Aall %>% filter(CHR==chr, POS>=start, POS<=end) %>% distinct(SNP,.keep_all=TRUE)
    if(nrow(Aloc0)<1){message("Skip ",locus,": window empty"); next}
    
    if(abf_only){abf_write(Aloc0,locus,position1); gc(); next}
    
    Aloc <- Aloc0
    go_next <- FALSE
    repeat{
      df2 <- Aloc %>% transmute(SNP, A1)
      snpfile <- tempfile(fileext=".txt"); write.table(df2$SNP, snpfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      allelefile <- tempfile(fileext=".txt"); write.table(df2, allelefile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
      pref <- file.path(position1, paste0("LD_",locus))
      st <- system2(plink_exec, c("--bfile","/mnt/d/LD/UK10K_1KG_qc_rsid","--extract",snpfile,"--a1-allele",allelefile,"--keep-allele-order","--write-snplist","--r","square","gz","--out",pref))
      if(!file.exists(paste0(pref,".ld.gz"))){abf_write(Aloc0,locus,position1); go_next <- TRUE; break}
      
      kept <- readLines(paste0(pref,".snplist")); mat <- as.matrix(data.table::fread(paste0(pref,".ld.gz")))
      if(length(kept)!=ncol(mat) || nrow(mat)!=ncol(mat)){abf_write(Aloc0,locus,position1); go_next <- TRUE; break}
      colnames(mat) <- kept; rownames(mat) <- kept
      
      Aloc <- Aloc %>% filter(SNP %in% kept) %>% slice(match(kept,SNP))
      bad <- rowSums(is.na(mat))>0; if(any(bad)){keep <- !bad; mat <- mat[keep,keep,drop=FALSE]; Aloc <- Aloc[keep,,drop=FALSE]}
      if(nrow(Aloc)==0){abf_write(Aloc0,locus,position1); go_next <- TRUE; break}
      
      if(!isSymmetric(mat)) mat <- (mat+t(mat))/2; mat <- make_psd(mat)
      if(nrow(Aloc)<4){abf_write(Aloc0,locus,position1); go_next <- TRUE; break}
      
      dz <- LD_diagnosis_susie_rss(Aloc$Z, R=mat, n=n)
      det <- which(dz$conditional_dist$logLR>2 & abs(dz$conditional_dist$z)>2)
      if(length(det)==0) break else {Aloc <- Aloc[-det,,drop=FALSE]; mat <- mat[-det,-det,drop=FALSE]; if(nrow(Aloc)==0){abf_write(Aloc0,locus,position1); go_next <- TRUE; break}}
    }
    if(go_next) {gc(); next}
    
    if(nrow(Aloc)<4){abf_write(Aloc0,locus,position1); gc(); next}
    
    Aloc2 <- Aloc %>% transmute(CHR,POS,Ref=other_allele,Alt=effect_allele,SNP,Z,Pval=P)
    annot <- dplyr::left_join(Aloc2 %>% select(SNP), B, by="SNP")
    if(is.null(annot) || nrow(annot)!=nrow(Aloc2)){abf_write(Aloc0,locus,position1); gc(); next}
    annot[is.na(annot)] <- 0; X <- as.matrix(cbind(Intercept=1, annot %>% select(-SNP)))
    res <- tryCatch(withCallingHandlers(CARMA(list(Aloc2$Z), list(mat), w.list=list(X), lambda.list=list(1), outlier.switch=FALSE), warning=function(w){invokeRestart("muffleWarning")}), error=function(e) e)
    if(inherits(res,"error") || is.null(res) || is.null(res[[1]]$PIPs)){abf_write(Aloc0,locus,position1); gc(); next}
    
    Aloc2 <- Aloc2 %>% mutate(PIP=res[[1]]$PIPs, CS=0L); cs <- safe_cs_idx(res); if(length(cs)) Aloc2$CS[cs] <- 1L
    readr::write_tsv(Aloc2, file.path(position1, paste0("CARMA_",locus,".txt"))); saveRDS(res, file.path(position1, paste0("CARMA_",locus,".rds.gz")), compress="gzip"); gc()
  }, error=function(e){message("Skip locus ", value, ": ", conditionMessage(e))})
}
