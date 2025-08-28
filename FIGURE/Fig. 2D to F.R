######
setwd("D:/pan_cancer/CCGWAS")
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
  library(vroom); library(dplyr); library(ggplot2); library(ggrepel)
  library(ieugwasr); library(plinkbinr); library(tibble)
})

plink_bin <- "G:/plink_Windows.exe"
bfile     <- "G:/Database/1000G/UK10K_1KG_qc_rsid"
clump_kb  <- 500L
clump_r2  <- 0.1

cc   <- vroom("LUNG_PRCA.out.results.gz", show_col_types=FALSE)
LUNG <- vroom("LUNG_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_LUNG=EA, NEA_LUNG=NEA, Neff_LUNG=Neff, beta_LUNG=log(OR), p_LUNG=P)
PRCA <- vroom("PRCA_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_PRCA=EA, NEA_PRCA=NEA, Neff_PRCA=Neff, beta_PRCA=log(OR), p_PRCA=P)

need_cols <- c("SNP","CHR","BP","OLS_beta","OLS_pval")
stopifnot(all(need_cols %in% names(cc)))
exact_cands <- c("EXACT_pval","Exact_pval","P_exact","P_EXACT","CC_exact","CCEXACT_pval","EXACT_P","exact_p")
exact_col <- intersect(exact_cands, names(cc)); if(length(exact_col)==0) stop("未找到 CC-exact p 值列"); exact_col <- exact_col[1]

cc <- cc[, c("SNP","CHR","BP","OLS_beta","OLS_pval", exact_col, intersect("CCGWAS_signif", names(cc)))]
dat <- cc |> inner_join(LUNG, "SNP") |> inner_join(PRCA, "SNP")
flip0 <- with(dat, EA_LUNG==NEA_PRCA & NEA_LUNG==EA_PRCA)
same0 <- with(dat, EA_LUNG==EA_PRCA & NEA_LUNG==NEA_PRCA)
dat <- dat[same0 | flip0, ]
flip <- with(dat, EA_LUNG==NEA_PRCA & NEA_LUNG==EA_PRCA)
dat$beta_PRCA[flip] <- -dat$beta_PRCA[flip]

signedZ <- function(beta,p) sign(beta) * qnorm(pmax(p, .Machine$double.xmin)/2, lower.tail=FALSE)
dat <- dat |>
  mutate(Z_LUNG = signedZ(beta_LUNG, p_LUNG),
         Z_PRCA = signedZ(beta_PRCA, p_PRCA),
         P_exact = .data[[exact_col]],
         is_ccsig = (OLS_pval < 5e-8) & (P_exact < 1e-4),
         tag = case_when(is_ccsig & p_LUNG>=5e-8 & p_PRCA>=5e-8 ~ "CC-specific",
                         is_ccsig ~ "CC&GWAS", TRUE ~ "NS"),
         dir = ifelse(Z_LUNG*Z_PRCA<0,"Opposite","Same"),
         tilt= ifelse(OLS_beta>0,"LUNG>PRCA","PRCA>LUNG"))

trim_mean <- function(v, trim=0.1){
  qs <- quantile(v, probs=c(trim, 1-trim), na.rm=TRUE, names=FALSE)
  mean(v[v>=qs[1] & v<=qs[2]], na.rm=TRUE)
}
Neff_LUNG_ref_all <- trim_mean(dat$Neff_LUNG, 0.10)
Neff_PRCA_ref_all <- trim_mean(dat$Neff_PRCA, 0.10)
rm(cc, LUNG, PRCA); gc()
to_clump <- dat |> filter(is_ccsig, !is.na(SNP), is.finite(OLS_pval)) |>
  transmute(rsid=SNP, pval=OLS_pval)
if(nrow(to_clump)==0) stop("无满足 CC-GWAS 双阈值的 SNP，无法 clump。")

old_wd <- getwd()
tmpdir <- file.path(tempdir(), paste0("ldclump_", Sys.getpid()))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)
setwd(tmpdir)
clumped <- tryCatch(
  ieugwasr::ld_clump(
    tibble(to_clump),
    plink_bin = plink_bin,
    bfile     = bfile,
    clump_kb  = clump_kb,
    clump_r2  = clump_r2,
    clump_p   = 1
  ),
  error=function(e) { setwd(old_wd); unlink(tmpdir, recursive=TRUE, force=TRUE); stop(e) }
)
setwd(old_wd); try(unlink(tmpdir, recursive=TRUE, force=TRUE), silent=TRUE); gc()
lead_rsids <- unique(clumped$rsid)
lead <- dat |> semi_join(tibble(SNP=lead_rsids), by="SNP") |>
  filter(is_ccsig, is.finite(OLS_pval), is.finite(Z_LUNG), is.finite(Z_PRCA))
lead <- lead |>
  mutate(x = Z_LUNG / sqrt(Neff_LUNG_ref_all),
         y = Z_PRCA / sqrt(Neff_PRCA_ref_all))
zthr <- qnorm(5e-8/2, lower.tail=FALSE)
thrX <- zthr / sqrt(Neff_LUNG_ref_all)
thrY <- zthr / sqrt(Neff_PRCA_ref_all)
lead$col <- ifelse(lead$tag=="CC-specific","CC-specific","CC&GWAS")
n_ccspec <- sum(lead$col=="CC-specific", na.rm=TRUE)
n_ccgwas <- sum(lead$col=="CC&GWAS",   na.rm=TRUE)
lab_ccspec <- bquote(N["CC-specific"]==.(n_ccspec))
lab_ccgwas <- bquote(N["CC&GWAS"]==.(n_ccgwas))
p <- ggplot(lead, aes(x, y))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_hline(yintercept=c(-thrY,thrY), linetype=2)+
  geom_vline(xintercept=c(-thrX,thrX), linetype=2)+
  geom_point(aes(color=col), size=2, stroke=0, alpha=.9)+
  scale_color_manual(values=c("CC-specific"="#d62728","CC&GWAS"="#000000"),
                     breaks=c("CC-specific","CC&GWAS"),
                     labels=list(lab_ccspec, lab_ccgwas))+
  labs(title="LC vs PRCA",
       x="LC standardized effect",
       y="PRCA standardized effect")+
  theme_classic(14)+
  theme(legend.position="top", legend.title=element_blank(),
        plot.title=element_text(hjust=0.5))
write_tsv(lead,file = "PRCA_LC_CC-GWAS.txt")

setwd("D:/pan_cancer/FIGURE 3")
ggsave("PRCA_LC_CC-GWAS_scatter_clumped_globalNeff.png", p, width=5, height=5, dpi=600)
gc()
#######
setwd("D:/pan_cancer/CCGWAS")
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
  library(vroom); library(dplyr); library(ggplot2); library(ggrepel)
  library(ieugwasr); library(plinkbinr); library(tibble)
})
plink_bin <- "G:/plink_Windows.exe"
bfile     <- "G:/Database/1000G/UK10K_1KG_qc_rsid"
clump_kb  <- 500L
clump_r2  <- 0.1
cc   <- vroom("BRCA_PRCA.out.results.gz", show_col_types=FALSE)
BRCA <- vroom("BRCA_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_BRCA=EA, NEA_BRCA=NEA, Neff_BRCA=Neff, beta_BRCA=log(OR), p_BRCA=P)
PRCA <- vroom("PRCA_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_PRCA=EA, NEA_PRCA=NEA, Neff_PRCA=Neff, beta_PRCA=log(OR), p_PRCA=P)

need_cols <- c("SNP","CHR","BP","OLS_beta","OLS_pval")
stopifnot(all(need_cols %in% names(cc)))
exact_cands <- c("EXACT_pval","Exact_pval","P_exact","P_EXACT","CC_exact","CCEXACT_pval","EXACT_P","exact_p")
exact_col <- intersect(exact_cands, names(cc)); if(length(exact_col)==0) stop("未找到 CC-exact p 值列"); exact_col <- exact_col[1]

cc <- cc[, c("SNP","CHR","BP","OLS_beta","OLS_pval", exact_col, intersect("CCGWAS_signif", names(cc)))]
dat <- cc |> inner_join(BRCA, "SNP") |> inner_join(PRCA, "SNP")
flip0 <- with(dat, EA_BRCA==NEA_PRCA & NEA_BRCA==EA_PRCA)
same0 <- with(dat, EA_BRCA==EA_PRCA & NEA_BRCA==NEA_PRCA)
dat <- dat[same0 | flip0, ]
flip <- with(dat, EA_BRCA==NEA_PRCA & NEA_BRCA==EA_PRCA)
dat$beta_PRCA[flip] <- -dat$beta_PRCA[flip]

signedZ <- function(beta,p) sign(beta) * qnorm(pmax(p, .Machine$double.xmin)/2, lower.tail=FALSE)
dat <- dat |>
  mutate(Z_BRCA = signedZ(beta_BRCA, p_BRCA),
         Z_PRCA = signedZ(beta_PRCA, p_PRCA),
         P_exact = .data[[exact_col]],
         is_ccsig = (OLS_pval < 5e-8) & (P_exact < 1e-4),
         tag = case_when(is_ccsig & p_BRCA>=5e-8 & p_PRCA>=5e-8 ~ "CC-specific",
                         is_ccsig ~ "CC&GWAS", TRUE ~ "NS"),
         dir = ifelse(Z_BRCA*Z_PRCA<0,"Opposite","Same"),
         tilt= ifelse(OLS_beta>0,"BRCA>PRCA","PRCA>BRCA"))

trim_mean <- function(v, trim=0.1){
  qs <- quantile(v, probs=c(trim, 1-trim), na.rm=TRUE, names=FALSE)
  mean(v[v>=qs[1] & v<=qs[2]], na.rm=TRUE)
}
Neff_BRCA_ref_all <- trim_mean(dat$Neff_BRCA, 0.10)
Neff_PRCA_ref_all <- trim_mean(dat$Neff_PRCA, 0.10)
rm(cc, BRCA, PRCA); gc()
to_clump <- dat |> filter(is_ccsig, !is.na(SNP), is.finite(OLS_pval)) |>
  transmute(rsid=SNP, pval=OLS_pval)
if(nrow(to_clump)==0) stop("无满足 CC-GWAS 双阈值的 SNP，无法 clump。")

old_wd <- getwd()
tmpdir <- file.path(tempdir(), paste0("ldclump_", Sys.getpid()))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)
setwd(tmpdir)
clumped <- tryCatch(
  ieugwasr::ld_clump(
    tibble(to_clump),
    plink_bin = plink_bin,
    bfile     = bfile,
    clump_kb  = clump_kb,
    clump_r2  = clump_r2,
    clump_p   = 1
  ),
  error=function(e) { setwd(old_wd); unlink(tmpdir, recursive=TRUE, force=TRUE); stop(e) }
)
setwd(old_wd); try(unlink(tmpdir, recursive=TRUE, force=TRUE), silent=TRUE); gc()
lead_rsids <- unique(clumped$rsid)
lead <- dat |> semi_join(tibble(SNP=lead_rsids), by="SNP") |>
  filter(is_ccsig, is.finite(OLS_pval), is.finite(Z_BRCA), is.finite(Z_PRCA))
lead <- lead |>
  mutate(x = Z_BRCA / sqrt(Neff_BRCA_ref_all),
         y = Z_PRCA / sqrt(Neff_PRCA_ref_all))
zthr <- qnorm(5e-8/2, lower.tail=FALSE)
thrX <- zthr / sqrt(Neff_BRCA_ref_all)
thrY <- zthr / sqrt(Neff_PRCA_ref_all)
lead$col <- ifelse(lead$tag=="CC-specific","CC-specific","CC&GWAS")
n_ccspec <- sum(lead$col=="CC-specific", na.rm=TRUE)
n_ccgwas <- sum(lead$col=="CC&GWAS",   na.rm=TRUE)
lab_ccspec <- bquote(N["CC-specific"]==.(n_ccspec))
lab_ccgwas <- bquote(N["CC&GWAS"]==.(n_ccgwas))
p <- ggplot(lead, aes(x, y))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_hline(yintercept=c(-thrY,thrY), linetype=2)+
  geom_vline(xintercept=c(-thrX,thrX), linetype=2)+
  geom_point(aes(color=col), size=2, stroke=0, alpha=.9)+
  scale_color_manual(values=c("CC-specific"="#d62728","CC&GWAS"="#000000"),
                     breaks=c("CC-specific","CC&GWAS"),
                     labels=list(lab_ccspec, lab_ccgwas))+
  labs(title="BRCA vs PRCA",
       x="BRCA standardized effect",
       y="PRCA standardized effect")+
  theme_classic(14)+
  theme(legend.position="top", legend.title=element_blank(),
        plot.title=element_text(hjust=0.5))
write_tsv(lead,file = "PRCA_BRCA_CC-GWAS.txt")
setwd("D:/pan_cancer/FIGURE 3")
ggsave("PRCA_BRCA_CC-GWAS_scatter_clumped_globalNeff.png", p, width=5, height=5, dpi=600)
gc()
######
setwd("D:/pan_cancer/CCGWAS")
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
  library(vroom); library(dplyr); library(ggplot2); library(ggrepel)
  library(ieugwasr); library(plinkbinr); library(tibble)
})
plink_bin <- "G:/plink_Windows.exe"
bfile     <- "G:/Database/1000G/UK10K_1KG_qc_rsid"
clump_kb  <- 500L
clump_r2  <- 0.1
cc   <- vroom("LUNG_BRCA.out.results.gz", show_col_types=FALSE)
BRCA <- vroom("BRCA_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_BRCA=EA, NEA_BRCA=NEA, Neff_BRCA=Neff, beta_BRCA=log(OR), p_BRCA=P)
LUNG <- vroom("LUNG_CC_GWAS.txt.gz",
              col_select=c(SNP,EA,NEA,Neff,OR,P), show_col_types=FALSE) |>
  transmute(SNP, EA_LUNG=EA, NEA_LUNG=NEA, Neff_LUNG=Neff, beta_LUNG=log(OR), p_LUNG=P)

need_cols <- c("SNP","CHR","BP","OLS_beta","OLS_pval")
stopifnot(all(need_cols %in% names(cc)))
exact_cands <- c("EXACT_pval","Exact_pval","P_exact","P_EXACT","CC_exact","CCEXACT_pval","EXACT_P","exact_p")
exact_col <- intersect(exact_cands, names(cc)); if(length(exact_col)==0) stop("未找到 CC-exact p 值列"); exact_col <- exact_col[1]

cc <- cc[, c("SNP","CHR","BP","OLS_beta","OLS_pval", exact_col, intersect("CCGWAS_signif", names(cc)))]
dat <- cc |> inner_join(BRCA, "SNP") |> inner_join(LUNG, "SNP")
flip0 <- with(dat, EA_BRCA==NEA_LUNG & NEA_BRCA==EA_LUNG)
same0 <- with(dat, EA_BRCA==EA_LUNG & NEA_BRCA==NEA_LUNG)
dat <- dat[same0 | flip0, ]
flip <- with(dat, EA_BRCA==NEA_LUNG & NEA_BRCA==EA_LUNG)
dat$beta_LUNG[flip] <- -dat$beta_LUNG[flip]

signedZ <- function(beta,p) sign(beta) * qnorm(pmax(p, .Machine$double.xmin)/2, lower.tail=FALSE)
dat <- dat |>
  mutate(Z_BRCA = signedZ(beta_BRCA, p_BRCA),
         Z_LUNG = signedZ(beta_LUNG, p_LUNG),
         P_exact = .data[[exact_col]],
         is_ccsig = (OLS_pval < 5e-8) & (P_exact < 1e-4),
         tag = case_when(is_ccsig & p_BRCA>=5e-8 & p_LUNG>=5e-8 ~ "CC-specific",
                         is_ccsig ~ "CC&GWAS", TRUE ~ "NS"),
         dir = ifelse(Z_BRCA*Z_LUNG<0,"Opposite","Same"),
         tilt= ifelse(OLS_beta>0,"BRCA>LUNG","LUNG>BRCA"))

trim_mean <- function(v, trim=0.1){
  qs <- quantile(v, probs=c(trim, 1-trim), na.rm=TRUE, names=FALSE)
  mean(v[v>=qs[1] & v<=qs[2]], na.rm=TRUE)
}
Neff_BRCA_ref_all <- trim_mean(dat$Neff_BRCA, 0.10)
Neff_LUNG_ref_all <- trim_mean(dat$Neff_LUNG, 0.10)
rm(cc, BRCA, LUNG); gc()
to_clump <- dat |> filter(is_ccsig, !is.na(SNP), is.finite(OLS_pval)) |>
  transmute(rsid=SNP, pval=OLS_pval)
if(nrow(to_clump)==0) stop("无满足 CC-GWAS 双阈值的 SNP，无法 clump。")

old_wd <- getwd()
tmpdir <- file.path(tempdir(), paste0("ldclump_", Sys.getpid()))
dir.create(tmpdir, recursive=TRUE, showWarnings=FALSE)
setwd(tmpdir)
clumped <- tryCatch(
  ieugwasr::ld_clump(
    tibble(to_clump),
    plink_bin = plink_bin,
    bfile     = bfile,
    clump_kb  = clump_kb,
    clump_r2  = clump_r2,
    clump_p   = 1
  ),
  error=function(e) { setwd(old_wd); unlink(tmpdir, recursive=TRUE, force=TRUE); stop(e) }
)
setwd(old_wd); try(unlink(tmpdir, recursive=TRUE, force=TRUE), silent=TRUE); gc()
lead_rsids <- unique(clumped$rsid)
lead <- dat |> semi_join(tibble(SNP=lead_rsids), by="SNP") |>
  filter(is_ccsig, is.finite(OLS_pval), is.finite(Z_BRCA), is.finite(Z_LUNG))
lead <- lead |>
  mutate(x = Z_BRCA / sqrt(Neff_BRCA_ref_all),
         y = Z_LUNG / sqrt(Neff_LUNG_ref_all))
zthr <- qnorm(5e-8/2, lower.tail=FALSE)
thrX <- zthr / sqrt(Neff_BRCA_ref_all)
thrY <- zthr / sqrt(Neff_LUNG_ref_all)
lead$col <- ifelse(lead$tag=="CC-specific","CC-specific","CC&GWAS")
n_ccspec <- sum(lead$col=="CC-specific", na.rm=TRUE)
n_ccgwas <- sum(lead$col=="CC&GWAS",   na.rm=TRUE)
lab_ccspec <- bquote(N["CC-specific"]==.(n_ccspec))
lab_ccgwas <- bquote(N["CC&GWAS"]==.(n_ccgwas))
p <- ggplot(lead, aes(x, y))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_hline(yintercept=c(-thrY,thrY), linetype=2)+
  geom_vline(xintercept=c(-thrX,thrX), linetype=2)+
  geom_point(aes(color=col), size=2, stroke=0, alpha=.9)+
  scale_color_manual(values=c("CC-specific"="#d62728","CC&GWAS"="#000000"),
                     breaks=c("CC-specific","CC&GWAS"),
                     labels=list(lab_ccspec, lab_ccgwas))+
  labs(title="BRCA vs LC",
       x="BRCA standardized effect",
       y="LC standardized effect")+
  theme_classic(14)+
  theme(legend.position="top", legend.title=element_blank(),
        plot.title=element_text(hjust=0.5))
write_tsv(lead,file = "LC_BRCA_CC-GWAS.txt")

setwd("D:/pan_cancer/FIGURE 3")
ggsave("LC_BRCA_CC-GWAS_scatter_clumped_globalNeff.png", p, width=5, height=5, dpi=600)
gc()
######
cd /mnt/g/Linux/mixer/
  ID1=BRCA
ID2=PRCA
ID=${ID1}_${ID2}
WKDIR=/mnt/g/Linux/mixer/${ID}
MFIGURE=/mnt/g/Linux/mixer/precimed/mixer_figures.py

python ${MFIGURE} two \
--json-fit  ${WKDIR}/${ID}.fit.json \
--json-test ${WKDIR}/${ID}.test.json \
--trait1 ${ID1} --trait2 ${ID2} \
--venn-colors "#eef0ac" "#E64B35FF" \
--venn-alpha 0.95 \
--label-fontsize 12 \
--label-autoadjust \
--ext pdf png \
--out ${WKDIR}/${ID} \
--statistic mean std

cd /mnt/g/Linux/mixer/
  ID1=LUNG
ID2=PRCA
ID=${ID1}_${ID2}
WKDIR=/mnt/g/Linux/mixer/${ID}
MFIGURE=/mnt/g/Linux/mixer/precimed/mixer_figures.py

python ${MFIGURE} two \
--json-fit  ${WKDIR}/${ID}.fit.json \
--json-test ${WKDIR}/${ID}.test.json \
--trait1 LC --trait2 ${ID2} \
--venn-colors "#F39B7FFF" "#E64B35FF" \
--venn-alpha 0.95 \
--label-fontsize 12 \
--label-autoadjust \
--ext pdf png \
--out ${WKDIR}/${ID} \
--statistic mean std



cd /mnt/g/Linux/mixer/
  ID1=LUNG
ID2=BRCA
ID=${ID1}_${ID2}
WKDIR=/mnt/g/Linux/mixer/${ID}
MFIGURE=/mnt/g/Linux/mixer/precimed/mixer_figures.py

python ${MFIGURE} two \
--json-fit  ${WKDIR}/${ID}.fit.json \
--json-test ${WKDIR}/${ID}.test.json \
--trait1 LC --trait2 ${ID2} \
--venn-colors "#F39B7FFF" "#eef0ac" \
--venn-alpha 0.95 \
--label-fontsize 12 \
--label-autoadjust \
--ext pdf png \
--out ${WKDIR}/${ID} \
--statistic mean std
