cd /mnt/d/linux/fusion

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/BRCA_GCST010098.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out BRCA_GCST010098_sCCA1${CHR}.dat
done
> BRCA_GCST010098_sCCA1_TWAS.txt
cp BRCA_GCST010098_sCCA11.dat BRCA_GCST010098_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 BRCA_GCST010098_sCCA1${i}.dat >> BRCA_GCST010098_sCCA1_TWAS.txt
done
rm BRCA_GCST010098_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/BRCA_GCST010098.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out BRCA_GCST010098_sCCA2${CHR}.dat
done
> BRCA_GCST010098_sCCA2_TWAS.txt
cp BRCA_GCST010098_sCCA21.dat BRCA_GCST010098_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 BRCA_GCST010098_sCCA2${i}.dat >> BRCA_GCST010098_sCCA2_TWAS.txt
done
rm BRCA_GCST010098_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/BRCA_GCST010098.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out BRCA_GCST010098_sCCA3${CHR}.dat
done
> BRCA_GCST010098_sCCA3_TWAS.txt
cp BRCA_GCST010098_sCCA31.dat BRCA_GCST010098_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 BRCA_GCST010098_sCCA3${i}.dat >> BRCA_GCST010098_sCCA3_TWAS.txt
done
rm BRCA_GCST010098_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/CRC_GCST90255675.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out CRC_GCST90255675_sCCA1${CHR}.dat
done
> CRC_GCST90255675_sCCA1_TWAS.txt
cp CRC_GCST90255675_sCCA11.dat CRC_GCST90255675_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 CRC_GCST90255675_sCCA1${i}.dat >> CRC_GCST90255675_sCCA1_TWAS.txt
done
rm CRC_GCST90255675_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/CRC_GCST90255675.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out CRC_GCST90255675_sCCA2${CHR}.dat
done
> CRC_GCST90255675_sCCA2_TWAS.txt
cp CRC_GCST90255675_sCCA21.dat CRC_GCST90255675_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 CRC_GCST90255675_sCCA2${i}.dat >> CRC_GCST90255675_sCCA2_TWAS.txt
done
rm CRC_GCST90255675_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/CRC_GCST90255675.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out CRC_GCST90255675_sCCA3${CHR}.dat
done
> CRC_GCST90255675_sCCA3_TWAS.txt
cp CRC_GCST90255675_sCCA31.dat CRC_GCST90255675_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 CRC_GCST90255675_sCCA3${i}.dat >> CRC_GCST90255675_sCCA3_TWAS.txt
done
rm CRC_GCST90255675_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/EC_GCST006464.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out EC_GCST006464_sCCA1${CHR}.dat
done
> EC_GCST006464_sCCA1_TWAS.txt
cp EC_GCST006464_sCCA11.dat EC_GCST006464_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 EC_GCST006464_sCCA1${i}.dat >> EC_GCST006464_sCCA1_TWAS.txt
done
rm EC_GCST006464_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/EC_GCST006464.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out EC_GCST006464_sCCA2${CHR}.dat
done
> EC_GCST006464_sCCA2_TWAS.txt
cp EC_GCST006464_sCCA21.dat EC_GCST006464_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 EC_GCST006464_sCCA2${i}.dat >> EC_GCST006464_sCCA2_TWAS.txt
done
rm EC_GCST006464_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/EC_GCST006464.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out EC_GCST006464_sCCA3${CHR}.dat
done
> EC_GCST006464_sCCA3_TWAS.txt
cp EC_GCST006464_sCCA31.dat EC_GCST006464_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 EC_GCST006464_sCCA3${i}.dat >> EC_GCST006464_sCCA3_TWAS.txt
done
rm EC_GCST006464_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/ESCA_GCST003739.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out ESCA_GCST003739_sCCA1${CHR}.dat
done
> ESCA_GCST003739_sCCA1_TWAS.txt
cp ESCA_GCST003739_sCCA11.dat ESCA_GCST003739_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 ESCA_GCST003739_sCCA1${i}.dat >> ESCA_GCST003739_sCCA1_TWAS.txt
done
rm ESCA_GCST003739_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/ESCA_GCST003739.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out ESCA_GCST003739_sCCA2${CHR}.dat
done
> ESCA_GCST003739_sCCA2_TWAS.txt
cp ESCA_GCST003739_sCCA21.dat ESCA_GCST003739_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 ESCA_GCST003739_sCCA2${i}.dat >> ESCA_GCST003739_sCCA2_TWAS.txt
done
rm ESCA_GCST003739_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/ESCA_GCST003739.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out ESCA_GCST003739_sCCA3${CHR}.dat
done
> ESCA_GCST003739_sCCA3_TWAS.txt
cp ESCA_GCST003739_sCCA31.dat ESCA_GCST003739_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 ESCA_GCST003739_sCCA3${i}.dat >> ESCA_GCST003739_sCCA3_TWAS.txt
done
rm ESCA_GCST003739_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/LC_GCST004748.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out LC_GCST004748_sCCA1${CHR}.dat
done
> LC_GCST004748_sCCA1_TWAS.txt
cp LC_GCST004748_sCCA11.dat LC_GCST004748_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 LC_GCST004748_sCCA1${i}.dat >> LC_GCST004748_sCCA1_TWAS.txt
done
rm LC_GCST004748_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/LC_GCST004748.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out LC_GCST004748_sCCA2${CHR}.dat
done
> LC_GCST004748_sCCA2_TWAS.txt
cp LC_GCST004748_sCCA21.dat LC_GCST004748_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 LC_GCST004748_sCCA2${i}.dat >> LC_GCST004748_sCCA2_TWAS.txt
done
rm LC_GCST004748_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/LC_GCST004748.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out LC_GCST004748_sCCA3${CHR}.dat
done
> LC_GCST004748_sCCA3_TWAS.txt
cp LC_GCST004748_sCCA31.dat LC_GCST004748_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 LC_GCST004748_sCCA3${i}.dat >> LC_GCST004748_sCCA3_TWAS.txt
done
rm LC_GCST004748_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/OC_GCST004462.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out OC_GCST004462_sCCA1${CHR}.dat
done
> OC_GCST004462_sCCA1_TWAS.txt
cp OC_GCST004462_sCCA11.dat OC_GCST004462_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 OC_GCST004462_sCCA1${i}.dat >> OC_GCST004462_sCCA1_TWAS.txt
done
rm OC_GCST004462_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/OC_GCST004462.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out OC_GCST004462_sCCA2${CHR}.dat
done
> OC_GCST004462_sCCA2_TWAS.txt
cp OC_GCST004462_sCCA21.dat OC_GCST004462_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 OC_GCST004462_sCCA2${i}.dat >> OC_GCST004462_sCCA2_TWAS.txt
done
rm OC_GCST004462_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/OC_GCST004462.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out OC_GCST004462_sCCA3${CHR}.dat
done
> OC_GCST004462_sCCA3_TWAS.txt
cp OC_GCST004462_sCCA31.dat OC_GCST004462_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 OC_GCST004462_sCCA3${i}.dat >> OC_GCST004462_sCCA3_TWAS.txt
done
rm OC_GCST004462_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/RCC_GCST90320057.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out RCC_GCST90320057_sCCA1${CHR}.dat
done
> RCC_GCST90320057_sCCA1_TWAS.txt
cp RCC_GCST90320057_sCCA11.dat RCC_GCST90320057_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 RCC_GCST90320057_sCCA1${i}.dat >> RCC_GCST90320057_sCCA1_TWAS.txt
done
rm RCC_GCST90320057_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/RCC_GCST90320057.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out RCC_GCST90320057_sCCA2${CHR}.dat
done
> RCC_GCST90320057_sCCA2_TWAS.txt
cp RCC_GCST90320057_sCCA21.dat RCC_GCST90320057_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 RCC_GCST90320057_sCCA2${i}.dat >> RCC_GCST90320057_sCCA2_TWAS.txt
done
rm RCC_GCST90320057_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/RCC_GCST90320057.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out RCC_GCST90320057_sCCA3${CHR}.dat
done
> RCC_GCST90320057_sCCA3_TWAS.txt
cp RCC_GCST90320057_sCCA31.dat RCC_GCST90320057_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 RCC_GCST90320057_sCCA3${i}.dat >> RCC_GCST90320057_sCCA3_TWAS.txt
done
rm RCC_GCST90320057_sCCA3{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/THCA_GCST90399736.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA1.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out THCA_GCST90399736_sCCA1${CHR}.dat
done
> THCA_GCST90399736_sCCA1_TWAS.txt
cp THCA_GCST90399736_sCCA11.dat THCA_GCST90399736_sCCA1_TWAS.txt
for i in {2..22}; do
tail -n +2 THCA_GCST90399736_sCCA1${i}.dat >> THCA_GCST90399736_sCCA1_TWAS.txt
done
rm THCA_GCST90399736_sCCA1{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/THCA_GCST90399736.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA2.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out THCA_GCST90399736_sCCA2${CHR}.dat
done
> THCA_GCST90399736_sCCA2_TWAS.txt
cp THCA_GCST90399736_sCCA21.dat THCA_GCST90399736_sCCA2_TWAS.txt
for i in {2..22}; do
tail -n +2 THCA_GCST90399736_sCCA2${i}.dat >> THCA_GCST90399736_sCCA2_TWAS.txt
done
rm THCA_GCST90399736_sCCA2{1..22}.dat

for CHR in {1..22}; do
Rscript FUSION.assoc_test.R \
--sumstats ./A/THCA_GCST90399736.${CHR}.sumstats \
--weights ./sCCA_weights_v8/sCCA3.pos \
--weights_dir ./sCCA_weights_v8/ \
--ref_ld_chr ./LDREF/UK10K_1KG_qc_rsid_ \
--chr ${CHR} \
--out THCA_GCST90399736_sCCA3${CHR}.dat
done
> THCA_GCST90399736_sCCA3_TWAS.txt
cp THCA_GCST90399736_sCCA31.dat THCA_GCST90399736_sCCA3_TWAS.txt
for i in {2..22}; do
tail -n +2 THCA_GCST90399736_sCCA3${i}.dat >> THCA_GCST90399736_sCCA3_TWAS.txt
done
rm THCA_GCST90399736_sCCA3{1..22}.dat

library(GenomicSEM)
#######sCCA----
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA1_TWAS.txt","BRCA_sCCA1_TWAS.txt","ESCA_sCCA1_TWAS.txt","EC_sCCA1_TWAS.txt", "OC_sCCA1_TWAS.txt", "RCC_sCCA1_TWAS.txt", "CRC_sCCA1_TWAS.txt","THCA_sCCA1_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_TSEM_sCCA<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)
gc()
write_tsv(gfactor_TSEM_sCCA,file = "TSEM_sCCA1_common_factor.txt")

library(GenomicSEM)
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA1_TWAS.txt","BRCA_sCCA1_TWAS.txt","ESCA_sCCA1_TWAS.txt","EC_sCCA1_TWAS.txt", "OC_sCCA1_TWAS.txt", "RCC_sCCA1_TWAS.txt", "CRC_sCCA1_TWAS.txt","THCA_sCCA1_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
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
               model = model1, printwarn = TRUE, cores = 10,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 10,sub = c("e ~ Gene"),
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
write_tsv(A,file = "sCCA1_e_factor.txt")
write_tsv(F1,file = "sCCA1_F1.txt")
write_tsv(F2,file = "sCCA1_F2.txt")
write_tsv(F3,file = "sCCA1_F3.txt")
#######
A<-vroom("sCCA1_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
B<-subset(A, A$Q_SNP_pval_F1 >= (0.05/37920/5) & A$Q_SNP_pval_F2 >= (0.05/37920/5) & A$Q_SNP_pval_F3 >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "e_factor_sCCA1_T_SEM.txt")

A<-vroom("sCCA1_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F1_sCCA1_T_SEM.txt")

A<-vroom("sCCA1_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F2_sCCA1_T_SEM.txt")

A<-vroom("sCCA1_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F3_sCCA1_T_SEM.txt")

A<-vroom("TSEM_sCCA1_common_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_pval >= (0.05/37920/5))

write_tsv(A,file = "CF_sCCA1_T_SEM.txt")
gc()
#######sCCA----
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA2_TWAS.txt","BRCA_sCCA2_TWAS.txt","ESCA_sCCA2_TWAS.txt","EC_sCCA2_TWAS.txt", "OC_sCCA2_TWAS.txt", "RCC_sCCA2_TWAS.txt", "CRC_sCCA2_TWAS.txt","THCA_sCCA2_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_TSEM_sCCA<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)
gc()
write_tsv(gfactor_TSEM_sCCA,file = "TSEM_sCCA2_common_factor.txt")

library(GenomicSEM)
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA2_TWAS.txt","BRCA_sCCA2_TWAS.txt","ESCA_sCCA2_TWAS.txt","EC_sCCA2_TWAS.txt", "OC_sCCA2_TWAS.txt", "RCC_sCCA2_TWAS.txt", "CRC_sCCA2_TWAS.txt","THCA_sCCA2_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
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
               model = model1, printwarn = TRUE, cores = 10,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 10,sub = c("e ~ Gene"),
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
write_tsv(A,file = "sCCA2_e_factor.txt")
write_tsv(F1,file = "sCCA2_F1.txt")
write_tsv(F2,file = "sCCA2_F2.txt")
write_tsv(F3,file = "sCCA2_F3.txt")
#######
A<-vroom("sCCA2_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
B<-subset(A, A$Q_SNP_pval_F1 >= (0.05/37920/5) & A$Q_SNP_pval_F2 >= (0.05/37920/5) & A$Q_SNP_pval_F3 >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "e_factor_sCCA2_T_SEM.txt")

A<-vroom("sCCA2_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F1_sCCA2_T_SEM.txt")

A<-vroom("sCCA2_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F2_sCCA2_T_SEM.txt")

A<-vroom("sCCA2_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F3_sCCA2_T_SEM.txt")

A<-vroom("TSEM_sCCA2_common_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_pval >= (0.05/37920/5))

write_tsv(A,file = "CF_sCCA2_T_SEM.txt")
gc()
#######sCCA----
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA3_TWAS.txt","BRCA_sCCA3_TWAS.txt","ESCA_sCCA3_TWAS.txt","EC_sCCA3_TWAS.txt", "OC_sCCA3_TWAS.txt", "RCC_sCCA3_TWAS.txt", "CRC_sCCA3_TWAS.txt","THCA_sCCA3_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
load("D:/pan_cancer/ldsc.covstruct.Rdata")
gfactor_TSEM_sCCA<-commonfactorGWAS(covstruc=ldsc.covstruct, toler= 1e-60,SNPs=gfactor_genes,parallel=T,TWAS=TRUE)
gc()
write_tsv(gfactor_TSEM_sCCA,file = "TSEM_sCCA3_common_factor.txt")

library(GenomicSEM)
setwd("D:/pan_cancer/T-SEM")
files<-list("LC_sCCA3_TWAS.txt","BRCA_sCCA3_TWAS.txt","ESCA_sCCA3_TWAS.txt","EC_sCCA3_TWAS.txt", "OC_sCCA3_TWAS.txt", "RCC_sCCA3_TWAS.txt", "CRC_sCCA3_TWAS.txt","THCA_sCCA3_TWAS.txt")
trait.names=c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA")
N=c(84804,270181,14595,50773,59713,100080,92392,26347)
binary=c(T,T,T,T,T,T,T,T)
gfactor_genes<- read_fusion(files=files,trait.names=trait.names,N=N,binary=binary)
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
               model = model1, printwarn = TRUE, cores = 10,sub = c("F1 ~ Gene","F2 ~ Gene","F3 ~ Gene"),
               toler = 1e-60,  parallel = TRUE,TWAS=TRUE,Q_SNP=TRUE,smooth_check=TRUE,fix_measurement=TRUE,std.lv = FALSE)
model1 <- '
  F1 =~ LC + ESCA + CRC
  F2 =~ RCC + THCA 
  F3 =~ OC + EC + BRCA
  e  =~ F1 + F2 + F3
  e ~ Gene
'

fit1<-userGWAS(covstruc = ldsc.covstruct, SNPs = gfactor_genes, estimation = "DWLS",
               model = model1, printwarn = TRUE, cores = 10,sub = c("e ~ Gene"),
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
write_tsv(A,file = "sCCA3_e_factor.txt")
write_tsv(F1,file = "sCCA3_F1.txt")
write_tsv(F2,file = "sCCA3_F2.txt")
write_tsv(F3,file = "sCCA3_F3.txt")
#######
setwd("./新建文件夹")
A<-vroom("sCCA3_e_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval_F1 >= (0.05/37920/5) & A$Q_SNP_pval_F2 >= (0.05/37920/5) & A$Q_SNP_pval_F3 >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "e_factor_sCCA3_T_SEM.txt")

A<-vroom("sCCA3_F1.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F1_sCCA3_T_SEM.txt")

A<-vroom("sCCA3_F2.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F2_sCCA3_T_SEM.txt")

A<-vroom("sCCA3_F3.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_SNP_pval >= (0.05/37920/5))
A<-subset(A, A$HSQ >= 0.01)

write_tsv(A,file = "F3_sCCA3_T_SEM.txt")

A<-vroom("TSEM_sCCA3_common_factor.txt")
A<-subset(A, A$Pval_Estimate <= (0.05/37920/5))
A<-subset(A, A$Q_pval >= (0.05/37920/5))

write_tsv(A,file = "CF_sCCA3_T_SEM.txt")
gc()
