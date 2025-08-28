cd /mnt/g/linux/FLAMES
source activate /home/huang/anaconda3/envs/FLAMES
./magma \
--bfile /mnt/g/Database/1000G/UK10K_1KG_qc_rsid \
--gene-annot ./magma_0kb.genes.annot \
--pval pan_cancer_magma.txt ncol=N \
--gene-model snp-wise=mean \
--out pan_cancer

./magma/magma \
--gene-results ./pan_cancer.genes.raw \
--gene-covar ./gtex_v8_ts_avg_log2TPM.txt \
--out tissue_pan_cancer

export PATH=$PATH:/mnt/g/linux/FLAMES/
  python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/pan_cancer \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/pan_cancer


python FLAMES.py annotate \
-o /mnt/g/linux/FLAMES/ \
-a /mnt/g/linux/FLAMES/Annotation_data/ \
-p /mnt/g/linux/FLAMES/pan_cancer.preds \
-m /mnt/g/linux/FLAMES/pan_cancer.genes.out \
-mt /mnt/g/linux/FLAMES/tissue_pan_cancer.gsa.out \
-id /mnt/g/linux/FLAMES/locus_files/common_factor_indexfile.txt

python FLAMES.py FLAMES \
-id /mnt/g/linux/FLAMES/locus_files/common_factor_indexfile.txt \
-o /mnt/g/linux/FLAMES/
  
  
  ######
./magma \
--bfile /mnt/g/Database/1000G/UK10K_1KG_qc_rsid \
--gene-annot ./magma_0kb.genes.annot \
--pval e_factor_magma.txt ncol=N \
--gene-model snp-wise=mean \
--out e_factor

./magma/magma \
--gene-results ./e_factor.genes.raw \
--gene-covar ./gtex_v8_ts_avg_log2TPM.txt \
--out tissue_e_factor

export PATH=$PATH:/mnt/g/linux/FLAMES/
  python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/e_factor \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/e_factor


python FLAMES.py annotate \
-o /mnt/g/linux/FLAMES/ \
-a /mnt/g/linux/FLAMES/Annotation_data/ \
-p /mnt/g/linux/FLAMES/e_factor.preds \
-m /mnt/g/linux/FLAMES/e_factor.genes.out \
-mt /mnt/g/linux/FLAMES/tissue_e_factor.gsa.out \
-id /mnt/g/linux/FLAMES/locus_files/e_factor_indexfile.txt

python FLAMES.py FLAMES \
-id /mnt/g/linux/FLAMES/locus_files/e_factor_indexfile.txt \
-o /mnt/g/linux/FLAMES/
  
  #####
./magma \
--bfile /mnt/g/Database/1000G/UK10K_1KG_qc_rsid \
--gene-annot ./magma_0kb.genes.annot \
--pval F1_magma.txt ncol=N \
--gene-model snp-wise=mean \
--out F1

./magma/magma \
--gene-results ./F1.genes.raw \
--gene-covar ./gtex_v8_ts_avg_log2TPM.txt \
--out tissue_F1

export PATH=$PATH:/mnt/g/linux/FLAMES/
  python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/F1 \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/F1


python FLAMES.py annotate \
-o /mnt/g/linux/FLAMES/ \
-a /mnt/g/linux/FLAMES/Annotation_data/ \
-p /mnt/g/linux/FLAMES/F1.preds \
-m /mnt/g/linux/FLAMES/F1.genes.out \
-mt /mnt/g/linux/FLAMES/tissue_F1.gsa.out \
-id /mnt/g/linux/FLAMES/locus_files/F1_indexfile.txt

python FLAMES.py FLAMES \
-id /mnt/g/linux/FLAMES/locus_files/F1_indexfile.txt \
-o /mnt/g/linux/FLAMES/
  #######
./magma \
--bfile /mnt/g/Database/1000G/UK10K_1KG_qc_rsid \
--gene-annot ./magma_0kb.genes.annot \
--pval F2_magma.txt ncol=N \
--gene-model snp-wise=mean \
--out F2

./magma/magma \
--gene-results ./F2.genes.raw \
--gene-covar ./gtex_v8_ts_avg_log2TPM.txt \
--out tissue_F2

export PATH=$PATH:/mnt/g/linux/FLAMES/
  python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/F2 \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/F2


python FLAMES.py annotate \
-o /mnt/g/linux/FLAMES/FLAMES_annotated/ \
-a /mnt/g/linux/FLAMES/Annotation_data/ \
-p /mnt/g/linux/FLAMES/F2.preds \
-m /mnt/g/linux/FLAMES/F2.genes.out \
-mt /mnt/g/linux/FLAMES/tissue_F2.gsa.out \
-id /mnt/g/linux/FLAMES/locus_files/F2_indexfile.txt

python FLAMES.py FLAMES \
-id /mnt/g/linux/FLAMES/locus_files/F2_indexfile.txt \
-o /mnt/g/linux/FLAMES/
  ######
./magma \
--bfile /mnt/g/Database/1000G/UK10K_1KG_qc_rsid \
--gene-annot ./magma_0kb.genes.annot \
--pval F3_magma.txt ncol=N \
--gene-model snp-wise=mean \
--out F3

./magma/magma \
--gene-results ./F3.genes.raw \
--gene-covar ./gtex_v8_ts_avg_log2TPM.txt \
--out tissue_F3

export PATH=$PATH:/mnt/g/linux/FLAMES/
  python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/F3 \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/F3


python FLAMES.py annotate \
-o /mnt/g/linux/FLAMES/ \
-a /mnt/g/linux/FLAMES/Annotation_data/ \
-p /mnt/g/linux/FLAMES/F3.preds \
-m /mnt/g/linux/FLAMES/F3.genes.out \
-mt /mnt/g/linux/FLAMES/tissue_F3.gsa.out \
-id /mnt/g/linux/FLAMES/locus_files/F3_indexfile.txt

python FLAMES.py FLAMES \
-id /mnt/g/linux/FLAMES/locus_files/F3_indexfile.txt \
-o /mnt/g/linux/FLAMES/
  
  