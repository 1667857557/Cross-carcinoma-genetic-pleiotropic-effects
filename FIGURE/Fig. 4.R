######Manhattan-------
library(topr)
setwd("D:/pan_cancer")
pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils","Seurat")
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())

common_factor<-vroom("common_factor_GWAS.txt.gz")
common_factor$P = pchisq((common_factor$BETA/common_factor$SE)^2, df=1, lower.tail=FALSE)
common_factor<-select(common_factor,c("CHR","POS","SNP","other_allele","effect_allele","P","BETA"))
colnames(common_factor)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
common_factor$P[common_factor$P < 2.2e-308] <- 1e-300

e_factor<-vroom("e_factor.txt.gz")
e_factor<-select(e_factor,c("CHR","POS","SNP","other_allele","effect_allele","P","BETA"))
colnames(e_factor)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
e_factor$P[e_factor$P < 2.2e-308] <- 1e-300

F1<-vroom("./F1.txt.gz")
F1<-select(F1,c("CHR","POS","SNP","other_allele","effect_allele","P","BETA"))
colnames(F1)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
F1$P[F1$P < 2.2e-308] <- 1e-300

F2<-vroom("./F2.txt.gz")
F2<-select(F2,c("CHR","POS","SNP","other_allele","effect_allele","P","BETA"))
colnames(F2)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
F2$P[F2$P < 2.2e-308] <- 1e-300

F3<-vroom("./F3.txt.gz")
F3<-select(F3,c("CHR","POS","SNP","other_allele","effect_allele","P","BETA"))
colnames(F3)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
F3$P[F3$P < 2.2e-308] <- 1e-300


gc()
# pdf('common_factor.pdf', width = 10, height = 3.5,compress = TRUE)
png("common_factor.png", width = 10, height = 3.5, units = "in", res = 600, bg = "white")
lead_snps <- get_lead_snps(common_factor)
manhattan(list(common_factor, lead_snps), 
          annotate=c(5e-100, 1e-8), 
          color=c("#996699","grey"), 
          legend_labels=NULL,
          size=c(1,2),
          shape=c(19,8), label_size = 3,build = 37,
          sign_thresh_color = "black",
          sign_thresh_size = 1,
          angle=60, 
          show_legend = FALSE,
          label_fontface = "bold.italic", 
          nudge_y = 5, 
          scale=0.8, 
          legend_position = "top", 
          ymax=22,  
          sign_thresh_label_size = 4,
          even_no_chr_lightness = c(0.8,0.5),
          verbose=F,
          title="") 
dev.off()
# pdf('e_factor.pdf', width = 10, height = 3,compress = TRUE)

png("e_factor.png", width = 10, height = 3, units = "in", res = 600, bg = "white")
lead_snps <- get_lead_snps(e_factor)
manhattan(list(e_factor, lead_snps), 
          annotate=c(5e-100, 5e-10), 
          color=c("#008B8B","grey"), 
          legend_labels=NULL,
          size=c(1,2),
          shape=c(19,8), 
          sign_thresh_color = "black",label_size = 3,build = 37,
          sign_thresh_size = 1,
          angle=60, 
          show_legend = FALSE,
          label_fontface = "bold.italic", 
          nudge_y = 5, 
          scale=0.8, 
          legend_position = "top", 
          ymax=20,  
          sign_thresh_label_size = 4,
          even_no_chr_lightness = c(0.8,0.5),
          verbose=F,
          title="") 
dev.off()

df<-F1
df$logP_original <- -log10(df$P)
print("原始数据和-log10转换:")
print(df)
df$logP_adjusted <- df$logP_original
indices_to_scale <- which(df$logP_original > 20)
if (length(indices_to_scale) > 0) {
  
  old_values <- df$logP_original[indices_to_scale]
  
  old_min <- 20
  old_max <- max(old_values)
  
  new_min <- 20.0
  new_max <- 24.9
  if (old_max > old_min) {
    scaled_values <- new_min + (old_values - old_min) * (new_max - new_min) / (old_max - old_min)
  } else {
    scaled_values <- rep(new_max, length(old_values))
  }
  
  df$logP_adjusted[indices_to_scale] <- scaled_values
}
df$P_adjusted <- 10^(-df$logP_adjusted)
df<-select(df,c("CHROM","POS","ID","REF","ALT","P_adjusted","BETA"))
colnames(df)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
print("调整后的最终数据框:")
png("F1.png", width = 10, height = 4.5, units = "in", res = 600, bg = "white")

# pdf('F1.pdf', width = 10, height = 5,compress = TRUE)
lead_snps <- get_lead_snps(df)
manhattan(list(df, lead_snps), 
          annotate=c(5e-100, 5e-10), 
          color=c("#B22222","grey"), 
          legend_labels=NULL,
          size=c(1,2),
          shape=c(19,8), 
          sign_thresh_color = "black",
          sign_thresh_size = 1,label_size = 3,build = 37,
          angle=60, 
          show_legend = FALSE,
          label_fontface = "bold.italic", 
          nudge_y = 5, 
          scale=0.8, 
          legend_position = "top", 
          ymax=32,  
          sign_thresh_label_size = 4,
          even_no_chr_lightness = c(0.8,0.5),
          verbose=F,
          title="") 
dev.off()

df<-F2
df$logP_original <- -log10(df$P)
print("原始数据和-log10转换:")
print(df)
df$logP_adjusted <- df$logP_original
indices_to_scale <- which(df$logP_original > 20)
if (length(indices_to_scale) > 0) {
  
  old_values <- df$logP_original[indices_to_scale]
  
  old_min <- 20
  old_max <- max(old_values)
  
  new_min <- 20.0
  new_max <- 24.9
  if (old_max > old_min) {
    scaled_values <- new_min + (old_values - old_min) * (new_max - new_min) / (old_max - old_min)
  } else {
    scaled_values <- rep(new_max, length(old_values))
  }
  
  df$logP_adjusted[indices_to_scale] <- scaled_values
}
df$P_adjusted <- 10^(-df$logP_adjusted)
df<-select(df,c("CHROM","POS","ID","REF","ALT","P_adjusted","BETA"))
colnames(df)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
print("调整后的最终数据框:")
# pdf('F2.pdf', width = 10, height = 3.5,compress = TRUE)
png("F2.png", width = 10, height = 4.5, units = "in", res = 600, bg = "white")

lead_snps <- get_lead_snps(df)
manhattan(list(df, lead_snps), 
          annotate=c(5e-100, 5e-10), 
          color=c("#4169E1","grey"), 
          legend_labels=NULL,
          size=c(1,2),
          shape=c(19,8), label_size = 3,build = 37,
          sign_thresh_color = "black",
          sign_thresh_size = 1,
          angle=60, 
          show_legend = FALSE,
          label_fontface = "bold.italic", 
          nudge_y = 5, 
          scale=0.8, 
          legend_position = "top", 
          ymax=32,  
          sign_thresh_label_size = 4,
          even_no_chr_lightness = c(0.8,0.5),
          verbose=F,
          title="") 
dev.off()


df<-F3
df$logP_original <- -log10(df$P)
print("原始数据和-log10转换:")
print(df)
df$logP_adjusted <- df$logP_original
indices_to_scale <- which(df$logP_original > 20)
if (length(indices_to_scale) > 0) {
  
  old_values <- df$logP_original[indices_to_scale]
  
  old_min <- 20
  old_max <- max(old_values)
  
  new_min <- 20.0
  new_max <- 24.9
  if (old_max > old_min) {
    scaled_values <- new_min + (old_values - old_min) * (new_max - new_min) / (old_max - old_min)
  } else {
    scaled_values <- rep(new_max, length(old_values))
  }
  
  df$logP_adjusted[indices_to_scale] <- scaled_values
}
df$P_adjusted <- 10^(-df$logP_adjusted)
df<-select(df,c("CHROM","POS","ID","REF","ALT","P_adjusted","BETA"))
colnames(df)<-c("CHROM","POS","ID","REF","ALT","P","BETA")
print("调整后的最终数据框:")
png("F3.png", width = 10, height = 4.5, units = "in", res = 600, bg = "white")

# pdf('F3.pdf', width = 10, height = 4.5,compress = TRUE)
lead_snps <- get_lead_snps(df)
manhattan(list(df, lead_snps), 
          annotate=c(5e-100, 5e-10), 
          color=c("#FF8C00","grey"), 
          legend_labels=NULL,
          size=c(1,2),
          shape=c(19,8), label_size = 3,build = 37,
          sign_thresh_color = "black",
          sign_thresh_size = 1,
          angle=60, 
          show_legend = FALSE,
          label_fontface = "bold.italic", 
          nudge_y = 5, 
          scale=0.8, 
          legend_position = "top", 
          ymax=35,  
          sign_thresh_label_size = 4,
          even_no_chr_lightness = c(0.8,0.5),
          verbose=F,
          title="") 
dev.off()
######UMAP------
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.neighbors import KernelDensity
import scanpy as sc
import os

output_dir = "umap_plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

traits = ['F1_GWAS', 'F2_GWAS', 'F3_GWAS', 'CF_GWAS', 'EF_GWAS']

for trait in traits:
    fig, ax = plt.subplots(figsize=(6, 6))
    
    sc.pl.umap(
        adata,
        color=trait,
        color_map="RdBu_r",
        vmin=-5,
        vmax=5,
        s=20,
        ax=ax,
        show=False,
        frameon=False
    )
    
    mask = adata.obs['curated_anno'].isin(['Malignant', 'Epithelial'])
    
    if mask.sum() > 10:
        x = adata.obsm["X_umap"][mask, 0]
        y = adata.obsm["X_umap"][mask, 1]
        
        x_min, x_max = x.min() - 1.0, x.max() + 1.0
        y_min, y_max = y.min() - 1.0, y.max() + 1.0
        xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
        
        xy_train = np.vstack([x, y]).T
        xy_test = np.vstack([xx.ravel(), yy.ravel()]).T
        
        kde = KernelDensity(bandwidth=0.5, metric='euclidean')
        kde.fit(xy_train)
        
        Z = np.exp(kde.score_samples(xy_test))
        Z = Z.reshape(xx.shape)
        
        ax.contour(
            xx, yy, Z,
            levels=[Z.max() * 0.05],
            colors='#ea3433',
            linestyles='-',
            linewidths=3.0,
            alpha=0.9
        )
    
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='#ea3433', lw=3, label='Epithelial (Malignant and Normal)')]
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(0.0, -0.05), frameon=False, fontsize=14)
    
    plt.tight_layout()
    
    plt.savefig(f"{output_dir}/{trait}_umap.pdf", format='pdf', bbox_inches='tight')
    plt.close()

print(f"Finish")
