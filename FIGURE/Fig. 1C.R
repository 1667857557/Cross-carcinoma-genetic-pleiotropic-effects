library(corrplot)

# 1. 构建相关性矩阵 cor
vars <- c("LC","BRCA","ESCA","EC","OC","RCC","CRC","THCA","PRCA")
vals <- c(
  1,
  0.2197,1,
  0.2598,0.0045,1,
  0.1909,0.2228,0.0998,1,
  0.1450,0.1963,0.1174,0.4386,1,
  0.4112,0.2056,0.0774,0.2937,0.1597,1,
  0.2890,0.1423,0.2250,0.1875,0.0519,0.1748,1,
  0.1206,0.1009,0.1156,0.0912,0.1247,0.3510,0.0494,1,
  -0.0371,0.0942,0.0573,0.0133,-0.0118,0.1257,0.0407,0.0909,1
)
m <- matrix(NA, 9, 9)
m[upper.tri(m, diag=TRUE)] <- vals
m[lower.tri(m)] <- t(m)[lower.tri(m)]
rownames(m) <- colnames(m) <- vars
cor <- m

# 2. 构建显著 locus 计数矩阵 n_locus
n_locus <- matrix(0, 9, 9, dimnames = list(vars, vars))
n_locus["LC","BRCA"] <- 1;    n_locus["LC","RCC"]  <- 2
n_locus["LC","CRC"]  <- 2;    n_locus["LC","PRCA"] <- 6
n_locus["BRCA","EC"] <- 1;    n_locus["BRCA","RCC"]  <- 1
n_locus["BRCA","CRC"]  <- 1;    n_locus["BRCA","THCA"] <- 1
n_locus["BRCA","PRCA"] <- 8
n_locus["EC","PRCA"] <- 1
n_locus["RCC","CRC"]   <- 1;    n_locus["RCC","THCA"] <- 2
n_locus["RCC","PRCA"]  <- 4
n_locus["CRC","THCA"]  <- 1;    n_locus["CRC","PRCA"] <- 3
n_locus["THCA","PRCA"] <- 1
n_locus[lower.tri(n_locus)] <- t(n_locus)[lower.tri(n_locus)]
n_locus[n_locus == 0] <- NA

# 3. 同步做层次聚类重排
ord  <- corrMatOrder(cor, order = "hclust")
cor2 <- cor[ord, ord]
n2   <- n_locus[ord, ord]
gc()
# 4. 绘图并叠加数字（所有字体都放大 1 号，标签往上挪）
col2 <- colorRampPalette(c("#053061","white","#67001F"), alpha = TRUE)
pdf("heatmap_LDSC_LAVA.pdf", width = 11, height = 11)


corrplot(cor2,
         order     = "original",
         col       = col2(100),
         method    = "square",
         cl.length = 5,
         type      = "lower",
         diag      = FALSE,
         tl.col    = "black",
         tl.cex    = 1.4,    # 标签字体增大
         tl.srt    = 0,      # 水平放置
         tl.offset = 1,    # 往上偏移
         cl.cex    = 1.4,    # 右侧图例文字增大
         cl.pos    = "r",
         cl.ratio  = 0.2)

# 5. 叠加 locus 数字（字体增大）
n <- length(vars)
for(i in seq_len(n)) {
  for(j in seq_len(i)) {
    cnt <- n2[i, j]
    if(!is.na(cnt)) {
      text(j, n - i + 1,
           labels = cnt,
           cex    = 1.4)  # 数字字体增大
    }
  }
}
## 输出 ------------------------------------------------------------------------
dev.off()

