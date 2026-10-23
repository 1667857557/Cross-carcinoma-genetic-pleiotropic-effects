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
n_locus["BRCA","CRC"]  <- 1;  n_locus["BRCA","THCA"] <- 1
n_locus["BRCA","PRCA"] <- 8
n_locus["EC","PRCA"] <- 1
n_locus["RCC","CRC"]   <- 1;  n_locus["RCC","THCA"] <- 2
n_locus["RCC","PRCA"]  <- 4
n_locus["CRC","THCA"]  <- 1;  n_locus["CRC","PRCA"] <- 3
n_locus["THCA","PRCA"] <- 1
n_locus[lower.tri(n_locus)] <- t(n_locus)[lower.tri(n_locus)]
n_locus[n_locus == 0] <- NA

# 3. 同步做层次聚类重排
ord  <- corrMatOrder(cor, order = "hclust")
cor2 <- cor[ord, ord]
n2   <- n_locus[ord, ord]
gc()

# 4. 绘图并叠加数字，右侧拼接方块大小图注
col2 <- colorRampPalette(c("#0571b0","white","#ca0020"), alpha = TRUE)

dev_w <- 13
dev_h <- 11
pdf("heatmap_LDSC_LAVA.pdf", width = dev_w, height = dev_h)

# ---------- 面板 1：相关性热图 ----------
par(fig = c(0, 0.78, 0, 1), mar = c(0.5, 0.5, 0.5, 0.5), new = FALSE)

corrplot(cor2,
         order     = "original",
         col       = col2(100),
         method    = "square",
         cl.length = 5,
         type      = "lower",
         diag      = FALSE,
         tl.col    = "black",
         tl.cex    = 1.4,
         tl.srt    = 0,
         tl.offset = 1,
         cl.cex    = 1.4,
         cl.pos    = "r",
         cl.ratio  = 0.2)

# 叠加 locus 数字
n <- length(vars)
for(i in seq_len(n)) {
  for(j in seq_len(i)) {
    cnt <- n2[i, j]
    if(!is.na(cnt)) {
      text(j, n - i + 1, labels = cnt, cex = 1.4)
    }
  }
}

# ---- 读取主图 y 方向“每格物理英寸数” ----
usr_main <- par("usr")
pin_main <- par("pin")
ky <- pin_main[2] / (usr_main[4] - usr_main[3])

# ---------- 面板 2：方块大小图注（左移 1 英寸）----------
shift_in <- 1
fig_x <- c(0.80 - shift_in/dev_w, 1.00 - shift_in/dev_w)   # c(0.723, 0.923)
fig_y <- c(0.15, 0.85)

par(fig = c(fig_x, fig_y), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
pin_leg <- par("pin")

# 标题下移，紧贴最上方方块
text(0.5, 0.69, "Square size", cex = 1.4, font = 2)

# 图例档位：0.1, 0.2, 0.4；缩小上下间距并整体居中
r_vals <- c(0.1, 0.2, 0.4)
y_pos  <- c(0.60, 0.47, 0.34)
cx     <- 0.30

size_mode <- "sqrt"   # corrplot squares 为面积映射：边长 ∝ sqrt(|r|)

for(k in seq_along(r_vals)) {
  rv <- r_vals[k]
  s_in <- switch(size_mode,
                 sqrt   = sqrt(rv) * ky,
                 linear = rv * ky)
  hw <- (s_in / 2) / pin_leg[1]
  hh <- (s_in / 2) / pin_leg[2]
  rect(cx - hw, y_pos[k] - hh,
       cx + hw, y_pos[k] + hh,
       border = "black", lwd = 1)
  text(0.62, y_pos[k],
       sprintf("|r| = %.1f", rv),
       adj = 0, cex = 1.1)
}

dev.off()
