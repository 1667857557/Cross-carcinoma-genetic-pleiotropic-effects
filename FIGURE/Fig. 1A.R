# ──────────────────── 依赖 ────────────────────
library(dplyr);  library(igraph);  library(Hmisc);  library(scales)
setwd("D:/pan_cancer/FIGURE 3")
# ──────────────────── 1. 输入 ────────────────────
rg_fname        <- "input_rg_mine.txt"
traitlist_fname <- "mine_traitlist.txt"

# ──────────────────── 2. 常量 & 读数据 ─────────────
P_CUTOFF  <- 0.05;  Q_CUTOFF <- 0.05
MIN_VSIZE <- 6;     MAX_VSIZE <- 20
MIN_EWIDTH <- 0.5;  MAX_EWIDTH <- 10

rg <- read.table(rg_fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
traitlist <- read.table(traitlist_fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        quote = "", fileEncoding = "utf-8", comment.char = "")

# ──────────────────── 3. 构图 ────────────────────
g <- graph_from_data_frame(
  transmute(rg, from = p1, to = p2, rg = rg, p = p, q = q,
            weight = abs(rg)),
  directed = FALSE
)

# ──────────────────── 4. 布局 ────────────────────
l <- norm_coords(layout_with_fr(g, weights = E(g)$weight))
V(g)$x <- l[,1];  V(g)$y <- l[,2]

# ──────────────────── 5. 节点样式 ────────────────────
V(g)$color <- traitlist$COLOR[match(V(g)$name, traitlist$TRAIT)]
h2_vec     <- traitlist$H2[match(V(g)$name, traitlist$TRAIT)]
V(g)$size  <- scales::rescale(h2_vec, to = c(MIN_VSIZE, MAX_VSIZE),
                      from = range(h2_vec, na.rm = TRUE))

## —— 文字固定大小 & 上方 —— ##
V(g)$label         <- V(g)$name
V(g)$label.cex     <- 2            # ← 固定字号
V(g)$label.color   <- "black"
V(g)$label.font    <- 2
V(g)$label.family  <- "sans"
V(g)$label.degree  <- pi/2         # 上方
V(g)$label.dist    <- 1          # 与节点边缘间距
## ———————————————— ##

# ──────────────────── 6. 边样式 ────────────────────
blue_grad <- colorRampPalette(c("#FFFFFF","#D1E5F0","#92C5DE",
                                "#4393C3","#2166AC","#053061"))(100)
red_grad  <- colorRampPalette(c("#FFFFFF","#FDDBC7","#F4A582",
                                "#D6604D","#B2182B","#67001F"))(100)

abs_rg_scaled <- ceiling(scales::rescale(abs(E(g)$rg), to = c(1,100)))

E(g)$color <- ifelse(E(g)$q >= Q_CUTOFF, "grey70",
                     ifelse(E(g)$rg >= 0, red_grad[abs_rg_scaled], blue_grad[abs_rg_scaled]))

E(g)$lty   <- ifelse(E(g)$p < Q_CUTOFF, 1, 2)
E(g)$width <- scales::rescale(abs(E(g)$rg), to = c(MIN_EWIDTH, MAX_EWIDTH),
                      from = range(abs(E(g)$rg), na.rm = TRUE))
E(g)$curved <- 0.2
E(g)$arrow.size <- 0.5

# ──────────────────── 7. 输出 ────────────────────
if (!dir.exists("output")) dir.create("output")
out_prefix <- file.path("output", paste0("network_rg_"))

png(paste0(out_prefix, ".png"), width = 10, height = 10, units = "in", res = 300)
plot(g,
     rescale = FALSE,
     vertex.frame.color = NA,
     margin = 0,
     xlim = c(-1.05, 1.05),   # ← 四周额外放 5%
     ylim = c(-1.05, 1.05))
dev.off()

pdf(paste0(out_prefix, ".pdf"), width = 12, height = 12)
plot(g,
     rescale = FALSE,
     vertex.frame.color = NA,
     margin = 0,
     xlim = c(-1.05, 1.05),   # ← 四周额外放 5%
     ylim = c(-1.05, 1.05))
dev.off()

save(g, l, file = paste0(out_prefix, ".RData"))

