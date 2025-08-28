# ──────────────────── 依赖 ────────────────────
library(dplyr)
library(igraph)
library(Hmisc)
library(scales)

setwd("D:/pan_cancer/FIGURE 3")   # 自行修改

# ──────────────────── 1. 输入 ────────────────────
rg_fname        <- "input_GNA_mine.txt"
traitlist_fname <- "mine_traitlist.txt"

# ──────────────────── 2. 常量 ────────────────────
P_CUTOFF   <- 0.05          # 显著性阈值（保留，可调）
MIN_VSIZE  <- 6;  MAX_VSIZE <- 20
MIN_EWIDTH <- 0.5;MAX_EWIDTH<- 10
RG_LIMIT   <- 0.5           # 相关性绝对值理论上界

# ──────────────────── 3. 读数据 ────────────────────
rg <- read.table(rg_fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rg$q <- p.adjust(rg$p, method = "BH")          # 如已有可省略

traitlist <- read.table(traitlist_fname, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "",
                        fileEncoding = "utf-8", comment.char = "")

# 若发现 rg 超出 ±0.5，给出提示并截断
if (any(abs(rg$rg) > RG_LIMIT)) {
  warning("发现 |rg| > 0.5，已被截断至 ±0.5 以匹配线宽标度。")
  rg$rg <- pmax(pmin(rg$rg,  RG_LIMIT), -RG_LIMIT)
}

# ──────────────────── 4. 构图 ────────────────────
g <- graph_from_data_frame(
  transmute(rg,
            from   = p1,
            to     = p2,
            rg     = rg,
            p      = p,
            q      = q,
            weight = abs(rg)          # 供布局用
  ),
  directed = FALSE
)

# ──────────────────── 5. 布局 ────────────────────
l <- norm_coords(layout_with_fr(g, weights = E(g)$weight))
V(g)$x <- l[,1]; V(g)$y <- l[,2]

# ──────────────────── 6. 节点样式 ────────────────────
V(g)$color <- traitlist$COLOR[match(V(g)$name, traitlist$TRAIT)]
h2_vec     <- traitlist$H2[match(V(g)$name, traitlist$TRAIT)]
V(g)$size  <- scales::rescale(h2_vec, to = c(MIN_VSIZE, MAX_VSIZE),
                      from = range(h2_vec, na.rm = TRUE))

V(g)$label        <- V(g)$name
V(g)$label.cex    <- 2
V(g)$label.color  <- "black"
V(g)$label.font   <- 2
V(g)$label.family <- "sans"
V(g)$label.degree <- pi/2   # 上方
V(g)$label.dist   <- 1

# ──────────────────── 7. 边样式 ────────────────────
## 颜色渐变
blue_grad <- colorRampPalette(c("#FFFFFF","#D1E5F0","#92C5DE",
                                "#4393C3","#2166AC","#053061"))(100)
red_grad  <- colorRampPalette(c("#FFFFFF","#FDDBC7","#F4A582",
                                "#D6604D","#B2182B","#67001F"))(100)

## 依据 |rg| 将值映射到 1–100；固定 from = c(0,0.5)
abs_rg_scaled <- ceiling(scales::rescale(abs(E(g)$rg), to = c(1,100),
                                 from = c(0, RG_LIMIT)))

## 颜色：显著性可选；正红负蓝
E(g)$color <- ifelse(E(g)$q >= 0.05, "grey70",
                     ifelse(E(g)$rg >= 0, red_grad[abs_rg_scaled], blue_grad[abs_rg_scaled]))


## 线宽：固定标度
E(g)$width <- scales::rescale(abs(E(g)$rg), to = c(MIN_EWIDTH, MAX_EWIDTH),
                      from = c(0, RG_LIMIT))

E(g)$lty        <- 1
E(g)$curved     <- 0.2
E(g)$arrow.size <- 0.5  # 无向图时无效，可留作备份

# ──────────────────── 8. 绘图 ────────────────────
plot(g, rescale = FALSE, vertex.frame.color = NA)

# 如需保存：
# png("output/network_GNA.png", 9, 9, units = "in", res = 300)
# plot(g, rescale = FALSE, vertex.frame.color = NA)
# dev.off()


png("output/network_GNA.png", width = 10, height = 10, units = "in", res = 300)
plot(g,
     rescale = FALSE,
     vertex.frame.color = NA)
dev.off()

pdf("output/network_GNA.pdf", width = 12, height = 12)
plot(g,
     rescale = FALSE,
     vertex.frame.color = NA)
dev.off()
