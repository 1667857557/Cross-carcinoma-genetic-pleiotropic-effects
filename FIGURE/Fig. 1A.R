library(dplyr)
library(igraph)
library(Hmisc)
library(scales)

setwd("G:/pan_cancer/FIGURE 3")
rg_fname <- "input_rg_mine.txt"
traitlist_fname <- "mine_traitlist.txt"

P_CUTOFF <- 0.05; Q_CUTOFF <- 0.05
MIN_VSIZE <- 6; MAX_VSIZE <- 20
MIN_EWIDTH <- 0.5; MAX_EWIDTH <- 10
NODE_LEGEND_DIAM_MM <- c(`0.1` = 6, `0.3` = 12, `0.5` = 19)
LEGEND_TITLE_CEX <- 2.05
LEGEND_TICK_CEX <- 1.55
LEGEND_TEXT_CEX <- 1.60

dir.create("output", showWarnings = FALSE, recursive = TRUE)
out_prefix <- file.path("output", "network_rg_")

rg <- read.table(rg_fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
traitlist <- read.table(traitlist_fname, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "",
                        fileEncoding = "utf-8", comment.char = "")
if (any(!is.finite(rg$rg))) stop("rg 含有 NA/Inf，请核对输入数据")

g <- graph_from_data_frame(
  transmute(rg, from = p1, to = p2, rg = rg, p = p, q = q, weight = abs(rg)),
  directed = FALSE
)
set.seed(123)
l <- norm_coords(layout_with_fr(g, weights = E(g)$weight))
V(g)$x <- l[, 1]; V(g)$y <- l[, 2]

V(g)$color <- traitlist$COLOR[match(V(g)$name, traitlist$TRAIT)]
h2_vec <- traitlist$H2[match(V(g)$name, traitlist$TRAIT)]
V(g)$size <- scales::rescale(h2_vec, to = c(MIN_VSIZE, MAX_VSIZE),
                              from = range(h2_vec, na.rm = TRUE))
V(g)$label <- V(g)$name
V(g)$label.cex <- 2
V(g)$label.color <- "black"
V(g)$label.font <- 2
V(g)$label.family <- "sans"
V(g)$label.degree <- pi/2
V(g)$label.dist <- 1

blue_grad <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(101)
red_grad <- colorRampPalette(c("#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(101)
# 颜色固定按 |rg|=0..1 映射，保证渐变尺 -1/0/1 的刻度含义准确
if (any(abs(E(g)$rg) > 1)) warning("部分 |rg| > 1：仅绘图颜色截断至 ±1")
edge_rg_clipped <- pmax(-1, pmin(1, E(g)$rg))
edge_idx <- pmin(101L, 1L + round(100 * abs(edge_rg_clipped)))
E(g)$color <- ifelse(edge_rg_clipped >= 0, red_grad[edge_idx], blue_grad[edge_idx])
E(g)$lty <- ifelse(E(g)$q < Q_CUTOFF, 1, 2)
E(g)$width <- scales::rescale(abs(E(g)$rg), to = c(MIN_EWIDTH, MAX_EWIDTH),
                              from = range(abs(E(g)$rg), na.rm = TRUE))
E(g)$curved <- 0.2
E(g)$arrow.size <- 0.5

# 左侧 10x10 英寸：完全沿用当前脚本的网络绘图参数
plot_network <- function() {
  par(mar = c(0, 0, 0, 0))
  plot(g, rescale = FALSE, vertex.frame.color = NA,
       margin = 0, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05))
}

plot_thin_taper_bar <- function(x0 = 0.60, x1 = 4.16, y = 7.50,
                               end_half_height = 0.065, middle_half_height = 0.017) {
  n <- 400L
  xx <- seq(x0, x1, length.out = n + 1L)
  hh <- middle_half_height + (end_half_height - middle_half_height) *
    abs(seq(-1, 1, length.out = n + 1L))
  cmap <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                             "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(n)
  for (i in seq_len(n)) {
    polygon(c(xx[i], xx[i + 1L], xx[i + 1L], xx[i]),
            c(y - hh[i], y - hh[i + 1L], y + hh[i + 1L], y + hh[i]),
            col = cmap[i], border = NA)
  }
}

# 右侧 5x10 英寸，坐标单位即英寸；三个圆采用固定实际直径
plot_legend <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 5), ylim = c(0, 10), xaxs = "i", yaxs = "i", asp = 1)
  bar_center_x <- (0.60 + 4.16)/2
  text(bar_center_x, 8.22, "Genetic correlation", font = 2,
       adj = c(0.5, 0.5), cex = LEGEND_TITLE_CEX)
  plot_thin_taper_bar()
  tick_values <- c(-1, -0.5, 0, 0.5, 1)
  tick_x <- 0.60 + (tick_values + 1)/2 * (4.16 - 0.60)
  segments(tick_x, 7.40, tick_x, 7.31, lwd = 0.85)
  text(tick_x, 7.11, c("-1", "-0.5", "0", "0.5", "1"), cex = LEGEND_TICK_CEX)
  segments(0.95, 6.60, 1.85, 6.60, lwd = 2.5, lty = 1, col = "black")
  text(2.08, 6.60, "q < 0.05", adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)
  segments(0.95, 6.13, 1.85, 6.13, lwd = 2.5, lty = 2, col = "grey65")
  text(2.08, 6.13, "q > 0.05", adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)
  text(bar_center_x, 5.55, "Heritability", font = 2,
       adj = c(0.5, 0.5), cex = LEGEND_TITLE_CEX)
  radius_in <- as.numeric(NODE_LEGEND_DIAM_MM)/(2 * 25.4)
  y <- c(4.83, 3.95 + 0.37, 2.79 + 0.8)
  symbols(rep(1.88, 3), y, circles = radius_in, inches = FALSE,
          bg = "black", fg = NA, add = TRUE)
  text(rep(2.53, 3), y, names(NODE_LEGEND_DIAM_MM),
       adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)
}

plot_all <- function() {
  par(fig = c(0, 10/15, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
  plot_network()
  par(fig = c(10/15, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot_legend()
}

png(paste0(out_prefix, ".png"), width = 15, height = 10, units = "in", res = 300)
plot_all(); dev.off()
pdf(paste0(out_prefix, ".pdf"), width = 15, height = 10)
plot_all(); dev.off()
save(g, l, file = paste0(out_prefix, ".RData"))
