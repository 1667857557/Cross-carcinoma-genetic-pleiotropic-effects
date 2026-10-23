library(dplyr)
library(igraph)
library(Hmisc)
library(scales)

setwd("G:/pan_cancer/FIGURE 3")
rg_fname <- "input_GNA_mine.txt"
traitlist_fname <- "mine_traitlist.txt"
P_CUTOFF <- 0.05
MIN_VSIZE <- 6; MAX_VSIZE <- 20
MIN_EWIDTH <- 0.5; MAX_EWIDTH <- 10
RG_LIMIT <- 0.5
NODE_LEGEND_DIAM_MM <- c(`0.1` = 6, `0.3` = 12, `0.5` = 19)
LEGEND_TITLE_CEX <- 2.05
LEGEND_TICK_CEX <- 1.55
LEGEND_TEXT_CEX <- 1.60
OUTDIR <- "output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

rg <- read.table(rg_fname, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rg$q <- p.adjust(rg$p, method = "BH")
traitlist <- read.table(traitlist_fname, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "",
                        fileEncoding = "utf-8", comment.char = "")

if (any(!is.finite(rg$rg))) stop("rg 含有非有限值，请检查输入数据。")
if (any(abs(rg$rg) > 1)) warning("输入存在 |rg| > 1：绘图颜色将限制到 ±1，原始 rg 不修改。")
rg_clip <- pmax(pmin(rg$rg, RG_LIMIT), -RG_LIMIT)

g <- graph_from_data_frame(
  transmute(rg, from = p1, to = p2, rg = rg, p = p, q = q, weight = abs(rg_clip)),
  directed = FALSE
)

l <- norm_coords(layout_with_fr(g, weights = E(g)$weight))
V(g)$x <- l[, 1]; V(g)$y <- l[, 2]
V(g)$color <- traitlist$COLOR[match(V(g)$name, traitlist$TRAIT)]
h2_vec <- traitlist$H2[match(V(g)$name, traitlist$TRAIT)]
if (any(!is.finite(h2_vec))) stop("部分节点未匹配到有效 H2，检查 traitlist$TRAIT 与 V(g)$name。")
h2_rng <- range(h2_vec)
if (diff(h2_rng) == 0) stop("所有节点 H2 相同，原脚本的 rescale() 无法确定 6–20 的大小映射。")
V(g)$size <- scales::rescale(h2_vec, to = c(MIN_VSIZE, MAX_VSIZE), from = h2_rng)
V(g)$label <- V(g)$name
V(g)$label.cex <- 2
V(g)$label.color <- "black"
V(g)$label.font <- 2
V(g)$label.family <- "sans"
V(g)$label.degree <- pi / 2
V(g)$label.dist <- 1

blue_grad <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(101)
red_grad <- colorRampPalette(c("#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(101)
edge_val <- pmax(-1, pmin(1, E(g)$rg))
edge_idx <- pmin(101L, 1L + round(100 * abs(edge_val)))
E(g)$color <- ifelse(E(g)$q >= P_CUTOFF, "grey70",
                     ifelse(edge_val >= 0, red_grad[edge_idx], blue_grad[edge_idx]))
E(g)$width <- scales::rescale(pmin(abs(E(g)$rg), RG_LIMIT),
                              to = c(MIN_EWIDTH, MAX_EWIDTH), from = c(0, RG_LIMIT))
E(g)$lty <- ifelse(E(g)$q < P_CUTOFF, 1, 2)
E(g)$curved <- 0.2
E(g)$arrow.size <- 0.5

measure_network <- function() {
  plot(g, rescale = FALSE, vertex.frame.color = NA)
  list(usr = par("usr"), pin = par("pin"))
}

legend_circle_radii_in <- function() {
  as.numeric(NODE_LEGEND_DIAM_MM) / (2 * 25.4)
}

plot_thin_taper_bar <- function(x0 = 0.60, x1 = 4.16, y = 7.50,
                               end_half_height = 0.065, middle_half_height = 0.017) {
  xx <- seq(x0, x1, length.out = 402)
  f <- abs(seq(-1, 1, length.out = length(xx)))
  hh <- middle_half_height + (end_half_height - middle_half_height) * f
  cmap <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF",
                             "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(401)
  for (i in seq_len(401)) {
    x_left <- if (i == 1L) xx[i] else xx[i] - 0.001
    x_right <- if (i == 401L) xx[i + 1] else xx[i + 1] + 0.001
    polygon(c(x_left, x_right, x_right, x_left),
            c(y - hh[i], y - hh[i + 1], y + hh[i + 1], y + hh[i]),
            col = cmap[i], border = NA)
  }
}

draw_legend <- function(network_state) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 5), ylim = c(0, 10), xaxs = "i", yaxs = "i", asp = 1)

  bar_center_x <- (0.60 + 4.16) / 2
  text(bar_center_x, 8.22, "Partial genetic correlation", adj = c(0.5, 0.5), font = 2, cex = LEGEND_TITLE_CEX)
  plot_thin_taper_bar()
  tick_values <- c(-1, -0.5, 0, 0.5, 1)
  tick_positions <- 0.60 + (tick_values + 1) / 2 * (4.16 - 0.60)
  segments(tick_positions, 7.40, tick_positions, 7.31, lwd = 0.85)
  text(tick_positions, 7.11, c("-1", "-0.5", "0", "0.5", "1"), cex = LEGEND_TICK_CEX)

  segments(0.95, 6.60, 1.85, 6.60, col = "black", lty = 1, lwd = 2.5)
  text(2.08, 6.60, "q < 0.05", adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)
  segments(0.95, 6.13, 1.85, 6.13, col = "grey65", lty = 2, lwd = 2.5)
  text(2.08, 6.13, "q > 0.05", adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)

  text(bar_center_x, 5.55, "Heritability", adj = c(0.5, 0.5), font = 2, cex = LEGEND_TITLE_CEX)
  radius_in <- legend_circle_radii_in()
  y <- c(4.83, 3.95 + 0.37, 2.79 + 0.8)
  symbols(x = rep(1.88, 3), y = y, circles = radius_in,
          inches = FALSE, bg = "black", fg = NA, add = TRUE)
  text(2.53, y, labels = names(NODE_LEGEND_DIAM_MM), adj = c(0, 0.5), cex = LEGEND_TEXT_CEX)
}

plot_composite <- function() {
  par(fig = c(0, 10 / 15, 0, 1))
  state <- measure_network()
  par(fig = c(10 / 15, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  draw_legend(state)
  invisible(state)
}

png(file.path(OUTDIR, "network_GNA_network_only.png"), width = 10, height = 10, units = "in", res = 300)
network_png_state <- measure_network(); dev.off()
pdf(file.path(OUTDIR, "network_GNA_network_only.pdf"), width = 10, height = 10)
network_pdf_state <- measure_network(); dev.off()
png(file.path(OUTDIR, "network_GNA_legend_only.png"), width = 5, height = 10, units = "in", res = 300)
draw_legend(network_png_state); dev.off()
pdf(file.path(OUTDIR, "network_GNA_legend_only.pdf"), width = 5, height = 10)
draw_legend(network_pdf_state); dev.off()
png(file.path(OUTDIR, "network_GNA.png"), width = 15, height = 10, units = "in", res = 300)
combined_png_state <- plot_composite(); dev.off()
pdf(file.path(OUTDIR, "network_GNA.pdf"), width = 15, height = 10)
combined_pdf_state <- plot_composite(); dev.off()
if (!isTRUE(all.equal(network_png_state$pin, combined_png_state$pin, tolerance = 1e-8)) ||
    !isTRUE(all.equal(network_png_state$usr, combined_png_state$usr, tolerance = 1e-8)) ||
    !isTRUE(all.equal(network_pdf_state$pin, combined_pdf_state$pin, tolerance = 1e-8)) ||
    !isTRUE(all.equal(network_pdf_state$usr, combined_pdf_state$usr, tolerance = 1e-8))) {
  warning("独立网络与拼接网络的绘图区参数不一致，请检查图形设备或 R 图形参数。")
}
