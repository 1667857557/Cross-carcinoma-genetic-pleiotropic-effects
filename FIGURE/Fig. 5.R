setwd("G:/pan_cancer/FIGURE 3")
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(grid)

old_opt <- ht_opt(HEATMAP_LEGEND_PADDING = unit(10, "mm"))
on.exit(ht_opt(old_opt), add = TRUE)

dat <- read.delim("heapmap.txt", check.names = FALSE)
mat <- as.matrix(dat[, c("CF", "EF", "F1", "F2", "F3")])
rownames(mat) <- dat$Trait
split <- factor(dat$Category, levels = c("Input set", "UKB_Finngen_meta", "EAS"))

col_fun <- colorRamp2(c(0, 0.4, 1), c("#0571b0", "white", "#ca0020"))
factor_col <- setNames(
  c("#996699", "#008B8B", "#B22222", "#4169E1", "#FF8C00"),
  c("CF", "EF", "F1", "F2", "F3")
)

top <- HeatmapAnnotation(
  f = factor(colnames(mat), levels = names(factor_col)),
  col = list(f = factor_col),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

row_anno <- rowAnnotation(
  block = anno_block(
    gp = gpar(fill = "#636363"),
    labels = levels(split),
    labels_gp = gpar(col = "white", fontsize = 7)
  ),
  show_annotation_name = FALSE,
  width = unit(8, "mm")
)

ht <- Heatmap(
  mat,
  name = "Genetic Correlation",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 1.2),
  width = unit(45, "mm"),
  top_annotation = top,
  left_annotation = row_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = split,
  row_gap = unit(4, "mm"),
  height = unit(nrow(mat) * 7, "mm"),
  column_split = factor(colnames(mat), levels = colnames(mat)),
  column_gap = unit(2, "mm"),
  show_column_names = FALSE,
  column_title_side = "top",
  column_title_rot = 0,
  column_title_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 6),
  row_title = NULL,
  column_names_side = "top",
  row_names_side = "left",
  heatmap_legend_param = list(
    title_position = "topcenter",
    direction = "vertical"
  ),
  show_heatmap_legend = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
  }
)

lgd_corr <- Legend(
  col_fun = col_fun,
  title = "Genetic Correlation",
  title_gp = gpar(fontsize = 8),
  title_position = "topcenter",
  direction = "vertical"
)

lgd_factor <- Legend(
  labels = names(factor_col),
  legend_gp = gpar(fill = factor_col),
  title = "Carcinomagenesis\nLatent Factor",
  title_gp = gpar(fontsize = 8),
  title_position = "topcenter",
  direction = "vertical"
)

legend_all <- packLegend(
  lgd_corr,
  lgd_factor,
  direction = "vertical",
  gap = unit(4, "mm")
)

pdf("heatmap_replication.pdf", width = 8, height = 11)

draw(
  ht,
  annotation_legend_list = NULL,
  heatmap_legend_list = list(legend_all),
  heatmap_legend_side = "right"
)

dev.off()
