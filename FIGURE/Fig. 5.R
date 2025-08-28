setwd("D:/pan_cancer/FIGURE 3")
library(ComplexHeatmap) 
library(circlize)
library(ggsci)
library(grid)          # packLegend() 需用到
old_opt <- ht_opt(HEATMAP_LEGEND_PADDING = unit(10, "mm"))  # 默认约 4 mm
on.exit(ht_opt(old_opt), add = TRUE)  # 画完恢复默认

dat  <- read.delim("heapmap.txt", check.names = FALSE)
mat  <- as.matrix(dat[, c("CF","EF","F1","F2","F3")])
rownames(mat) <- dat$Trait
split <- factor(dat$Category, levels = c("Input set","UKB_Finngen_meta","EAS"))

## 颜色 -----------------------------------------------------------------------
col_fun    <- colorRamp2(c(0, 0.4, 1), c("#0571b0", "white", "#ca0020"))
factor_col <- setNames(pal_npg("nrc")(5), c("CF","EF","F1","F2","F3"))

## 顶端列注（这里关闭自带图例，稍后手动打包） -------------------------------
top <- HeatmapAnnotation(
  f  = factor(colnames(mat), levels = names(factor_col)),
  col = list(f = factor_col),
  show_annotation_name = FALSE,
  show_legend = FALSE      ## <— 不要自动生成图注
)

## 行分块注 --------------------------------------------------------------------
row_anno <- rowAnnotation(
  block = anno_block(
    gp = gpar(fill = "#636363"),
    labels = levels(split),
    labels_gp = gpar(col = "white", fontsize = 7)  ## <— 分类文字大小
  ),
  show_annotation_name = FALSE,
  width = unit(8, "mm")                           ## <— 分类色块宽度
)
## 主热图 ----------------------------------------------------------------------
ht <- Heatmap(
  mat,
  name = "Genetic Correlation",
  col  = col_fun,
  rect_gp = gpar(col = "black", lwd = 1.2),
  width = unit(45, "mm"),
  top_annotation  = top,
  left_annotation = row_anno,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  row_split       = split,
  row_gap         = unit(4, "mm"),
  height = unit(nrow(mat)*7, "mm"),
  column_split    = factor(colnames(mat), levels = colnames(mat)),
  column_gap      = unit(2, "mm"),
  show_column_names = FALSE,
  column_title_side = "top",
  column_title_rot  = 0,
  column_title_gp   = gpar(fontsize = 9),
  row_names_gp      = gpar(fontsize = 6),
  row_title = NULL,
  column_names_side = "top",
  row_names_side    = "left",
  heatmap_legend_param = list(title_position = "topcenter",
                              direction = "vertical"),
  show_heatmap_legend = FALSE,               ## <— 关闭自动图例，避免重复
  cell_fun = function(j, i, x, y, w, h, fill)
    grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
)

## ---------------- 手动打包两个图例 ----------------
lgd_corr <- Legend(
  col_fun = col_fun,
  title   = "Genetic Correlation",
  title_gp = gpar(fontsize = 8),             ## <— 缩小标题字号
  title_position = "topcenter",
  direction = "vertical"
)

lgd_factor <- Legend(
  labels = names(factor_col),
  legend_gp = gpar(fill = factor_col),
  title  = "Carcinomagenesis\nLatent Factor",## <— 分两行并居中
  title_gp = gpar(fontsize = 8),             ## <— 缩小标题字号
  title_position = "topcenter",
  direction = "vertical"
)

legend_all <- packLegend(
  lgd_corr, lgd_factor,
  direction = "vertical",
  gap = unit(4, "mm")
)
pdf("heatmap_replication.pdf", width = 8, height = 11)

draw(
  ht,
  annotation_legend_list = NULL,
  heatmap_legend_list    = list(legend_all),
  heatmap_legend_side    = "right"
)
## 输出 ------------------------------------------------------------------------
dev.off()


