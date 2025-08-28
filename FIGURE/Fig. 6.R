library(data.table)
library(dplyr)
library(circlize)

data<-fread("mutiomics.txt")
data<-as.data.frame(data)
rownames(data)<-data$symbol
data<-data[,-1]
data <- data %>% mutate(sum = rowSums(across(1:7, ~ .x == "Y"), na.rm = TRUE))

col_fun = list(
  col_1 = c("Y"="#374E55", "N"="#E7EAEB"),
  col_2 = c("Y"="#DF8F44", "N"="#FBF2E9"),
  col_3 = c("Y"="#00A1D5", "N"="#E0F4FA"),
  col_4 = c("Y"="#B24745", "N"="#F6E9E9"),
  col_5 = c("Y"="#79AF97", "N"="#EFF5F3"),
  col_6 = c("Y"="#6A6599", "N"="#EDEDF3"),
  col_7 = c("Y"="#80796B", "N"="#F0EFED")
)

{
pdf("plot.pdf", height = 10, width = 10)
circos.par$gap.degree <- 12
circos.par$start.degree <- 85
circos.par$track.margin <- c(0.001, 0.001)

data_tmp <- as.matrix(data[,8])
rownames(data_tmp) <- rownames(data)
circos.heatmap.initialize(data_tmp, cluster = F)

circos.track(
  ylim = c(0, 6),
  bg.border = NA,
  panel.fun = function(x, y){
    circos.barplot(data_tmp,
      border = NA,
      1:nrow(data_tmp) - 0.5,
      col = "#b7a085")
  })

for (i in 1:7) {
  data_tmp <- as.matrix(data[,i])
  if (i == 1) {
    rownames(data_tmp) <- rownames(data)
  }
  colnames(data_tmp) <- colnames(data)[i]
  circos.heatmap(data_tmp,
    col = col_fun[[i]],
    rownames.side = "outside",
    cluster = FALSE,
    cell.border = "white",
    track.height = 0.05)
}

circos.track(track.index = get.current.track.index(),
  panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 1) {
      circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
        CELL_META$cell.xlim[2] + convert_x(6, "mm"), 7,
        col = "#e5e8ed", border = NA)
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 0.5,
        "7", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 1.5,
        "6", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 2.5,
        "5", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 3.5,
        "4", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 4.5,
        "3", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 5.5,
        "2", cex = 0.5, facing = "inside")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 6.5,
        "1", cex = 0.5, facing = "inside")
    }
  },
  bg.border = NA)
circos.clear()
dev.off()
}
