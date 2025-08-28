term_volcano<-function (data, term_data, x = "log2FoldChange", y = "padj", 
    normal_point_color = "#999999", normal_point_size = 1, deg_point_color = "black", 
    deg_point_fill = c(dendritic = "#49c2c6", `ion transport.` = "#fbcbcc", 
        metabolic = "#eef0ac", myelin = "#b1daa7", synaptic = "#d0d0a0"), 
    deg_point_size = 2, legend_background_fill = "#fefde2", legend_title = NULL, 
    legend_position = "UL", add_line = TRUE, log2FC_cut = NULL, 
    FDR_cut = 0.05, add_label = TRUE, label = "row", label_number = 10, 
    custom_label = NULL, x_lab = NULL, y_lab = NULL, output = TRUE, 
    filename = "volcano_plot") 
{
  library(ggVolcano)
  library(latex2exp)
    colnames(term_data) <- c("geneName", "term")
    colnames(data)[colnames(data) == x] <- "x"
    colnames(data)[colnames(data) == y] <- "y"
    colnames(data)[colnames(data) == label] <- "geneName"
    data$GO_term <- "others"
    term_data <- term_data[term_data$geneName %in% data$geneName, 
        ]
    data[term_data$geneName, ]$GO_term <- term_data$term
    Down_num <- length(which(data$y < 0.05 & data$x < 0))
    Up_num <- length(which(data$y < 0.05 & data$x > 0))
    color <- rep(normal_point_color, nrow(data))
    if (is.null(custom_label)) {
        if (label_number != 0) {
            data$label <- rep("", nrow(data))
            data$label[order(data$y)[1:label_number]] <- data$geneName[order(data$y)[1:label_number]]
        }
        else {
            data$label <- rep("", nrow(data))
        }
    }
    else {
        data$label <- rep("", nrow(data))
        data$label[match(custom_label, data$geneName)] <- custom_label
    }
    p <- ggplot(data[which(data$GO_term != "others"), ], aes(x, 
        -log10(y), fill = GO_term)) + geom_point(data = data[which(data$GO_term == 
        "others"), ], aes(x, -log10(y)), size = normal_point_size, 
        color = normal_point_color) + geom_point(size = deg_point_size, 
        shape = 21, color = deg_point_color) + scale_fill_manual(values = deg_point_fill) + 
        labs(x = x_lab %||% TeX("rho"), y = y_lab %||% 
            TeX("$-Log_{10} \\textit{P} $")) + theme(title = element_text(size = 15), 
        text = element_text(size = 15)) + theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.background = element_rect(fill = legend_background_fill, 
            colour = "black", size = 0.2), legend.key = element_rect(fill = legend_background_fill), 
        legend.key.size = unit(12, "pt")) + guides(fill = guide_legend(title = legend_title %||% 
        "Cancer"))
    if (add_line == TRUE) {
        if (is.null(log2FC_cut)) {
            p <- p + geom_vline(xintercept = 0, linetype = "longdash") + 
                geom_hline(yintercept = -log10(FDR_cut), linetype = "longdash")
        }
        else {
            p <- p + geom_vline(xintercept = c(-log2FC_cut, log2FC_cut), 
                linetype = "dashed") + geom_hline(yintercept = -log10(FDR_cut), 
                linetype = "dashed")
        }
    }
    if (add_label == TRUE) {
        p <- p + geom_text_repel(data = data, aes(x, -log10(y), 
            label = label), size = 3, max.overlaps = 100, key_glyph = draw_key_point)
    }
    if (legend_position == "UL") {
        p <- p + theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 
            1), )
    }
    else if (legend_position == "UR") {
        p <- p + theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 
            1), )
    }
    else if (legend_position == "DL") {
        p <- p + theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 
            0), )
    }
    else {
        p <- p + theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 
            0), )
    }
    if (output == TRUE) {
        ggsave(paste0(filename, ".pdf"), plot = p, height = 5, 
            width = 6)
    }
    return(p)
}
library(ggsci)
library(scales)   
pal_npg_10 <- pal_npg("nrc")(10)
# show_col(pal_npg_10)
library(vroom)


setwd("D:/pan_cancer/FIGURE 3")
A<-vroom("PRCA_LAVA.txt")
B<-vroom("PRCA_LAVA_term.txt")
df       <- as.data.frame(A)
term_df  <- as.data.frame(B)
colnames(term_df) <- c("geneName", "term")
deg_point_fill = c(PRCA = "#E64B35FF",BRCA = "#eef0ac", LC = "#F39B7FFF", THCA ="#00A087FF", 
                   CRC = "#4DBBD5FF", RCC = "#91D1C2FF", EC = "#8491B4FF",
                   ESCA = "#3C5488FF", OC = "#7E6148FF")
rownames(df) <- df$Gene
pdf("volcano_PRCA_LAVA.pdf", width = 5, height = 4)

term_volcano(
  df, term_df,legend_title = "PRCA vs",
  x              = "rho",
  y              = "p",
  FDR_cut = 0.05/334,
  deg_point_fill = deg_point_fill,
  add_label = FALSE,
  label          = "Gene",
  label_number   = 10,
  output         = FALSE
)
dev.off()


setwd("D:/pan_cancer/FIGURE 3")
A<-vroom("BRCA_LAVA.txt")
B<-vroom("BRCA_LAVA_term.txt")
df       <- as.data.frame(A)
term_df  <- as.data.frame(B)
colnames(term_df) <- c("geneName", "term")
rownames(df) <- df$Gene
pdf("volcano_BRCA_LAVA.pdf", width = 5, height = 4)

term_volcano(
  df, term_df, term_df,legend_title = "BRCA vs",
  x              = "rho",
  y              = "p",
  FDR_cut = 0.05/334,
  deg_point_fill = deg_point_fill,
  add_label = FALSE,
  label          = "Gene",
  label_number   = 10,
  output         = FALSE
)
dev.off()

setwd("D:/pan_cancer/FIGURE 3")
A<-vroom("LUNG_LAVA.txt")
B<-vroom("LUNG_LAVA_term.txt")
df       <- as.data.frame(A)
term_df  <- as.data.frame(B)
colnames(term_df) <- c("geneName", "term")
rownames(df) <- df$Gene
pdf("volcano_LUNG_LAVA.pdf", width = 5, height = 4)

term_volcano(
  df, term_df, term_df,legend_title = "LC vs",
  x              = "rho",
  y              = "p",
  FDR_cut = 0.05/334,
  deg_point_fill = deg_point_fill,
  add_label = FALSE,
  label          = "Gene",
  label_number   = 10,
  output         = FALSE
)
dev.off()

