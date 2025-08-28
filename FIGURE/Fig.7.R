library(dplyr)
library(data.table)
library(geni.plots)

######
A<-fread("G:/MR/EUR_MR/PhWas_phenotypes_e_factor_plot.txt")
min(A$pvalue)

A$sign <- ifelse(A$b > 0, 1, -1)
library(ggsci)
groups_vec <- sort(unique(A$group))
make_cols <- function(n){
  banks <- list(
    pal_npg("nrc")(10),
    pal_lancet()(9),
    pal_jco()(10),
    pal_nejm()(8),
    pal_aaas()(10),
    pal_d3("category20")(20),
    pal_tron("legacy")(7)
  )
  cols <- unlist(banks, use.names=FALSE)
  if(length(cols) < n) cols <- rep(cols, length.out=n)
  cols[seq_len(n)]
}
col_vec <- make_cols(length(groups_vec))
names(col_vec) <- groups_vec
p1 <- fig_phewas(
  data = A,
  groups = groups_vec,
  colours = col_vec,
  label_thresh = 6.06e-6,thresh = 6.06e-6,title_center = TRUE,trunc = 1e-30,title = "Phenome-wide mendelian randomisation in MVP — outcome: EF latent factor",
  axis_text_angle = -70,
  axis_text_size = 10,
  label_size = 3)
# p1
ggplot2::ggsave("phewas_EF.png",p1, width=8, height=5)

######
A<-fread("G:/MR/EUR_MR/PhWas_phenotypes_common_factor_plot.txt")
min(A$pvalue)

A$sign <- ifelse(A$b > 0, 1, -1)
library(ggsci)
groups_vec <- sort(unique(A$group))
make_cols <- function(n){
  banks <- list(
    pal_npg("nrc")(10),
    pal_lancet()(9),
    pal_jco()(10),
    pal_nejm()(8),
    pal_aaas()(10),
    pal_d3("category20")(20),
    pal_tron("legacy")(7)
  )
  cols <- unlist(banks, use.names=FALSE)
  if(length(cols) < n) cols <- rep(cols, length.out=n)
  cols[seq_len(n)]
}
col_vec <- make_cols(length(groups_vec))
names(col_vec) <- groups_vec
p2 <- fig_phewas(
  data = A,
  groups = groups_vec,
  colours = col_vec,
  label_thresh = 6.06e-6,thresh = 6.06e-6,title_center = TRUE,trunc = 1e-30,title = "Phenome-wide mendelian randomisation in MVP — outcome: CF latent factor",
  axis_text_angle = -70,
  axis_text_size = 10,
  label_size = 3)
# p2
ggplot2::ggsave("phewas_CF.png",p2, width=8, height=5)

######
A<-fread("G:/MR/EUR_MR/PhWas_phenotypes_F1_plot.txt")
min(A$pvalue)

A$sign <- ifelse(A$b > 0, 1, -1)
library(ggsci)
groups_vec <- sort(unique(A$group))
make_cols <- function(n){
  banks <- list(
    pal_npg("nrc")(10),
    pal_lancet()(9),
    pal_jco()(10),
    pal_nejm()(8),
    pal_aaas()(10),
    pal_d3("category20")(20),
    pal_tron("legacy")(7)
  )
  cols <- unlist(banks, use.names=FALSE)
  if(length(cols) < n) cols <- rep(cols, length.out=n)
  cols[seq_len(n)]
}
col_vec <- make_cols(length(groups_vec))
names(col_vec) <- groups_vec
p3 <- fig_phewas(
  data = A,
  groups = groups_vec,thresh = 6.06e-6,title_center = TRUE,trunc = 1e-30,title = "Phenome-wide mendelian randomisation in MVP — outcome: F1 latent factor",
  colours = col_vec,
  label_thresh = 6.06e-6,
  axis_text_angle = -70,
  axis_text_size = 10,
  label_size = 3)
# p3
ggplot2::ggsave("phewas_F1.png",p3, width=8, height=5)

######
A<-fread("G:/MR/EUR_MR/PhWas_phenotypes_F2_plot.txt")
min(A$pvalue)

A$sign <- ifelse(A$b > 0, 1, -1)
library(ggsci)
groups_vec <- sort(unique(A$group))
make_cols <- function(n){
  banks <- list(
    pal_npg("nrc")(10),
    pal_lancet()(9),
    pal_jco()(10),
    pal_nejm()(8),
    pal_aaas()(10),
    pal_d3("category20")(20),
    pal_tron("legacy")(7)
  )
  cols <- unlist(banks, use.names=FALSE)
  if(length(cols) < n) cols <- rep(cols, length.out=n)
  cols[seq_len(n)]
}
col_vec <- make_cols(length(groups_vec))
names(col_vec) <- groups_vec
p4 <- fig_phewas(
  data = A,
  groups = groups_vec,thresh = 6.06e-6,title_center = TRUE,trunc = 1e-30,title = "Phenome-wide mendelian randomisation in MVP — outcome: F2 latent factor",
  colours = col_vec,
  label_thresh = 6.06e-6,
  axis_text_angle = -70,
  axis_text_size = 10,
  label_size = 3)
# p4
ggplot2::ggsave("phewas_F2.png",p4, width=8, height=5)

######
A<-fread("G:/MR/EUR_MR/PhWas_phenotypes_F3_plot.txt")
min(A$pvalue)
A$sign <- ifelse(A$b > 0, 1, -1)
library(ggsci)
groups_vec <- sort(unique(A$group))
make_cols <- function(n){
  banks <- list(
    pal_npg("nrc")(10),
    pal_lancet()(9),
    pal_jco()(10),
    pal_nejm()(8),
    pal_aaas()(10),
    pal_d3("category20")(20),
    pal_tron("legacy")(7)
  )
  cols <- unlist(banks, use.names=FALSE)
  if(length(cols) < n) cols <- rep(cols, length.out=n)
  cols[seq_len(n)]
}
col_vec <- make_cols(length(groups_vec))
names(col_vec) <- groups_vec
p5 <- fig_phewas(
  data = A,
  groups = groups_vec,
  colours = col_vec,thresh = 6.06e-6,
  label_thresh = 6.06e-6,title_center = TRUE,trunc = 1e-50,title = "Phenome-wide mendelian randomisation in MVP — outcome: F3 latent factor",
  axis_text_angle = -70,
  axis_text_size = 10,
  label_size = 3)
# p5
ggplot2::ggsave("phewas_F3.png",p5, width=8, height=5)


gc()



