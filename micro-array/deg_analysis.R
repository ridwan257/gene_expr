
setwd('/Users/ridwan204/Desktop/micro-array/gse149507/')


library(tibble)
library(dplyr)
library(pheatmap)
source('~/my_lib/bioR/qq_plot.R')


edata <- read.table(file='./normalized_expression.tsv', sep='\t', header = T) %>% as_tibble()
top_table <- read.table(file='./tumor_vs_normal.tsv', sep='\t', header = T)
query_genes <- read.table("../query_set/all_genes.txt", sep = '\t', header = TRUE)

query_genes <- list(
  "metabolism" = query_genes$xenobiotics_metabolism[query_genes$xenobiotics_metabolism != ""],
  "reactive" = query_genes$oxidative_stress[query_genes$oxidative_stress != ""],
  "transporter" = query_genes$xenobiotics_transporter[query_genes$xenobiotics_transporter != ""]
)


genes <- top_table %>% 
  filter(Gene %in% query_genes$metabolism, abs(logFC)<1)


data <- edata %>%
  filter(Gene %in% genes$Gene) %>% tibble::column_to_rownames("Gene")

pheatmap(data, main = "Xenibiotic Metabolism", fontsize_row = 5,
         cluster_cols = F)
  