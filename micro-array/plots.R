setwd('/Users/ridwan204/Desktop/micro-array/gse149507/')

library(dplyr)
library(pheatmap)

edata <- read.table('./normalized_expression.tsv', sep='\t', header = TRUE)

edata <- edata %>% tibble::as_tibble() %>% tibble::column_to_rownames('Gene')


# annotation sample
ann_df <- data.frame(condition=if_else(grepl("^t\\w+", colnames(edata)), "Tumor", "Healthy"),
                    row.names = colnames(edata))
ann_color <- list(
  condition = c("Tumor"="red", 'Healthy'='green')
)

library(ggplot2)

pheatmap(edata, main="GSE149507", show_rownames=FALSE,
         #clustering_method = 'ward.D',
         border_color = F, fontsize_col = 8,
         annotation_names_col = F,
         annotation_col = ann_df, annotation_colors = ann_color
         
         )

# ggsave("heatmap_plot_all.png", heatmap_plot, dpi = 1200)

# hierarcical clustering
dist1 <- edata %>% t() %>% dist()
hclust1 = hclust(dist1)
plot(hclust1, hang = -1)

library(dendextend)
dend = as.dendrogram(hclust1)
# dend = color_labels(hclust1, 2, 1:2)
labels_colors(dend) = c(2, rep(3,6), 2, rep(3,6), rep(2,10))
hclustplot <- plot(dend, main='GSE149507')

# ggsave("hclust.png", hclustplot, dpi = 1200)


# kmeans clustering
k1 = kmeans(edata, centers=3)
table(k1$cluster)
pheatmap(edata[order(k1$cluster),], cluster_cols = F, cluster_rows=F,
         show_rownames = F)
  


