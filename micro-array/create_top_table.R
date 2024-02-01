
setwd('/Users/ridwan204/Desktop/micro-array/gse149507/')


library(tibble)
library(dplyr)
source('~/my_lib/bioR/qq_plot.R')

edata <- read.table(file='./normalized_expression.tsv', sep='\t', header = T)

# edata <- edata %>% 
#    tibble::column_to_rownames('Gene')

ncol(edata)

results <- data.frame(
  'Gene' = edata$Gene,
  't' = rep(0, nrow(edata)),
  'pval' = rep(0, nrow(edata)),
  'logFC' = rep(0, nrow(edata))
)



# logFC = tumor / normal
for (i in 1:nrow(edata)){
  values <- edata[i,2:ncol(edata)] %>% as.numeric()
  results$logFC[i] <- (mean(values[1:12]) - mean(values[13:24]))
  t_stat <- t.test(values[1:12], values[13:24])
  results$t[i] <- t_stat$statistic
  results$pval[i] <- t_stat$p.value
  # break
}

# write.table(results, file='./tumor_vs_normal.tsv', sep = '\t', row.names = F, quote = F)




















