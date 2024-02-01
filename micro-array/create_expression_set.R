
# setwd('/Users/ridwan204/Desktop/micro-array/gse149507')

library(GEOquery)
library(tibble)
library(dplyr)

gse <- getGEO("GSE149507")
length(gse)
gse <- gse[[1]]

pdata <- pData(gse)
fdata <- fData(gse)
edata <- exprs(gse)


# format the patient data 
meta_data <- pdata %>% 
  filter(`smoking_status:ch1` == 'Ever') %>% 
  select('smoking_status:ch1', 'tissue:ch1') %>% 
  tibble::rownames_to_column('P_ID') %>% 
  rename(smoker='smoking_status:ch1', tissue='tissue:ch1')

status <- c()
i <- j <- 1
for(value in grepl('^sm.*', meta_data$tissue)){
  if (value){
    status <- c(status, paste0("tumor", i))
    i <- i + 1
  } else {
    status <- c(status, paste0("normal", j))
    j <- j + 1
  }
}
meta_data$tissue <- status

# entrez gene id converter
gene_mapper <- read.table('../entrez_gene_annonation.tsv', 
                          sep='\t', header = FALSE, quote="",
                          col.names = c("ENTREZ_GENE_ID", "Symbol"))

# spot id to symbol convert
spotID_mapper <- fdata %>% 
  select('ID', 'ENTREZ_GENE_ID') %>% 
  inner_join(gene_mapper, by='ENTREZ_GENE_ID') %>% 
  rename('SPOT_ID'='ID') %>% select(SPOT_ID, 'Symbol')

# change the column name as meta_data condition
edata <- edata %>% as.data.frame() %>% 
  rownames_to_column("SPOT_ID") %>% as_tibble() %>% 
  select(SPOT_ID, meta_data$P_ID) %>% 
  rename_at(meta_data$P_ID, ~meta_data$tissue)

edata$SPOT_ID <- spotID_mapper$Symbol[match(edata$SPOT_ID, spotID_mapper$SPOT_ID)]

edata <- edata %>% rename(Gene = SPOT_ID) 
edata <- edata %>% select(Gene, names(edata)[seq(2,25,2)], everything())

edata %>% is.na() %>% colSums()
edata <- edata %>% na.omit()

# edata %>% write.table(file='normalized_expression.tsv', sep='\t', row.names = F)

