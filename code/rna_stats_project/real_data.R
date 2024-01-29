library(DESeq2)

read_csv = read.csv('/Users/savely/evry/statistics/project/data/Tutorial_COUNTS.csv',sep='\t')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("genefilter","DESeq2"),type="source")

while (11 > 2){i <- 0}