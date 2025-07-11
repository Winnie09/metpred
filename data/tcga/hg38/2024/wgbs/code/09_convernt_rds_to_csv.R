## tcga wgbs and matched rna samples
## convert from rds to csv object
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38')
suppressMessages(library(GenomicRanges))

## read in wgbs
wgbs <- readRDS(paste0('wgbs/combine/me_cpg_by_sample.rds'))  
wgbscpg = rownames(wgbs)
write.csv(wgbs, 'wgbs/combine/me_cpg_by_sample.csv') 

ge = readRDS('wgbs/combine/ge.rds')
write.csv(ge, 'wgbs/combine/ge.csv')

