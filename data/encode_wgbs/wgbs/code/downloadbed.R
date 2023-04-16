k <- readRDS('/home/whou10/data/whou/metpred/data/encode/meta/combine/wgbs.rds')
library(data.table)
d <- fread('/home/whou10/data/whou/metpred/data/encode/meta/wgbs/metadata.tsv',data.table=F)
d <- d[d[,'File format']=='bed bed9+',]
d <- d[d[,'Experiment accession'] %in% k[,'Accession'],]
#identical(sort(d[,'Experiment accession']),sort(k[,'Accession']))
writeLines(paste0('wget ',d[,'File download URL'],' -P /home/whou10/data/whou/metpred/data/encode/wgbs/bed'),'/home/whou10/data/whou/metpred/data/encode/wgbs/code/downloadbed.sh')
