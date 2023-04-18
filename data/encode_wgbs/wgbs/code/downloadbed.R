k <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/wgbs.rds')
library(data.table)
d <- fread('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/wgbs/metadata.tsv',data.table=F)
d <- d[d[,'Audit ERROR']!='extremely low coverage',]
d <- d[d[,'Experiment accession'] %in% k[,'Accession'],]
#identical(sort(d[,'Experiment accession']),sort(k[,'Accession']))
writeLines(paste0('wget ',d[,'File download URL'],' -P /home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/bed'),'/home/whou10/data/whou/metpred/data/encode_wgbs/wgbs/code/downloadbed.sh')

