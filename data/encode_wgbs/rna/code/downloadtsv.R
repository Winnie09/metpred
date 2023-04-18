k <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/rna.rds')
library(data.table)
d <- fread('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/rna/metadata.tsv',data.table=F)
d <- d[grep('ENCODE4',d[,'File analysis title']),]
d <- d[nchar(d[,'Audit ERROR'])==0,]
d <- d[d[,'Experiment accession'] %in% k[,'Accession'],]
identical(unique(sort(d[,'Experiment accession'])),sort(k[,'Accession']))
writeLines(paste0('wget ',d[,'File download URL'],' -P /home/whou10/data/whou/metpred/data/encode_wgbs/rna/tsv'),'/home/whou10/data/whou/metpred/data/encode_wgbs/rna/code/downloadtsv.sh')


