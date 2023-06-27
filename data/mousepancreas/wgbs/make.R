load('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/data/adm_bsseq_smoothed_cov_filtered.rda')
library(bsseq)
gr <- bsseq_smoothed_sub@rowRanges
m <- getMeth(bsseq_smoothed_sub)
rownames(m) <- paste0(as.character(seqnames(gr)),':',start(gr))
m <- m[rowMeans(is.na(m))==0,]
colnames(m) <- sub('_ADM$',':adm',sub('_unrecovered_ADM',':adm',sub('_duct',':duct',sub('_acinar',':acinar',colnames(m)))))
saveRDS(m,file='/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')

chr <- sub(':.*', '', sub('_.*', '', rownames(m)))
m2 <- m[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y')), ]
saveRDS(m2,file='/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')

