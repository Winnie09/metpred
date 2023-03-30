ge <- readRDS('/home/whou10/data/whou/metpred/combine/rna/hg38_encode.rds')
me <- readRDS('/home/whou10/data/whou/metpred/combine/wgbs/hg38_encode.rds')
identical(colnames(ge), colnames(me))
ge2 <- ge[,1:20]
me2 <- me[,1:20]
me2 <- me2[1:1e3, ]
saveRDS(ge2, '/home/whou10/data/whou/metpred/software_test/rna_encode_sub_20.rds')
saveRDS(me2, '/home/whou10/data/whou/metpred/software_test/wgbs_encode_sub_20.rds')
