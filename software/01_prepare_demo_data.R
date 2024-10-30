## prepare data to test new function
ge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')

ge2 = ge[1:1e4, ]
me2 = me[1:200, ]
saveRDS(ge2, '/home/whou10/data/whou/metpred/software/demo_expr.rds')
saveRDS(me2, '/home/whou10/data/whou/metpred/software/demo_dnam.rds')
