me <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/wgbs.rds')
source('/home/whou10/data/whou/metpred/software/hg38_bin.R')
me <- binfunc(me)
saveRDS(me,'/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/wgbs_bin.rds')


