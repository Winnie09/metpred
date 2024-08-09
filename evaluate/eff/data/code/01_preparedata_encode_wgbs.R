ddir <- '/home/whou10/data/whou/metpred/evaluate/eff/data/processed/'
ge <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')
num.sample = as.numeric(commandArgs(trailingOnly = T)[1])
print(num.sample)
set.seed(12345)
id = sample(ncol(ge), num.sample, replace = T)
ge2 = ge[,id]
me2 = me[,id]
dir.create(paste0(ddir,num.sample),showWarnings = F)
saveRDS(ge2, paste0(ddir, num.sample,'/rna.rds'))      
saveRDS(me2, paste0(ddir, num.sample,'/wgbs.rds'))      

