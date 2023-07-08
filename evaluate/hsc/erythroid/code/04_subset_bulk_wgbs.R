w = readRDS('/home/whou10/data/whou/metpred/data/hsc/proc/wgbs.rds')
ct = sapply(colnames(w), function(i) strsplit(i, '_')[[1]][2])
ct2 =  c(ct[grep('HSC', ct)], ct[grep('MPP', ct)], ct[grep('CMP', ct)], ct[grep('MEP', ct)])
w2 = w[, names(ct2)]
saveRDS(w2, '/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/bulk/bulk_wgbs.rds')

pred2 = readRDS(file='/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/res/pred_rd1e5.rds')
w3 = w2[rownames(pred2), ]
saveRDS(w3, '/home/whou10/data/whou/metpred/evaluate/hsc/erythroid/bulk/bulk_wgbs.rds')
