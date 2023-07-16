rdir = '/home/whou10/data/whou/metpred/data/pbmc/res/'
w.tr = readRDS(paste0(rdir, 'wgbs_train.rds'))
source('/home/whou10/scratch16/whou10/resource/startup.R')
rs = rowsds(w.tr)
m = apply(w.tr,1,max) - apply(w.tr, 1, min)
selcpg = names(sort(rs[abs(m) > 0.1], decreasing = F)[1:1e4])
saveRDS(selcpg, '/home/whou10/data/whou/metpred/evaluate/pbmc/res/selcpg/top1e4sd.rds')


