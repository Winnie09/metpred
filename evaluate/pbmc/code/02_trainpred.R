ddir <- '/home/whou10/data/whou/metpred/data/pbmc/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
r.te = readRDS(paste0(ddir, 'rna_test.rds'))
selcpg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/res/selcpg/top1e4sd.rds')

source('/home/whou10/data/whou/metpred/software/trainpredict.R')

pred <- trainpredict(trainexpr=r.tr,testexpr=r.te,trainmeth=w.tr[selcpg,])
saveRDS(pred, '/home/whou10/data/whou/metpred/evaluate/pbmc/res/pred/pred_top1e4_cpg.rds')
