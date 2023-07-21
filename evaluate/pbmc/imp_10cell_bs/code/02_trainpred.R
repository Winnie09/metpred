ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp_10cell_bs/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))
selcpg = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell_bs/res/selcpg/top1e4_mono_t_diff_greater0.2_tr_te.rds')

source('/home/whou10/data/whou/metpred/software/trainpredict.R')

pred <- trainpredict(trainexpr=r.tr,testexpr=r.te,trainmeth=w.tr[selcpg,],
                     clunumlist = c(1e3),
                     lambdalist = c(10e-1))
saveRDS(pred, '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell_bs/res/pred/pred_top1e4_mono_t_diff_greater0.2_tr_te.rds')

