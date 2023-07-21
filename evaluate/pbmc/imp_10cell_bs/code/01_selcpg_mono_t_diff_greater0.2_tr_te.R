rdir = '/home/whou10/data/whou/metpred/data/pbmc/imp_10cell_bs/res/'
w.tr = readRDS(paste0(rdir, 'wgbs_train.rds'))
w.te = readRDS(paste0(rdir, 'wgbs_test.rds'))
colnames(w.te)
colnames(w.tr)
w.tr2 = w.tr[, grepl('alpha-beta_T_cell', colnames(w.tr)) | grepl('classical_monocyte', colnames(w.tr))]
w.te2 = w.te[, grepl('alpha-beta_T_cell', colnames(w.te)) | grepl('classical_monocyte', colnames(w.te))]

source('/home/whou10/scratch16/whou10/resource/startup.R')
w.te.agg = aggregatefunc2(d = w.te2, 
                          by = c(rep('monocyte', 2), rep('T_cell', 4)), 
                          fun = 'mean')
w.tr.agg = aggregatefunc2(d = w.tr2, 
                          by = c(rep('monocyte', 2), rep('T_cell', 4)), 
                          fun = 'mean')
diff.tr = w.tr.agg[,2] - w.tr.agg[,1]
diff.te = w.te.agg[,2] - w.te.agg[,1]
diff.tr2 = diff.tr[abs(diff.tr) > 0.2]
diff.te2 = diff.te[abs(diff.te) > 0.2]
int = intersect(names(diff.tr2), names(diff.te2))
set.seed(12345)
int2 = sample(int, 1e4)
cor(diff.tr[int2], diff.te[int2])
saveRDS(int2, '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell_bs/res/selcpg/top1e4_mono_t_diff_greater0.2_tr_te.rds')

