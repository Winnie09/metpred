ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp_10cell/res/'
r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
w.tr2 = w.tr[complete.cases(w.tr), ]
cpglist = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/res/cpggroup/cpggroup_list.rds')

i = as.numeric(commandArgs(trailingOnly = T)[1])
print(i)

source('/home/whou10/data/whou/metpred/software/trainpredict.R')
pred <- trainpredict(trainexpr=r.tr,testexpr=r.te,trainmeth=w.tr[cpglist[[i]],],
                     clunumlist = c(1e3),
                     lambdalist = c(10e-1))
saveRDS(pred, paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/res/predall/', names(cpglist)[i], '.rds'))


