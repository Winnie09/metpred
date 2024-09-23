source('/home/whou10/data/whou/metpred/software/trainpredict.R')
source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'

jid <- as.numeric(commandArgs(trailingOnly = T))
print(jid)
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
rgroup <- as.numeric(cut(1:nrow(w.tr),20))
w.tr = w.tr[rgroup==jid,]

r.tr = readRDS(paste0(ddir, 'rna_train.rds'))
r.te = readRDS(paste0(ddir, 'rna_test_sc.rds'))

n <- colnames(r.tr)
r.tr <- r.tr[,sample(1:ncol(r.tr))]
colnames(r.tr) <- n


ct = sub(':.*', '', colnames(r.te))


r.te = r.te[, ct %in% c('b_cells', 'cd14_monocytes', 'cd56_nk', 
                        'cytotoxic_t', 'memory_t', 'naive_t',
                        'regulatory_t') ]
ct = sub(':.*', '', colnames(r.te))
id <- unlist(sapply(unique(ct), function(uct){
  which(ct == uct)[1:300]
}, simplify = F))
r.te = r.te[,id]

pred <- trainpredict(trainexpr=r.tr,
                     testexpr=r.te,
                     trainmeth=w.tr,
                     clunumlist = c(1e3),
                     lambdalist = c(10e-1))
saveRDS(pred, paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/competing/res/permu_',jid,'.rds'))



