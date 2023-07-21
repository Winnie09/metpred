# Complete cases in w.tr. 
# Group CpGs according to their sd in w.tr.
source('/home/whou10/scratch16/whou10/resource/startup.R')
ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp_10cell/res/'
w.tr = readRDS(paste0(ddir, 'wgbs_train.rds'))
w.tr2 = w.tr[complete.cases(w.tr), ]
rs = rowsds(w.tr2)
rs2 <- cut(rs,seq(0,1,0.05))
names(rs2) <- names(rs)
dir.create(rdir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/imp_10cell/res/cpggroup/', recursive = T)
cpglist <- list()
for (rs2.s in unique(rs2)){
  print(rs2.s)
  cpglist[[sub('\\]', '', sub('\\(','',rs2.s))]] = names(rs2)[which(rs2 == rs2.s)]
}
cpglist = cpglist[sapply(cpglist, length) > 0]
saveRDS(cpglist, paste0(rdir, 'cpggroup_list.rds'))


