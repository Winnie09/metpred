
setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
w <- readRDS('data/mousepancreas/wgbs/bs.rds')
alls = sub(':.*','',colnames(w))
s <- commandArgs(trailingOnly = T)[1]
if (s == 'A46_F_P_K4_D4') stop('exists!')
w.sub = w[, !alls %in% s]
rs <- rowsds(w.sub)
rs2 <- cut(rs,seq(0,1,0.05))
names(rs2) <- names(rs)
rdir <- 'evaluate/mousepancreas/celltypespecific_cv/cpg_group/'
cpglist <- list()
for (rs2.s in unique(rs2)){
  print(rs2.s)
  cpglist[[sub('\\]', '', sub('\\(','',rs2.s))]] = names(rs2)[which(rs2 == rs2.s)]
}
saveRDS(cpglist, paste0(rdir, s, '.rds'))
