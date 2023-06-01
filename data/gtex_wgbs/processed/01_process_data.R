w <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/wgbs/proc.rds')
r <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/rnaseq/pb.rds')
w <- w[,colnames(r)]
wbrain <- rowMeans(w[,grep('Brain_',colnames(w))])
w <- cbind(Brain=wbrain,w[,!grepl('Brain_',colnames(w))])

rbrain <- rowMeans(r[,grep('Brain_',colnames(r))])
r <- cbind(Brain=rbrain,r[,!grepl('Brain_',colnames(r))])
identical(colnames(w),colnames(r))

saveRDS(w,file='/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/wgbs.rds')
saveRDS(r,file='/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/rna.rds')
