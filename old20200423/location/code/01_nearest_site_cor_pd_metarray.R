setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
library(GenomicRanges)
library(matrixStats)
set.seed(12345)
mat = readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
gr = readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
chrn <- as.character(commandArgs(trailingOnly = T)[[1]])
id = which(as.character(seqnames(gr))==chrn)
mat = mat[id,]
gr = gr[id] ###
saveRDS(mat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/data/mearray/beta_',chrn,'.rds'))
saveRDS(gr, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/data/mearray/basegr_',chrn,'.rds'))
        
s = start(gr)
identical(s, sort(s))
id <- order(s)
mat = mat[id,]
s <- s[id]
id <- which(rowSds(mat) > 0.1)
mat <- mat[id,]
s <- s[id]
sid <- which(s[-1]-s[-length(s)] < 100)
corv <- sapply(sid,function(i) {  
  cor(mat[i,], mat[i+1,], method='spearman')
})
summary(unlist(corv))

randomcorv <- sapply(1:length(sid),function(i) {  
  tid <- sample(1:length(s),2)
  #if (abs(s[tid[1]]-s[tid[2]])>100)
  cor(mat[tid[1],], mat[tid[2],], method='spearman')
})
summary(unlist(randomcorv))

df = data.frame(correlation=c(unlist(corv),unlist(randomcorv)), type=as.factor(c(rep('two sites(dist<100bp)',length(unlist(corv))), rep('two random sites',length(randomcorv)))))
saveRDS(df,paste0('./metpred/location/plot/array/nearest_site_cor_',chrn, '_pd.rds'))
