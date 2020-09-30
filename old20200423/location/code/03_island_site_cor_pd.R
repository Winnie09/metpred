setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
library(GenomicRanges)
library(matrixStats)
chrn = as.character(commandArgs(trailingOnly = T)[1])
set.seed(12345)
tb = read.table('./metpred/data/data/island/hg38/cpgIslandExt.txt',as.is=T)
tb[,5] = paste0(tb[,5],tb[,6])
tb = tb[,-6]
cn = read.table('./metpred/data/data/island/hg38/colname.txt',as.is=T)
colnames(tb) <- as.character(cn[1,])
tb = tb[tb$chrom==chrn,]
mat=readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/data/wgbs/beta_',chrn,'.rds'))
gr=readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/data/wgbs/basegr_',chrn,'.rds'))
if (!identical(gr, sort(gr))) gr <- sort(gr)
tb <- GRanges(seqnames=tb[,2],IRanges(start=tb[,3],end=tb[,4]))
o <- as.matrix(findOverlaps(gr,tb))
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
v <- rep(0,length(gr))
v[o[,1]] <- o[,2]
corv = sapply(1:100,function(i) { ## 1000
  set.seed(i)
  isd <- sample(1:length(tb),1)
  tmp = which(v==isd)
  tmp = sample(tmp, min(50, length(tmp)))
  a = corfunc(t(mat[tmp,]), t(mat[tmp,]))
  a[upper.tri(a)] = NA
  a = as.vector(a)
  a = a[!is.na(a)]
})
corv = unlist(corv)
summary(corv)

randomcorv = sapply(1:100,function(i) { ## 1000
  set.seed(i)
  isd <- sample(1:length(tb),2)
  tmp1 = which(v==isd[1])
  tmp2 = which(v==isd[2])
  tmp1 = sample(tmp1, min(50, length(tmp1)))
  tmp2 = sample(tmp2, min(50, length(tmp2)))
  a = corfunc(t(mat[tmp1,]), t(mat[tmp2,]))
  a[upper.tri(a)] = NA
  a = as.vector(a)
  a = a[!is.na(a)]
})
randomcorv = unlist(randomcorv)
summary(randomcorv)

df = data.frame(correlation=c(unlist(corv),unlist(randomcorv)), type=as.factor(c(rep('two sites(same island)',length(unlist(corv))), rep('two sites(different island)',length(randomcorv)))))
df$type = factor(df$type, levels = c('two sites(same island)', 'two sites(different island)'))
saveRDS(df,paste0('./metpred/location/plot/wgbs/island_site_cor_',chrn, '_pd.rds'))
