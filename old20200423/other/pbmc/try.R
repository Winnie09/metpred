cellnum <- 1
m <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/pbmc/data/proc/matrix/saver.rds')
ct <- readRDS('/home-4/zji4@jhu.edu/scratch/raisin/pbmc/data/proc/ct/sc.rds')
m <- sapply(c('B_cell','natural_killer_cell','CD14_monocyte','CD8_T_cell'),function(sct) {
  rowMeans(m[,sample(which(ct==sct),cellnum),drop=F])
})

d <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38mat/mat.rds')
load('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/metadata/RNAseq/processed/processed.rda')
hg19_experiment[,'Biosample summary'] <- gsub(' ','_',hg19_experiment[,'Biosample summary'])
b <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/hg19/RNAseq/matrix/FPKM.rds')
colnames(b) <- sub('.tsv','',colnames(b))
b <- log2(b + 1)
b <- sapply(colnames(d),function(ub) {
  af <- hg19_file[hg19_file[,'Experiment accession'] %in% hg19_experiment[hg19_experiment[,'Biosample summary']==ub,'Accession'],1]
  rowMeans(b[,af,drop=F])
})

testid <- c('B_cell_male_adult_(37_years)','natural_killer_cell_male_adult_(37_years)','CD14-positive_monocyte_male_adult_(37_years)','T-cell_male_adult_(37_years)')
trainid <- setdiff(colnames(b),testid)
testb <- b[,testid]
testd <- d[,testid]
b <- b[,trainid]
d <- d[,trainid]

load('/home-4/zji4@jhu.edu/scratch/resource/gn/res/hg19.rda')
m <- m[row.names(m) %in% genename,]
row.names(m) <- geneid[match(row.names(m),genename)]
intg <- intersect(row.names(m),row.names(b))
m <- m[intg,]
b <- b[intg,]

library(preprocessCore)
dn <- dimnames(b)
b <- normalize.quantiles(b)
dimnames(b) <- dn
bm <- rowMeans(apply(b,2,sort))
rn <- row.names(m)
m <- apply(m,2,function(i) bm[rank(i)])
row.names(m) <- rn

dn <- dimnames(b)
sb <- t(apply(b,1,scale))
dimnames(sb) <- dn

clu <- kmeans(sb,nrow(b)/50)$cluster 
ecl <- t(sapply(1:max(clu),function(i) {
  colMeans(b[names(clu)[clu==i],])
}))
scecl <- t(sapply(1:max(clu),function(i) {
  colMeans(m[names(clu)[clu==i],])
}))
testecl <- t(sapply(1:max(clu),function(i) {
  colMeans(testb[names(clu)[clu==i],])
}))

library(matrixStats)
sdv <- rowSds(d)
id <- which(sdv > 30)

d <- d[id,]
d <- d/100
d[d==0] <- 0.5/100
d[d==100] <- 99.5/100
d <- log(d/(1-d))
# dmv <- rowMeans(d)
# dsdv <- apply(d,1,sd)
# d <- (d-dmv)/dsdv

# cid <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/metclu/clu1000.rds')
# dm <- t(sapply(1:max(cid),function(i) {
#   colMeans(d[cid==i,])
# }))

set.seed(12345)

pred <- d %*% (solve(t(ecl) %*% ecl) %*% t(ecl) %*% scecl)
#clupred <- dm %*% solve(t(ecl) %*% ecl) %*% t(ecl) %*% scecl
pred <- singlepred
#pred <- (clupred[cid,] + singlepred)/2
#pred <- pred*dsdv + dmv
pred <- exp(pred)
pred <- pred/(1+pred)

testd <- testd[id,]
sapply(1:ncol(pred),function(i) {
  cor(pred[,i],testd[,i])
})

summary(sapply(sample(1:nrow(pred),10000),function(i) {
  cor(pred[i,],testd[i,])
}))


i <- 2
pred <- d[i,] %*% (solve(t(ecl) %*% ecl) %*% t(ecl) %*% scecl)
