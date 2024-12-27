setwd('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/')
# trainarray_predepic.rds  
pdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/plot/'

rowsds <- function(data) {
  cm <- rowMeans(data)
  sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
}

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}

################ WGBS 
pred = readRDS('trainarray_predwgbs.rds')
me = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs.rds')
str(me)
int = intersect(rownames(pred), rownames(me))
pred = pred[int, ]
me = me[int, ]

c = corfunc(pred, me)
pdf(paste0(pdir, 'wgbs_acrosssample_pcc.pdf'), width = 5, height = 3)
hist(c, xlab = 'across-sample PCC', ylab = '#samples', main='')
dev.off()

c2 = corfunc(t(pred), t(me))
pdf(paste0(pdir, 'wgbs_acrosscpg_pcc.pdf'), width = 5, height = 3)
hist(c2, xlab='across-CpG PCC', ylab='#CpGs', main = '')
dev.off()


################ EPIC 
pred = readRDS('trainarray_predepic.rds')
str(pred)
me = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/EPIC.rds')
cm = colMeans(is.na(me))
pdf(paste0(pdir, 'epic_prop_NA.pdf'), width = 5, height = 3)
hist(cm, xlab='Proportion of NAs in each sample', main='', ylab='#EPIC samples')
dev.off()

id = which(colMeans(is.na(me)) > 0.05)
me = me[, -id]

str(me)
int = intersect(rownames(pred), rownames(me))
pred = pred[int, colnames(me)]
me = me[int, ]

c = corfunc(pred, me)
pdf(paste0(pdir, 'EPIC_acrosssample_pcc.pdf'), width = 5, height = 3)
hist(c, xlab = 'across-sample PCC', ylab = '#samples', main='')
dev.off()

c2 = corfunc(t(pred), t(me))
pdf(paste0(pdir, 'EPIC_acrosscpg_pcc.pdf'), width = 5, height = 3)
hist(c2, xlab='across-CpG PCC', ylab='#CpGs', main = '')
dev.off()

