source('/hpc/group/jilab/zj/met/software/trainmodel.R')
source('/hpc/group/jilab/zj/met/software/predict.R')
d <- readRDS('/hpc/group/jilab/zj/met/final/procrna.rds')
m <- readRDS('/hpc/group/jilab/zj/met/final/nonawgbs_hg38.rds')

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
}

sfunc <- function(tmpm) {
  cm <- rowMeans(tmpm)
  csd <- sqrt((rowMeans(tmpm*tmpm) - cm^2) / (ncol(tmpm) - 1) * ncol(tmpm))
  names(head(sort(csd,decreasing = T),100000))  
}

set.seed(12345)
trainid <- sample(colnames(m),ncol(m)*0.8)
testid <- setdiff(colnames(m),trainid)
hsdid <- sfunc(m[,testid])
mod <- trainmodel(d[,trainid],m[,trainid])
pred <- predict(d[,testid],mod)
true <- m[,testid]
sampcv <- corfunc(pred[hsdid,],true[hsdid,])
citecv <- corfunc(t(pred),t(true))
summary(sampcv)
summary(citecv)

set.seed(12345)
mod <- trainmodel(d[,trainid],m[,sample(trainid)])
pred <- predict(d[,testid],mod)
true <- m[,testid]
persampcv <- corfunc(pred[hsdid,],true[hsdid,])
percitecv <- corfunc(t(pred),t(true))
summary(persampcv)
summary(percitecv)


pd <- data.frame(cv=c(sampcv,citecv,persampcv,percitecv),type=rep(c('across samples','across CpG sites','permutation\nacross samples','permutation\nacross CpG sites'),c(length(sampcv),length(citecv),length(persampcv),length(percitecv))),stringsAsFactors = F)

library(ggplot2)
pdf('/hpc/group/jilab/zj/met/res/plot/bulkcv/acrosssamplecor.pdf',width=2,height=3)
ggplot(pd[pd$type%in%c('across samples','permutation\nacross samples'),],aes(cv,type,col=type)) + geom_violin() + scale_color_manual(values=c('across samples'='royalblue','permutation\nacross samples'='grey')) + geom_boxplot(width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation')
dev.off()

pdf('/hpc/group/jilab/zj/met/res/plot/bulkcv/acrosscitecor.pdf',width=2,height=3)
ggplot(pd[pd$type%in%c('across CpG sites','permutation\nacross CpG sites'),],aes(cv,type,col=type)) + geom_violin() + scale_color_manual(values=c('across CpG sites'='orange','permutation\nacross CpG sites'='grey')) + geom_boxplot(width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation')
dev.off()


id <- which.max(sampcv)
pdf('/hpc/group/jilab/zj/met/res/plot/bulkcv/acrosssampleexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[id,],y=true[id,]),aes(x=x,y=y)) + geom_point(col='royalblue') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()

id <- which.max(citecv)
samppid <- sample(1:nrow(pred),200000)
pdf('/hpc/group/jilab/zj/met/res/plot/bulkcv/acrossciteexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[samppid,id],y=true[samppid,id]),aes(x=x,y=y)) + geom_point(alpha=0.1,size=0.1,col='orange') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()


# samp <- sub('_.*','',colnames(d))
# tab <- table(samp)
# acrosscv <- sapply(names(tab)[tab >= 10],function(s) {
#   print(s)
#   trainid <- colnames(d)[which(samp!=s)]
#   testid <- setdiff(colnames(d),trainid)
#   hsdid <- sfunc(m[,testid])
#   clud <- clugene(d)
#   coef <- lmcoef(t(clud[,trainid]),t(m[hsdid,trainid]))
#   pred <- t(t(clud[,testid]) %*% coef)
#   true <- m[hsdid,testid]
#   cv <- corfunc(pred,true)
# })
# 
# withincv <- sapply(names(tab)[tab >= 10],function(s) {
#   print(s)
#   sl <- colnames(d)[which(samp!=s)]
#   trainid <- sample(sl,length(sl)*0.8)
#   testid <- setdiff(sl,trainid)
#   hsdid <- sfunc(m[,testid])
#   clud <- clugene(d[,sl])
#   coef <- lmcoef(t(clud[,trainid]),t(m[hsdid,trainid]))
#   pred <- t(t(clud[,testid]) %*% coef)
#   true <- m[hsdid,testid]
#   cv <- corfunc(pred,true)
# })
# 





