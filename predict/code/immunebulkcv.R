source('/hpc/group/jilab/zj/met/software/trainmodel.R')
source('/hpc/group/jilab/zj/met/software/predict.R')
d <- readRDS('/hpc/group/jilab/zj/met/final/procrna.rds')
m <- readRDS('/hpc/group/jilab/zj/met/final/wgbs_hg38.rds')

nav <- rowMeans(is.na(m))
gid <- names(which(nav==0))
m <- m[gid,]

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

ct <- readRDS('/hpc/group/jilab/zj/met/quality/ct.rds')

testid <- sample(names(which(ct=='WBC')),sum(ct=='WBC',na.rm=T)/2)
#testid <- names(which(ct=='WBC'))
trainid <- setdiff(colnames(m),testid)
hsdid <- sfunc(m[,testid])
mod <- trainmodel(d[,trainid],m[,trainid])
pred <- predict(d[,testid],mod)
true <- m[,testid]
sampcv <- corfunc(pred[hsdid,],true[hsdid,])
citecv <- corfunc(t(pred),t(true))
summary(sampcv)
summary(citecv)

pd <- data.frame(cv=c(sampcv,citecv),type=rep(c('Across samples','Across CpG cites'),c(length(sampcv),length(citecv))),stringsAsFactors = F)

library(ggplot2)
pdf('/hpc/group/jilab/zj/met/res/plot/immunebulkcv/acrosssamplecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across samples',],aes(cv,type,col=type)) + geom_violin(col='royalblue') + geom_boxplot(col='royalblue',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation')
dev.off()

pdf('/hpc/group/jilab/zj/met/res/plot/immunebulkcv/acrosscitecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across CpG cites',],aes(cv,type,col=type)) + geom_violin(col='orange') + geom_boxplot(col='orange',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation') + scale_color_manual(values=c('Across CpG cites'='orange','Across samples'='royalblue'))
dev.off()


id <- names(which.max(sampcv))
pdf('/hpc/group/jilab/zj/met/res/plot/immunebulkcv/acrosssampleexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[id,],y=true[id,]),aes(x=x,y=y)) + geom_point(col='royalblue') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()

id <- which.max(citecv)
samppid <- sample(1:nrow(pred),200000)
pdf('/hpc/group/jilab/zj/met/res/plot/immunebulkcv/acrossciteexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[samppid,id],y=true[samppid,id]),aes(x=x,y=y)) + geom_point(alpha=0.1,size=0.1,col='orange') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()


