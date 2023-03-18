source('/hpc/group/jilab/zj/met/software/trainmodel.R')
source('/hpc/group/jilab/zj/met/software/predict.R')
d <- readRDS('/hpc/group/jilab/zj/met/final/procrna.rds')
m <- readRDS('/hpc/group/jilab/zj/met/final/wgbs_hg38.rds')

tid <- which(sub('_.*','',colnames(d))=='GSE121723')
d <- d[,tid]
m <- m[,tid]

nav <- rowMeans(is.na(m))
gid <- names(which(nav==0))
m <- m[gid,]


gd <- readRDS('/hpc/group/jilab/zj/met/tcga/hg38/combine/ge.rds')
gm <- readRDS('/hpc/group/jilab/zj/met/tcga/hg38/combine/me.rds')
pro <- readRDS('/hpc/group/jilab/zj/met/tcga/hg38/combine/project.rds')
pro <- pro[grep('GBM',pro)]
gd <- gd[,names(pro)]
gm <- gm[,names(pro)]
suppressMessages(library(GenomicRanges))
grgr <- readRDS('/hpc/group/jilab/zj/met/tcga/hg38/gr/gr.rds')
gr <- paste0(as.character(seqnames(grgr)),'_',start(grgr))
names(gr) <- names(grgr)
rownames(gm) <- gr[rownames(gm)]
gm <- gm[intersect(rownames(gm),rownames(m)),]
gm <- gm[rowSums(is.na(gm))==0,]
load('/hpc/group/jilab/zj/met/tcga/hg38/gn/grch38.rda')
geneid <- sub('\\..*','',geneid)
gd <- gd[rownames(gd) %in% geneid,]
gn <- genename[match(rownames(gd),geneid)]
gd <- gd[!duplicated(gn),]
rownames(gd) <- gn[!duplicated(gn)]


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
  names(head(sort(csd,decreasing = T),10000))  
}

set.seed(12345)
hsdid <- sfunc(gm)
mod <- trainmodel(d,m[rownames(gm),])
pred <- predict(gd,mod)
true <- gm
sampcv <- corfunc(pred[hsdid,],true[hsdid,])
citecv <- corfunc(t(pred),t(true))
#sampcv <- corfunc(t(apply(pred,1,frank)),t(apply(true,1,frank)))
#citecv <- corfunc(apply(pred,1,frank),apply(true,1,frank))
summary(sampcv)
summary(citecv)

pd <- data.frame(cv=c(sampcv,citecv),type=rep(c('Across samples','Across CpG cites'),c(length(sampcv),length(citecv))),stringsAsFactors = F)

library(ggplot2)
pdf('/hpc/group/jilab/zj/met/res/plot/tcga_gbm/acrosssamplecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across samples',],aes(cv,type,col=type)) + geom_violin(col='royalblue') + geom_boxplot(col='royalblue',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation')
dev.off()

pdf('/hpc/group/jilab/zj/met/res/plot/tcga_gbm/acrosscitecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across CpG cites',],aes(cv,type,col=type)) + geom_violin(col='orange') + geom_boxplot(col='orange',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation') + scale_color_manual(values=c('Across CpG cites'='orange','Across samples'='royalblue'))
dev.off()


id='chr1_201731358'
pdf('/hpc/group/jilab/zj/met/res/plot/tcga_gbm/acrosssampleexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[id,],y=true[id,]),aes(x=x,y=y)) + geom_point(col='royalblue') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()

id <- which.max(citecv)
pdf('/hpc/group/jilab/zj/met/res/plot/tcga_gbm/acrossciteexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[,id],y=true[,id]),aes(x=x,y=y)) + geom_point(alpha=0.1,size=0.1,col='orange') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none') + geom_abline(slope=1,intercept = 0,col='red',linetype=2)
dev.off()

