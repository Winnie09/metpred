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
source('/hpc/group/jilab/zj/met/software/trainmodel.R')
source('/hpc/group/jilab/zj/met/software/predict.R')
m <- readRDS('/hpc/group/jilab/zj/met/final/wgbs_hg38.rds')
nav <- rowMeans(is.na(m))
gid <- names(which(nav==0))
m <- m[gid,]

d <- readRDS('/hpc/group/jilab/zj/met/final/procrna.rds')
af <- list.files('/hpc/group/jilab/zj/met/scdata/HCL/saver')
names(af) <- sub('_','',sub('.rds','',af))
names(af)[names(af)=='AdultThyroid1'] <- 'AdultThyroidGland1'
names(af)[names(af)=='HESC1'] <- 'hESCs'
names(af)[names(af)=='AdultTransverseColon2'] <- 'AdultTransverseColon1'
names(af)[names(af)==''] <- 'PancreaticbetacellsderivedfromH9'
names(af)[names(af)==''] <- 'hiPSCs'
l <- read.table('/hpc/group/jilab/zj/met/scdata/HCL/list',as.is=T,sep='\t')
l[,1] <- gsub(' ','',l[,1])
af <- af[names(af) %in% l[,1]]

sc <- sapply(1:length(af),function(i) {
  tmp <- readRDS(paste0('/hpc/group/jilab/zj/met/scdata/HCL/saver/',af[i]))
  colnames(tmp) <- paste0(names(af)[i],'_',1:ncol(tmp))
  tmp
})
names(sc) <- names(af)
gene <- unique(unlist(sapply(sc,rownames)))
cell <- unique(unlist(sapply(sc,colnames)))
scm <- matrix(0,nrow=length(gene),ncol=length(cell),dimnames = list(gene,cell))
for (i in names(sc)) {
  print(i)
  scm[rownames(sc[[i]]),colnames(sc[[i]])] <- sc[[i]]
}
scm <- log2(scm + 1)

int <- intersect(rownames(scm),rownames(d))
d <- d[int,]
scm <- scm[int,]
#sampscm <- scm[,sample(1:ncol(scm),1000)]

testid <- l[,4]
trainid <- setdiff(colnames(m),testid)
hsdid <- sfunc(m[,testid])
mod <- trainmodel(d[,trainid],m[hsdid,trainid])
pred <- predict(scm,mod)
true <- m[hsdid,l[match(sub('_.*','',colnames(scm)),l[,1]),4]]
ct <- sub('_.*','',colnames(pred))

perf <- sapply(1:1000,function(sid) {
  sel <- as.vector(sapply(unique(ct),function(i) sample(colnames(pred)[which(ct==i)],50)))
  sel <- which(colnames(pred) %in% sel)
  sampcv <- corfunc(pred[,sel],true[,sel])
  citecv <- corfunc(t(pred[,sel]),t(true[,sel]))
  c(median(sampcv),median(citecv))
})

summary(sampcv)
summary(citecv)



bpred <- predict(d[,testid],mod)
btrue <- m[hsdid,testid]
sampcv <- corfunc(bpred,btrue)
citecv <- corfunc(t(bpred),t(btrue))
c(median(sampcv[1:10000]),median(citecv))

pd <- data.frame(cv=c(sampcv,citecv),type=rep(c('Across samples','Across CpG cites'),c(length(sampcv),length(citecv))),stringsAsFactors = F)

library(ggplot2)
pdf('/hpc/group/jilab/zj/met/res/plot/acrosssamplecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across samples',],aes(cv,type,col=type)) + geom_violin(col='royalblue') + geom_boxplot(col='royalblue',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation')
dev.off()

pdf('/hpc/group/jilab/zj/met/res/plot/acrosscitecor.pdf',width=2,height=3)
ggplot(pd[pd$type=='Across CpG cites',],aes(cv,type,col=type)) + geom_violin(col='orange') + geom_boxplot(col='orange',width=0.2) + coord_flip() + theme_classic() + theme(legend.position = 'none') + ylab('') + xlab('Correlation') + scale_color_manual(values=c('Across CpG cites'='orange','Across samples'='royalblue'))
dev.off()


id <- which.max(sampcv)
pdf('/hpc/group/jilab/zj/met/res/plot/acrosssampleexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[id,],y=true[id,]),aes(x=x,y=y)) + geom_point(col='royalblue') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none')
dev.off()

id <- which.max(citecv)
pdf('/hpc/group/jilab/zj/met/res/plot/acrossciteexp.pdf',width=3,height=3)
ggplot(data=data.frame(x=pred[,id],y=true[,id]),aes(x=x,y=y)) + geom_point(alpha=0.1,size=0.1,col='orange') + theme_classic() + xlab('Predicted methylation') + ylab('True methylation') + theme(legend.position = 'none')
dev.off()

