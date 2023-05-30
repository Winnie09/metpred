library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'sample'
print(type)
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/neargene/',type,'.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/permu/',type,'.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/ramp/',type,'.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)

apd <- tapply(pd$cor,list(pd$method,pd$cvid),median,na.rm=T)
library(reshape2)
apd <- melt(apd)
colnames(apd) <- c('method','cvid','cor')
apd$method <- factor(as.character(apd$method),levels=c('Ramp','Nearest gene','Permute'))

pdf('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/cpg.pdf',width=3,height=3)
ggplot(apd,aes(x=method,y=cor,fill=method)) + geom_boxplot() + theme_classic() + theme(legend.position = 'none') + xlab('') + ylab('PCC across CpGs')
dev.off()


library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'difcor'
print(type)
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/neargene/',type,'.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/permu/',type,'.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/ramp/',type,'.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)

apd <- tapply(pd$cor,list(pd$method,pd$cvid),median,na.rm=T)
library(reshape2)
apd <- melt(apd)
colnames(apd) <- c('method','cvid','cor')
apd$method <- factor(as.character(apd$method),levels=c('Ramp','Nearest gene','Permute'))

pdf('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/difcor.pdf',width=3,height=3)
ggplot(apd,aes(x=method,y=cor,fill=method)) + geom_boxplot() + theme_classic() + theme(legend.position = 'none') + xlab('') + ylab('PCC across CpGs')
dev.off()



library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'cpg'
print(type)
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/neargene/',type,'.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/permu/',type,'.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/ramp/',type,'.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)
pd$sdcut <- cut(pd$sd,seq(0,1,0.05))

pd$method <- factor(as.character(pd$method),levels=c('Ramp','Nearest gene','Permute'))

mpd <- tapply(pd$cor,list(pd$method,pd$sdcut,pd$cvid),median,na.rm=T)

library(reshape2)
apd <- do.call(rbind,sapply(1:10,function(i) {
  data.frame(cvid=i,melt(mpd[,,i]))
},simplify = F))
apd <- apd[!is.na(apd$value),]
apd <- apd[apd$Var2 %in% c('(0,0.05]','(0.05,0.1]','(0.1,0.15]','(0.15,0.2]','(0.2,0.25]','(0.25,0.3]'),]

pdf('/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/sample.pdf',width=3.7,height=3)
ggplot(apd,aes(x=Var2,y=value,color=Var1)) + geom_boxplot() + theme_classic() + theme(legend.position = 'bottom') + xlab('Standard deviation across samples') + ylab('PCC across samples') + scale_color_manual(values=pal) + theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) + scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75))
dev.off()



