# ## =================================================
# ## for each cv, plot histrogram of acrosssample pcc
# ## =================================================
# ## read in a cv perf
# key = 'neargene'
# rdir <- paste0('/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/res/', key, '/')
# sample = readRDS(paste0(rdir, 'sample.rds'))
# cpg = readRDS(paste0(rdir, 'cpg.rds'))
# 
# pdir = paste0('/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/',key,'/')
# pdf(paste0(pdir, 'across_sample_cor.pdf'), width = 20, height = 20)
# par(mfrow=c(5,4))
# for (i in 1:10){
#   tmp = cpg[cpg[,4]==i,,drop=F]
#   hist(tmp[,1], main = paste0(i, ': across-sample PCC'))
#   hist(tmp[,2], main = paste0(i, ': across-sample RMSE'))
# }
# dev.off()

## ==================================================
## plot across-sample pcc for all method for each cv
## ==================================================
library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL

cvid = 1
mpdlist <- list()
for (cvid in 1:10){
  type = 'acrosssample'
  rdir <- '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/res/'
  p1 <- data.frame(method='Nearest gene',type=type,readRDS(file=paste0(rdir, 'neargene/', type, '/cv', cvid, '.rds')),stringsAsFactors = F)
  p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir, 'permu/', type, '/cv', cvid, '.rds')),stringsAsFactors = F)
  p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir, 'ramp/', type, '/cv', cvid, '.rds')),stringsAsFactors = F)
  pd <- rbind(p1,p2,p3)
  
  # type = 'all'
  # rdir <- '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/res/'
  # p1 <- data.frame(method='Nearest gene',type=type,readRDS(file=paste0(rdir, 'neargene/all/cpg.rds')),stringsAsFactors = F)
  # p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir, 'permu/all/cpg.rds')),stringsAsFactors = F)
  # p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir, 'ramp/all/cpg.rds')),stringsAsFactors = F)
  # pd <- rbind(p1,p2,p3)
  # > str(pd)
  # 'data.frame':	84905217 obs. of  6 variables:
  #   $ method: chr  "Nearest gene" "Nearest gene" "Nearest gene" "Nearest gene" ...
  # $ type  : chr  "acrosssample" "acrosssample" "acrosssample" "acrosssample" ...
  # $ cor   : num  -0.704 -0.705 -0.717 -0.721 -0.724 ...
  # $ rmse  : num  0.224 0.223 0.219 0.218 0.216 ...
  # $ sd    : num  0.207 0.207 0.204 0.203 0.203 ...
  # $ cvid  : int  1 1 1 1 1 1 1 1 1 1 ...
  # 
  pd$sdcut <- cut(pd$sd,seq(0,1,0.05))
  
  pd$method <- factor(as.character(pd$method),levels=c('Ramp','Nearest gene','Permute'))
  
  # mpd <- tapply(pd$cor,list(pd$method,pd$sdcut,pd$cvid),median,na.rm=T)
  mpdlist[[cvid]] <- tapply(pd$cor,list(pd$method,pd$sdcut),median,na.rm=T)
}  
library(reshape2)
apd <- do.call(rbind,sapply(1:length(mpdlist),function(i) {
  data.frame(cvid=i,melt(mpdlist[[i]]))
},simplify = F))

# ## ===========================================
# ## option 1: if multiple cv all in mpd
# apd <- do.call(rbind,sapply(1:10,function(i) {
#   data.frame(cvid=i,melt(mpd[,,i,drop=F]))
# },simplify = F))
# ## option 2: else if only one cv in mpd
# # apd = data.frame(cvid=1,melt(mpd[,,1,drop=F]))
# ## ===========================================

apd <- apd[!is.na(apd$value),]
# apd <- apd[apd$Var2 %in% c('(0,0.05]','(0.05,0.1]','(0.1,0.15]','(0.15,0.2]','(0.2,0.25]','(0.25,0.3]'),]
apd <- apd[!apd$Var2 %in% 'NA',]

# pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
# pdf(paste0(pdir, '/acrosssample_cv', cvid, '.pdf'), width=3.7,height=3)
# print(ggplot(apd,aes(x=Var2,y=value,color=Var1)) + geom_boxplot() + theme_classic() + theme(legend.position = 'bottom') + xlab('Standard deviation across samples') + ylab('PCC across samples') + scale_color_manual(values=pal) + theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) + scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75)))
# dev.off()
pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
saveRDS(apd, paste0(pdir, 'acrosssample.rds'))

pdf(paste0(pdir, '/acrosssample.pdf'), width=3.7,height=3)
ggplot(apd,aes(x=Var2,y=value,fill = Var1)) + 
  geom_boxplot(alpha = 0.2, outlier.colour = NA) + 
  geom_point(position=position_dodge(width=0.75),aes(group=Var1), 
             alpha = 0.5, stroke = 0, size = 0.7) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Across-sample CpG variability') + 
  ylab('Across-sample PCC') + scale_color_manual(values=pal) + 
  scale_fill_manual(values=pal)+
  theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) 
  # scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75))
dev.off()

## ================================
################ only use top 5 cv
## ================================
apd <- do.call(rbind,sapply(1:5,function(i) {
  data.frame(cvid=i,melt(mpdlist[[i]]))
},simplify = F))


apd <- apd[!is.na(apd$value),]
# apd <- apd[apd$Var2 %in% c('(0,0.05]','(0.05,0.1]','(0.1,0.15]','(0.15,0.2]','(0.2,0.25]','(0.25,0.3]'),]
apd <- apd[!apd$Var2 %in% 'NA',]

saveRDS(apd, paste0(pdir, 'acrosssample_5cv.rds'))

# pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
# pdf(paste0(pdir, '/acrosssample_cv', cvid, '.pdf'), width=3.7,height=3)
# print(ggplot(apd,aes(x=Var2,y=value,color=Var1)) + geom_boxplot() + theme_classic() + theme(legend.position = 'bottom') + xlab('Standard deviation across samples') + ylab('PCC across samples') + scale_color_manual(values=pal) + theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) + scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75)))
# dev.off()
pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
pdf(paste0(pdir, '/acrosssample_5cv.pdf'), width=3.7,height=3)
ggplot(apd,aes(x=Var2,y=value,fill = Var1)) + 
  geom_boxplot(alpha = 0.2, outlier.colour = NA) + 
  geom_point(position=position_dodge(width=0.75),aes(group=Var1), 
             alpha = 0.5, stroke = 0, size = 0.7) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Across-sample CpG variability') + 
  ylab('Across-sample PCC') + 
  scale_color_manual(values=pal) + 
  scale_fill_manual(values=pal) + 
  theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) 
# scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75))
dev.off()



########################## =====================
library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL


rdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/res/'
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0(rdir,'neargene/all/sample.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir,'permu/all/sample.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir,'ramp/all/sample.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)

apd <- tapply(pd$cor,list(pd$method,pd$cvid),median,na.rm=T)
library(reshape2)
apd <- melt(apd)
colnames(apd) <- c('method','cvid','cor')
apd$method <- factor(as.character(apd$method),levels=c('Ramp','Nearest gene','Permute'))

pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
pdf(paste0(pdir, 'cpg.pdf'),width=3,height=2.2)
ggplot(apd,aes(x=method,y=cor,fill=method)) + 
  geom_boxplot(alpha = 0.2, width = 0.5) + 
  geom_point(stroke = 0, size = 0.6, alpha = 0.6)+
  scale_fill_manual(values=pal)+
  theme_classic() + theme(legend.position = 'none') + 
  xlab('Method') + 
  ylab('Across-CpG PCC')
dev.off()
