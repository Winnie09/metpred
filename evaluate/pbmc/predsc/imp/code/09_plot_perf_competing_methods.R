rm(list=ls())
library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'cpg'
print(type)
rdir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/'
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0(rdir, 'neargene/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir,'permu/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir, 'ramp/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)

pd$method <- factor(as.character(pd$method),levels=c('Ramp','Nearest gene','Permute'))

pdir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
pdf(paste0(pdir, 'across_cpg_pcc.pdf'),width=2.2,height=2.6)
# ggplot(pd, aes(x = method, y = cor, fill = method)) + 
#   geom_boxplot() + 
#   theme_classic() + 
#   theme(legend.position = 'none') + 
#   scale_fill_manual(values = pal) +
#   xlab('Method') + 
#   ylab('Across-CpG PCC')
ggplot(pd,aes(x=method,y=cor,fill=method)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 0.3) + 
  geom_jitter(alpha = 0.5, stroke = 0, size = 0.5, width = 0.3) +
  theme_classic() + 
  xlab('Method') + ylab('Across-CpG PCC') +
  scale_fill_manual(values = pal) +
  ggtitle('Compare DNAm') +
  theme(legend.position = 'none',
        title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7)) 
dev.off()



# library(ggplot2)
# pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
# pd <- NULL
# type <- 'difcor'
# print(type)
# p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/evalpred/perf/neargene/',type,'.rds')),stringsAsFactors = F)
# p2 <- data.frame(method='Permute',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/evalpred/perf/permu/',type,'.rds')),stringsAsFactors = F)
# p3 <- data.frame(method='Ramp',type=type,readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/evalpred/perf/ramp/',type,'.rds')),stringsAsFactors = F)
# pd <- rbind(p1,p2,p3)
# 
# pd$method <- factor(as.character(pd$method),levels=c('Ramp','Nearest gene','Permute'))
# 
# pdf('/home/whou10/data/whou/metpred/evaluate/encode_gtex/evalpred/plot/difcor.pdf',width=2.2,height=2.6)
# # ggplot(pd,aes(x=method,y=cor,fill=method)) + 
# #   geom_boxplot() + theme_classic() + theme(legend.position = 'none') + 
# #   scale_fill_manual(values = pal) +
# #   xlab('Method') + ylab('Across-CpG difference PCC')
# ggplot(pd,aes(x=method,y=cor,fill=method)) + 
#   geom_boxplot(alpha = 0.5, outlier.size = 0.3) + 
#   geom_jitter(alpha = 0.5, stroke = 0, size = 0.5, width = 0.3) +
#   theme_classic() + 
#   xlab('Method') + ylab('Across-CpG difference PCC') +
#   scale_fill_manual(values = pal) +
#   ggtitle('Compare DNAm difference') +
#   theme(legend.position = 'none',
#         title = element_text(size = 7),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 7)) 
# dev.off()
# 


library(ggplot2)
pal <- c('Ramp'='red','Nearest gene'='orange','Permute'='royalblue')
pd <- NULL
type <- 'sample'
print(type)
rdir <- '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/'
p1 <- data.frame(method='Nearest gene',type=type,readRDS(paste0(rdir, 'neargene/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p2 <- data.frame(method='Permute',type=type,readRDS(paste0(rdir,'permu/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
p3 <- data.frame(method='Ramp',type=type,readRDS(paste0(rdir, 'ramp/pred_var_10000_CpG_across',type,'_pcc.rds')),stringsAsFactors = F)
pd <- rbind(p1,p2,p3)
pd$sd <- as.factor(as.character(pd$sd))
str(pd)

pd$method <- factor(as.character(pd$method),levels=c('Ramp','Nearest gene','Permute'))

apd <- tapply(pd$cor,list(pd$method,pd$sd),median,na.rm=T)
library(reshape2)
apd <- melt(apd)
apd <- apd[!is.na(apd$value),]
saveRDS(apd, paste0(pdir, 'across_sample_plotdata.rds'))

pdf(paste0(pdir, 'across_sample_pcc.pdf'),width=3.5,height=2.8)
ggplot(apd,aes(x=Var2,y=value,color=Var1, group=Var1)) + 
  geom_line() +
  geom_point() + 
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Across-sample CpG variability') + 
  ylab('Across-sample PCC') + 
  scale_color_manual(values=pal) + 
  theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) + 
  scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75,1))
dev.off()

