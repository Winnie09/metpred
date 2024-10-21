rm(list=ls())
rdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/ramp/'
d = readRDS(paste0(rdir, 'predicted_cpg_sd0.3.rds'))
gs = readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/eval/eval_goldstandard_ctmean.rds')
int = intersect(rownames(d), rownames(gs))
print(str(int))
gs = gs[int, ]
d = d[int, ]
str(gs)
colnames(gs) = paste0('WGBS_', colnames(gs))
str(d)
pd = cbind(d, gs)

source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
pr = pca_lmFilter(genebycellmat = pd, topnum = 2e3)
u = UMAP(pr$x)

# pr = PCA(pd, findVariableGenes = FALSE)
# u = UMAP(pr)
saveRDS(pr, paste0(rdir, 'pred_dnam_and_bulk_pr.rds'))
saveRDS(u, paste0(rdir, 'pred_dnam_and_bulk_umap.rds'))

## ====================
## plot
## =====================
pdir = '/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/plot/'
u = readRDS(paste0(rdir, 'pred_dnam_and_bulk_umap.rds'))
type = sapply(rownames(u), function(i) ifelse(grepl('WGBS',i), 'WGBS', 'Predicted'))
rownames(u) = sub('WGBS_','',rownames(u))
pd = data.frame(u1 = u[,1], u2 = u[,2], 
                celltype = sub(':.*', '', rownames(u)),
                type = type)
table(pd[,3])
table(pd[,4])
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')
pdf(paste0(pdir, 'pred_dnam_and_bulk_umap.pdf'), width = 1.8, height = 1.8)
# plotumap(plotdata = pd)
# dev.off()

table(pd[,3])
plotdata1 = pd[pd[,4]=='Predicted',]
plotdata2 = pd[pd[,4]!='Predicted',]
plotdata = pd
tpd = data.frame(u1 = tapply(plotdata[,1], plotdata[,3], mean),
                 u2 = tapply(plotdata[,2], plotdata[,3], mean),
                 celltype = unique(plotdata[,3]))
ggplot() + 
  geom_point(data = plotdata1, aes(x = u1, y = u2, color = celltype),
             size = 0.1, stroke = 0, alpha = 0.7)+
  geom_point(data = plotdata2, aes(x = u1, y = u2),
             size = 1, stroke = 0, alpha = 0.7)+
  xlab('UMAP1') + 
  ylab('UMAP2') +
  geom_text(data = tpd, aes(x = u1, y = u2, color = celltype, label = celltype),
            size = 2.2, alpha = 0.8) +
  theme(legend.position = 'none') +
  # scale_color_manual(values = colorRampPalette(rev(
  #   brewer.pal(8, 'Set1')
  # ))(length(unique(
  #   plotdata[, 3]
  # ))))  +
  xlab('UMAP1')

dev.off()
u[2100:2103,]
