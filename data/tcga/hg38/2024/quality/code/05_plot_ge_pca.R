## check tcga gene expression batch effect
## across those corresponding to three DNA methylation tech
rdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/res/'
pdir = '/home/whou10/data/whou/metpred/data/tcga/hg38/2024/quality/plot/'

source('/home/whou10/data/whou/resource/myfunc/01_function.R')
prfull = readRDS(paste0(rdir, 'prfull.rds'))
pr = t(prfull$x[,1:10])
sdev = ((prfull$sdev)^2/sum((prfull$sdev)^2))[1:10]

pd = data.frame(pc1 = pr[1,],
                pc2 = pr[2,],
                project = sub('-.*','', colnames(pr)))
rownames(pd) = colnames(pr) ## error
pd$project[grepl('^C',pd$project)] <- 'C'
pd$project[grepl('^11LU',pd$project)] <- '11LU'

wgbssp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/wgbs_samples.txt')
epicsp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/epic_samples.txt')
arraysp = readLines('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/450k_samples.txt')

pd$technology = rep('WGBS', nrow(pd))
pd$technology[rownames(pd) %in% epicsp] = 'EPIC'
pd$technology[rownames(pd) %in% arraysp] = '450K'
pd$technology[rownames(pd) %in% intersect(arraysp, wgbssp)] = 'WGBS & 450K'

pro = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/combine/project.rds')
str(pro)
cancer = pro$Cancer
names(cancer) = pro$Sample
pd$cancer = cancer[rownames(pd)]

source('/home/whou10/scratch16/whou10/resource/myfunc/ggplot_theme.R')
library(ggplot2)
pdf(paste0(pdir, 'expr_pc_technology.pdf'), width = 4, height = 2.2)
ggplot(pd) + geom_point(aes(x = pc1, y = pc2, col = technology), stroke = 0, size = 0.8, alpha = 0.5) + xlab('PC1') + ylab('PC2')
dev.off()

pdf(paste0(pdir, 'expr_pc_project.pdf'), width = 4, height = 2.2)
ggplot(pd) + geom_point(aes(x = pc1, y = pc2, col = project), stroke = 0, size = 0.8, alpha = 0.5) + xlab('PC1') + ylab('PC2')
dev.off()

pdf(paste0(pdir, 'expr_pc_cancer.pdf'), width = 5, height = 2.2)
ggplot(pd) + geom_point(aes(x = pc1, y = pc2, col = cancer), stroke = 0, size = 0.8, alpha = 0.5) + xlab('PC1') + ylab('PC2')
dev.off()


