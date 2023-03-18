library(ggplot2)
library(gridExtra)
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/res/res/'
for (disease in c('cancer','normal')) {
  d <- readRDS(paste0(rdir, 'predicted_DNAm_on_testset_', disease,'.rds'))
  d <- d[apply(d,1,sd) > 0.1,]
  pr <- prcomp(t(d),scale=T)$x
  stu <- sub('_.*','',rownames(pr))
  tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/doc/ct_tissue_info_for_matched_wgbs_and_RNAseq.csv', header = TRUE, as.is = TRUE)
  tis <- tb[match(rownames(pr),tb[,1]),2]
  p1 <- ggplot(data.frame(PC1=pr[,1],PC2=pr[,2],tissue=tis),aes(x=PC1,y=PC2,col=tissue)) + geom_point() + theme_classic()
  p2 <- ggplot(data.frame(PC1=pr[,1],PC2=pr[,2],tissue=tis),aes(x=PC1,y=PC2,col=tissue)) + geom_point() + theme_classic() + facet_wrap(~tissue)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/res/plot/predictedvalue_pca/',disease,'.pdf'),width=14,height=7)
  grid.arrange(p1,p2,nrow=1)
  dev.off()
}

