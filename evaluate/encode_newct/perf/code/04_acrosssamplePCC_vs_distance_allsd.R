pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
apd = readRDS(paste0(pdir, 'acrosssample.rds'))

ddir1 = '/home/whou10/data/whou/metpred/evaluate/encode_newct/pc/rna_cosineSimilarity/'
ddir2 = '/home/whou10/data/whou/metpred/evaluate/encode_newct/pc/wgbs_cosineSimilarity/'
dis <- sapply(1:10, function(cvid){
  dis_rna = 1-readRDS(paste0(ddir1, cvid, '.rds'))
  dis_wgbs = 1-readRDS(paste0(ddir2, cvid, '.rds'))
  c(median(apply(dis_rna, 1, median, na.rm=T), na.rm=T),
  median(apply(dis_wgbs, 1, median, na.rm=T), na.rm=T))
})
dimnames(dis) = list(c('dis_rna', 'dis_wgbs'), 1:10)

pd = rbind(cbind(apd, distance = dis[1,match(apd[,1],colnames(dis))], type = 'RNA_dis'),
          cbind(apd, distance = dis[2,match(apd[,1],colnames(dis))], type = 'WGBS_dis'))


library(ggplot2)
pdir = '/home/whou10/data/whou/metpred/evaluate/encode_newct/perf/plot/compare/'
pdf(paste0(pdir, '/acrosssamplePCC_vs_distance_allsd.pdf'), width=4.5,height=3)
ggplot(pd[pd[,2]=='Ramp',],aes(x=distance,y=value,fill = type, group = distance)) + 
  geom_boxplot(alpha=0.3, outlier.colour = NA, lwd=0.3) + 
  geom_point(size = 0.5, stroke = 0, alpha = 0.8) + 
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Distance between training and testing') + 
  ylab('Across-sample PCC') + 
  #scale_color_manual(values=pal) + 
  theme(axis.text.x = element_text(angle = 30,vjust=0.5),legend.title = element_blank()) 
# scale_y_continuous(breaks = c(-0.25,0, 0.25, 0.5,0.75))
dev.off()

