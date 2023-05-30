me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')

#type <- commandArgs(trailingOnly = T)
type = 'ramp'
i <- 10
d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/pred/',type,'/',i,'.rds'))

c <- readRDS('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/ramp/cpg.rds')

## a sample across CpG
library(ggplot2)
pdir <- '/home/whou10/data/whou/metpred/evaluate/encode_cv/plot/'
for (s in 1:ncol(d)){
  data <- data.frame(x = me[,s], y = d[, s], stringsAsFactors = F)
  set.seed(1)
  data2 <- data[sample(1:nrow(data), 10^5), ]
  
  
  
  pdf(paste0(pdir, 'Example_sample_across_CpGs_', s, '.pdf'), width = 2.2, height = 1.6)
  print(ggplot(data2, aes(x=x, y=y) ) +
          geom_bin2d(bins = 100) +
          scale_fill_continuous(type = "viridis") +
          geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.2)+
          #theme_bw() +
          xlab('Measured DNAm')+
          ylab('Predicted DNAm'))
  dev.off()  
}


## a CpG acorss sample
c.sub <- c[1 > c[,1] & c[,1] > 0.98 & c[,4]==1, , drop = FALSE]
str(c.sub)
summary(c.sub[,1])

set.seed(1)
for (cpg in sample(rownames(c.sub), 10)){
  print(cpg)
  data <- data.frame(x = me[cpg, colnames(d)], y = d[cpg, ], stringsAsFactors = F)
  pdf(paste0(pdir, paste0('Example_CpG_across_samples_', cpg, '.pdf')), width = 1.6, height = 1.6)
  print(ggplot(data, aes(x=x, y=y) ) +
          geom_point(size = 1) +
          #theme_bw() +
          geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.2)+
          xlab('Measured DNAm')+
          ylab('Predicted DNAm')+
          ggtitle(paste0('PCC=',round(c.sub[cpg, 1], 2))))
  dev.off()
}

