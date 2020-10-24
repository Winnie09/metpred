library(ggplot2)
library(RColorBrewer)
library(gridExtra)
pdir <- '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/plot/plot/sens_fdr_curve/'
dir.create(pdir, recursive = T)

for (seed in 1:100){
  print(seed)
  plist <- list()
  rdir1 <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/perf/', seed, '/sensfdr/')
  afd <- list.files(rdir1)
  for (fd in afd){
    print(fd)
    ddir <- paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simuseed/perf/', seed, '/sensfdr/', fd, '/')## 200_200_0.2_0.8
    af <- list.files(ddir)
    if (length(af) > 0){
      allsens <- lapply(af, function(f){
        m <- gsub('_.*','',sub('.rds.rds','',f))
        if (m == 'glmm') m <- 'GLMM'
        if (m == 'wilcoxon') m <- 'wilcox'
        if (m == 't') m <- 't-test'
        data.frame(readRDS(paste0(ddir, f)), Method = m)
      })
      pd3 <- do.call(rbind, allsens)
      pd3[,4] <- factor(pd3[,4])
      plist[[paste0(seed,':',fd)]] <- ggplot(data = pd3, aes(x = Real_FDR, y = Sensitivity, color=Method, group = Method)) +
        geom_line(size = 1) +
        theme_classic()  +
        geom_abline(slope = 1, intercept = 0, color = 'red')+
        # geom_vline(xintercept = 0.25, color = 'blue')+
        theme(axis.text = element_text(size = 8, color = 'black')) +
        xlim(c(0,0.1)) + 
        ggtitle(paste0(seed,':',fd)) + 
        scale_color_brewer(palette='Set1')  
    }
    
  }
  pdf(paste0(pdir, seed, '.pdf'), width = 17.5, height = 9)
  grid.arrange(grobs = plist, nrow = 4)  
  dev.off()
}

