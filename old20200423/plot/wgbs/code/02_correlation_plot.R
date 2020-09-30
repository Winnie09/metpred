rm(list=ls())
setwd('/scratch/users/whou10@jhu.edu/Wenpin')
# setwd('/Users/wenpinhou/Dropbox/')
res = readRDS('./metpred/plot/wgbs/plot/correlation.pd.rds')
m = res[['m']]
e = res[['e']]
o = res[['o']]
lcv = res[['lcv']]

pdf('./metpred/plot/wgbs/plot/hist_acrossample.pdf')
hist(lcv)
dev.off()
scv <- sapply(1:ncol(m),function(i) cor(m[o[,1],i],e[o[,2],i]))
pdf('./metpred/plot/wgbs/plot/hist_acrossite.pdf')
hist(scv)
dev.off()

pd  = data.frame(y=m[o[,1],'GM12878'], x=e[o[,2],'GM12878'])
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
library(viridis)
pd$density = get_density(pd$x, pd$y, n = 20)
library(ggplot2)
p <- ggplot() + geom_point(data=pd,aes(x=x,y=y,col=density),alpha=0.2,size=0.2)+
  theme_classic() + xlab('Expression')+ylab('Methylation')+
  theme(legend.position = 'right',legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6,alpha=1)))+
  theme(legend.spacing.x = unit(-0.1, 'cm'),legend.spacing.y = unit(-0.3, 'cm'))+
  scale_color_viridis()

pdf('./metpred/plot/wgbs/plot/GM12878_acrosssite.pdf',width=5,height=5)
print(p)
dev.off()

pdf('./metpred/plot/wgbs/plot/GM12878_acrosssite2.pdf',width=3.5,height=3.5)
smoothScatter(pd$y~pd$x)
dev.off()

pdf('./metpred/plot/wgbs/plot/acrossample.pdf')
plot(m[o[,1],][100000,]~e[o[,2],][100000,],xlab='Expression',ylab='Methylation')
dev.off()

