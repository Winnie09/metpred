setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
res = readRDS('./metpred/expr_meth/plot/expr_meth_pd.rds')
g = res[['g']]
aggmeth = res[['aggmeth']]
expr = res[['expr']]
df = res[['df']]
gname = res[['gname']]

id = intersect(rownames(aggmeth), rownames(expr))
expr = expr[id,]
aggmeth = aggmeth[id,]
## plot expr vs. meth
library(ggplot2)
library(gridExtra)

expr = expr[rowMeans(expr>0)>0.5,]
aggmeth = aggmeth[rownames(expr),]
rownames(expr) <- rownames(aggmeth) <- gname[match(rownames(expr), g)]
expr = expr[unique(rownames(expr)),]
aggmeth = aggmeth[unique(rownames(expr)),]
for (i in 1:1000){
  png(paste0('./metpred/expr_meth/plot/expr_meth_',i,'.png'),width = 500,height = 500)
  plist = list()
  for (gene in rownames(expr)[16*(i-1) + 1:16]){
    print(gene)
    pd = data.frame(meth = aggmeth[which(rownames(aggmeth)==gene),], expr = expr[which(rownames(expr)==gene),], type = as.factor(df[,2]))
    str(pd)
    plist[[gene]] = ggplot() + geom_point(data=pd, aes(x= expr, y = meth, col=type)) +
      theme_classic() + ggtitle(gene) + theme(legend.position = 'none') + xlab('expressioin') + ylab('methylation')
  }
  grid.arrange(grobs=plist)
  dev.off()  
}


