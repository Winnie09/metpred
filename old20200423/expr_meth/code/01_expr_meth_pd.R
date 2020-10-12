setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/')
library(GenomicRanges)
library(matrixStats)
set.seed(12345)
mat = readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38filtermat/mat.rds')
gr = readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38filtermat/gr.rds')

tss <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/geneTSS.rds')

n = colnames(mat)
n = n[order(nchar(n))]
df = data.frame(celltype=colnames(mat), type='normal',stringsAsFactors = F)
df[match(c('A549', 'K562', 'HepG2', 'SK-N-SH','OCI-LY7','HeLa-S3'), df[,1]),2] <- 'cancer'


## gene expr
expr = readRDS('./metpred/data/data/expr/mat.rds')
expr = expr[rowSums(expr>1)>1, ]
cv = apply(expr,1,sd)/rowMeans(expr)
summary(cv)
expr = expr[cv >1, ]
rownames(expr) = sub('\\..*','',rownames(expr))
id = intersect(colnames(mat), colnames(expr))
expr = expr[, id]
mat = mat[,id]
## gene tss 200bp region
g = sub('.*;','',names(tss))
gname = sub(';.*','',names(tss))
g = intersect(g, rownames(expr))
gname = gname[match(g, sub('.*;','',names(tss)))]
up <- promoters(tss[match(g,rownames(expr))], upstream = 2000,downstream = 1000)


## gene tss 200bp region methylation sum
o = as.matrix(findOverlaps(gr, up))
v = rep(0,length(gr))
v[o[,1]] <- o[,2]
library(Matrix.utils)
agg = aggregate.Matrix(mat, v, fun='sum')
agg = as.matrix(agg)
rownames(agg) <- sapply(rownames(agg),function(i) {
  tmp  = as.numeric(i)
  ifelse(tmp==0,'none', g[tmp])
})
names(rownames(agg)) = NULL
id = intersect(rownames(agg), rownames(expr))
saveRDS(list(g = g, gname=gname, aggmeth=agg, expr=expr, df=df),'./metpred/expr_meth/plot/expr_meth_pd.rds')

