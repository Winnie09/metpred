suppressMessages(library(GenomicRanges))
d1 <- readRDS(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/me_cpg_by_sample.rds'))  
saveRDS(rownames(d1), '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/cpg_names.rds')

d = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/me.rds')
saveRDS(rownames(d), '/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/cpg_names.rds')
             
int = intersect(rownames(d1), rownames(d))
saveRDS(int, 'intersect.rds')


g <- readRDS(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/ge.rds'))
gtf <- data.table::fread('/home/whou10/data/whou/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gr) <- gn

pro <- promoters(gr,upstream=1000,downstream = 1000)
seq <- sub(':.*','',rownames(d))
pos <- as.numeric(sub('.*:','',rownames(d)))

dgr <- GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
o <- as.data.frame(findOverlaps(dgr,pro))
cb <- rowsum(d[o[,1],],names(pro)[o[,2]])
tab <- table(names(pro)[o[,2]])
cb <- cb/as.vector(tab[rownames(cb)])

int <- intersect(rownames(cb),rownames(g))
cb <- cb[int,]
g <- g[int,]

sd1 <- apply(cb,1,sd)
sd2 <- apply(g,1,sd)
int <- which(sd1 > 0.1 & sd2 > 0.1)
cv <- sapply(1:ncol(cb),function(i) cor(cb[int,i],g[int,i]))
saveRDS(cv, 'cv.rds')

r <- sapply(1:ncol(cb),function(i1) {
  v <- sapply(1:ncol(cb),function(i2) {
    cor(cb[int,i1],g[int,i2])
  }) 
  (v[i1]-mean(v[-i1]))/sd(v[-i1])
})
saveRDS(r, 'r.rds')
pr <- prcomp(t(g[int,]),scale=T)$x
saveRDS(pr, 'pr.rds')
