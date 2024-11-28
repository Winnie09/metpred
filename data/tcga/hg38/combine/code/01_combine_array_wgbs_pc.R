## check wgbs and array data 
## pc
## use intercept cpg
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38')
suppressMessages(library(GenomicRanges))

## read in wgbs
wgbs <- readRDS(paste0('wgbs/combine/me_cpg_by_sample.rds'))  
wgbscpg = rownames(wgbs)


## read in array data
array = readRDS('450k_2021/combine/Methylformer_data/me_rownamesloc.rds')
arraycpg = rownames(array)

## cbind wgbs and array, intersect cpg
int = intersect(wgbscpg, arraycpg) ##  chr [1:320132] "chr1_631826" "chr1_631932" 
wgbsint = wgbs[int, ]
colnames(wgbsint) = paste0('wgbs:', colnames(wgbsint))
arrayint = array[int, ]
colnames(arrayint) = paste0('450k:', colnames(arrayint))
d = cbind(wgbsint, arrayint)

## intersect genes
wgbs_rna <- readRDS(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/ge.rds'))
array_rna <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/ge.rds')
rownames(array_rna) <- sub(';.*', '', rownames(array_rna))
intg = intersect(rownames(wgbs_rna), rownames(array_rna))
wgbs_rna2 <- wgbs_rna[intg,]
array_rna2 <- array_rna[intg,]
colnames(wgbs_rna2) <- paste0('wgbs:', colnames(wgbs_rna2))
colnames(array_rna2) <- paste0('wgbs:', colnames(array_rna2))
expr = cbind(wgbs_rna2, cbind(array_rna2))


## promoters
gtf <- data.table::fread('/home/whou10/data/whou/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gr) <- gn

pro <- promoters(gr,upstream=1000,downstream = 1000)
seq <- sub('_.*','',rownames(d)) ## 
pos <- as.numeric(sub('.*_','',rownames(d))) ## 

dgr <- GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
o <- as.data.frame(findOverlaps(dgr,pro))
cb <- rowsum(d[o[,1],],names(pro)[o[,2]])
tab <- table(names(pro)[o[,2]])
cb <- cb/as.vector(tab[rownames(cb)])

## intersect promoter dnam and expr
int <- intersect(rownames(cb),rownames(expr))
cb <- cb[int,]
g <- g[int,]

saveRDS(list(dnam = d,
             expr = expr,
             promoterdnam_int = cb, 
             expr_int = g), '/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/array_wgbs_combine.rds')

## calculate statistics
rdir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/combine/res/'
sd1 <- apply(cb,1,sd)
sd2 <- apply(g,1,sd)
int <- which(sd1 > 0.1 & sd2 > 0.1)
cv <- sapply(1:ncol(cb),function(i) cor(cb[int,i],g[int,i]))
saveRDS(cv, paste0(rdir, 'cv.rds'))

r <- sapply(1:ncol(cb),function(i1) {
  v <- sapply(1:ncol(cb),function(i2) {
    cor(cb[int,i1],g[int,i2])
  }) 
  (v[i1]-mean(v[-i1]))/sd(v[-i1])
})
saveRDS(r, paste0(rdir, 'r.rds'))


## principal components
pr <- prcomp(t(g[int,]),scale=T)$x
saveRDS(pr, paste0(rdir, 'dnam_pr.rds'))

pr <- prcomp(t(cb[int,]),scale=T)$x
saveRDS(pr, paste0(rdir, 'rna_pr.rds'))

