## check tcga array data 
## promoter DNAm and gene expression pcc
## compared to unpaired array-rnaseq
setwd('/home/whou10/data/whou/metpred/data/tcga/hg38')
suppressMessages(library(GenomicRanges))

## read in wgbs
wgbs <- readRDS(paste0('wgbs/combine/me_cpg_by_sample.rds'))  
wgbscpg = rownames(wgbs)


## read in array data
array = readRDS('450k_2021/combine/Methylformer_data/me_rownamesloc.rds')
arraycpg = rownames(array)

## intersect cpg
int = intersect(wgbscpg, arraycpg) ##  chr [1:320132] "chr1_631826" "chr1_631932" 
d = array[int, ]


## load in array's rna-seq samples
## find promoters
## calculate promoter dnam
g <- readRDS(paste0('450k_2021/combine/Methylformer_data/ge.rds'))

## shuffle sample names, so that rnaseq and array have unpaired samples
set.seed(1)
sample <- sample(colnames(g))
colnames(g) <- sample

gtf <- data.table::fread('/home/whou10/data/whou/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gr) <- gn

pro <- promoters(gr,upstream=1000,downstream = 1000)
seq <- sub(':.*','',rownames(d)) ## 
pos <- as.numeric(sub('.*:','',rownames(d))) ## 

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
saveRDS(cv, '450k_2021/combine/Methylformer_data/unpaired_sample_cv.rds')

r <- sapply(1:ncol(cb),function(i1) {
  v <- sapply(1:ncol(cb),function(i2) {
    cor(cb[int,i1],g[int,i2])
  }) 
  (v[i1]-mean(v[-i1]))/sd(v[-i1])
})
saveRDS(r, '450k_2021/combine/Methylformer_data/unpaired_sample_r.rds')

