ddir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/raw/download/wgbs/'
rdir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/'
af <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/matched_samples.rds')

me <- sapply(af, function(f){
  print(f)
  d = data.table::fread(paste0(ddir, f,'.bed.gz'), data.table=F)
  v = d[,4]
  names(v) <- paste0(d[,1], ':', d[,2]+1)
  v
})
saveRDS(me, paste0(rdir, 'me_cpg_by_sample_list.rds'))

a <- names(me[[1]])
for (i in 2:length(me)) {
  a <- intersect(a,names(me[[i]]))
}

## check if the CpGs are actual CpGs
library(BSgenome.Hsapiens.UCSC.hg38)
seq <- sub(':.*','',a)
pos <- as.numeric(sub('.*:','',a))
k <- GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
ch <- Views(Hsapiens,k)
## if all are CpGs, then we only see "CG" in the table outputs
print(table(as.character(ch))) 

a <- a[seq %in% paste0('chr',1:22)]

d <- sapply(me,function(i) i[a])
saveRDS(d, paste0(rdir, 'me_cpg_by_sample.rds'))


