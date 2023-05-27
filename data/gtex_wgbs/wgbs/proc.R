suppressMessages(library(GenomicRanges))
library(rtracklayer)
af <- list.files('/hpc/group/jilab/zj/metpred/data/gtex_wgbs/wgbs/bw')
d <- sapply(af,function(f) {
  d <- import(paste0('/hpc/group/jilab/zj/metpred/data/gtex_wgbs/wgbs/bw/',f))
  v <- paste0(as.character(seqnames(d)),':',start(d))
  k <- unlist(mcols(d))
  names(k) <- v
  k
},simplify = F)

d <- do.call(cbind,d)
colnames(d) <- sub('.mCG.small_smooth.bw','',colnames(d))
saveRDS(d,file='/hpc/group/jilab/zj/metpred/data/gtex_wgbs/wgbs/proc.rds')
