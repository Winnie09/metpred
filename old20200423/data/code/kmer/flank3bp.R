suppressMessages(library(GenomicRanges))
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifRG)
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
end(gr) <- start(gr) + 4
start(gr) <- start(gr) - 3
res <- getSequence(gr, BSgenome.Hsapiens.UCSC.hg19)
res <- as.matrix(res)
res <- apply(res,1,paste0,collapse='')
tab <- table(res)
tab <- sort(tab,decreasing = T)
#id <- which(cumsum(tab)/length(gr) < 0.5)
keep <- names(tab)[1:1000]
mat <- sapply(keep,function(i) res==i)
saveRDS(mat,file='/home-4/zji4@jhu.edu/scratch/metpred/data/kmer/flank3bp.rds')


