suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
seq <- toupper(as.data.frame(getSeq(Hsapiens, as.character(seqnames(gr)), start = start(gr)-100, end = end(gr) + 101))[,1])
saveRDS(seq,file='/home-4/zji4@jhu.edu/scratch/metpred/data/data/seq/seq.rds')
chrmat <- do.call(rbind,strsplit(seq,''))
chrmat <- matrix(as.numeric(cbind(chrmat=='A',chrmat=='T',chrmat=='G',chrmat=='C')),nrow=nrow(chrmat))
saveRDS(chrmat,file='/home-4/zji4@jhu.edu/scratch/metpred/data/data/seq/onehotmat.rds')

suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
seq <- toupper(as.data.frame(getSeq(Hsapiens, as.character(seqnames(gr)), start = start(gr)-100, end = end(gr) + 101))[,1])
saveRDS(seq,file='/home-4/zji4@jhu.edu/scratch/metpred/data/data/seq/seq.rds')
chrmat <- do.call(rbind,strsplit(seq,''))
chrmat <- matrix(as.numeric(cbind(chrmat=='A',chrmat=='T',chrmat=='G',chrmat=='C')),nrow=nrow(chrmat))
saveRDS(chrmat,file='/home-4/zji4@jhu.edu/scratch/metpred/data/data/seq/onehotmat.rds')


