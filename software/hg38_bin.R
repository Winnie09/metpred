suppressMessages(library(GenomicRanges))
bin <- readRDS('/home/whou10/data/whou/metpred/data/grbin/hg38.rds')
binfunc <- function(dat) {
  nbin <- paste0(as.character(seqnames(bin)),':',start(bin),'-',end(bin))
  chr <- sub(':.*','',rownames(dat))
  pos <- as.numeric(sub('.*:','',rownames(dat)))
  gr <- GRanges(seqnames=chr,IRanges(start=pos,end=pos))
  o <- as.matrix(findOverlaps(gr,bin))
  d <- rowsum(dat[o[,1],],nbin[o[,2]])
  d <- d/as.vector(table(nbin[o[,2]])[rownames(d)])
}

