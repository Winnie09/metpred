d <- readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/wgbs/combine/me_cpg_by_sample.rds')
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
seq <- sub('_.*','',rownames(d))
pos <- as.numeric(sub('.*_','',rownames(d)))
dna <- as.data.frame(getSeq(Hsapiens, seq, start = pos, end = pos + 1))[,1]
print(table(dna))
