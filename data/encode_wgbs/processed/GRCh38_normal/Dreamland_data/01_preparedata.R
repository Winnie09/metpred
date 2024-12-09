ddir <- '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/'
ddir2 = '/home/whou10/data/whou/metpred/data/tcga/hg38/combine/'
rdir1 = '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/Dreamland_data/'
rdir2 = '/home/whou10/data/whou/metpred/data/tcga/hg38/combine/Dreamland_data/'

rna = readRDS(paste0(ddir, 'rna.rds'))
ge = readRDS(paste0(ddir2, 'ge.rds'))
int = intersect(rownames(rna), rownames(ge))
rna2 = rna[int, , drop = F]
ge2 = ge[int, , drop = F]
saveRDS(rna2, paste0(rdir1, 'ge.rds'))
write.csv(rna2, paste0(rdir1, 'ge.csv'))
saveRDS(ge2, paste0(rdir2, 'ge.rds'))
write.csv(ge2, paste0(rdir2, 'ge.csv'))

w = readRDS(paste0(ddir, 'wgbs.rds'))
rn = rownames(w)
rn2 = sub(':', '_', rn)
rownames(w) = rn2

## remove chrX, chrY. Focus on actual CpG
a = rownames(w)
library(BSgenome.Hsapiens.UCSC.hg38)
seq <- sub('_.*','',a)
pos <- as.numeric(sub('.*_','',a))
k <- GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
ch <- Views(Hsapiens,k)
## if all are CpGs, then we only see "CG" in the table outputs
print(table(as.character(ch))) 
a <- a[seq %in% paste0('chr',1:22)] ## 27078450
d2 = w[a, ]
saveRDS(a, paste0(rdir1, 'me_cpgnames.rds'))
saveRDS(d2, paste0(rdir1, 'me_rownamesloc.rds'))
write.csv(d2, paste0(rdir1,'me_rownamesloc.csv'))

