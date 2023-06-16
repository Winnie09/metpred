## very long to load
w <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')
r <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
source('/home/whou10/scratch16/whou10/resource/startup.R')
setwd('/home/whou10/data/whou/metpred/')
rdir <- 'evaluate/mousepancreas/spatialpred/validity/'

## select the CpGs in genes' promoters
suppressMessages(library(bsseq))
library(data.table)
gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('\\..*','',sub('".*','',sub('gene_id "','',gtf[,9])))
gene <- gene[intersect(names(gene), rownames(r))]
pro <- promoters(gene,upstream=1000,downstream=0)
wgbs.cpg = GRanges(seqnames = sub(':.*', '', rownames(w)), 
                   IRanges(start = sub('.*:', '', rownames(w))))
o <- as.matrix(findOverlaps(wgbs.cpg,pro))
cpg.id <- rownames(w)[o[,1]]
agg <- aggregate(w[cpg.id, ], list(names(pro)[o[,2]]), 'mean')
rownames(agg) <- agg[,1]
agg <- as.matrix(agg[,-1])
saveRDS(agg, paste0(rdir, 'res/encode_human_promoter_dnam.rds'))

int = intersect(rownames(r), rownames(agg))
agg = agg[int, ]
r = r[int,]
cc = corfunc(agg, r)
cc2 = corfunc(t(agg), t(r))
pdf(paste0(rdir, 'plot/encode_PCC.pdf'), width = 6, height = 3.3)
par(mfrow=c(1,2))
hist(cc, main='gene across-sample', xlab = 'PCC', breaks = 20)
hist(cc2, main='sample across-gene', xlab = 'PCC', breaks = 20)
dev.off()

