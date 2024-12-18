setwd('/home/whou10/data/whou/metpred/')
allcpg = readRDS('data/mousepancreas/wgbs/bs_cpgnames.rds')
str(allcpg)

## select the CpGs in gene promoters
suppressMessages(library(bsseq))
library(data.table)
gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))

pro <- promoters(gene,upstream=1000,downstream=0)
wgbs.cpg = GRanges(seqnames = sub(':.*', '', allcpg), 
                   IRanges(start = sub('.*:', '', allcpg)))
o <- as.matrix(findOverlaps(wgbs.cpg,pro))
cpg.id <- allcpg[o[,1]]
df = data.frame(cpg = allcpg[o[,1]], gene = names(pro)[o[,2]])
saveRDS(df, '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/promoterDNAm/res/df_of_cpg_gene.rds')        

