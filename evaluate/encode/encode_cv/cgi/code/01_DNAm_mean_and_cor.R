
## ==================================
ddir = '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/pred/ramp/'
rdir = '/home/whou10/data/whou/metpred/evaluate/encode/encode_cv/cgi/res/'
## read in CpG island (CGI) location information
tb = read.table('/home/whou10/data/whou/metpred/data/cpg_island/hgTables.txt', stringsAsFactors = FALSE)
library(GenomicRanges)
gr <- GRanges(seqnames=tb[,2],IRanges(start=tb[,3],end=tb[,4]))
# gr <- gr[as.character(seqnames(gr))  %in% paste0('chr', c(1:22, 'X'))]

library(data.table)
gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
# gene <- gene[as.character(seqnames(gene)) %in% paste0('chr', c(1:22, 'X'))]
pro <- promoters(gene,upstream=2e3,downstream=0)
o <- as.matrix(findOverlaps(pro, gr)) ## caution: order
tb.tmp = tb[o[,2],]
o.df <- data.frame(gene = names(pro)[o[,1]],
                   cgi = paste0(tb.tmp[,2], ':', tb.tmp[,3], '_', tb.tmp[,4]))

## cgi that overlap with promoter
cgi = gr[unique(o[,2])] 

## cpg that are within these cgi
rn = readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs_CpGnames.rds')
start <- as.numeric(sub('.*:', '', rn))
cpg <- GRanges(seqnames= sub(':.*', '', rn), IRanges(start=start,end=start))
o2 = as.matrix(findOverlaps(cgi, cpg))
cpg_in_cgi = cpg[unique(o2[,2])]

o2.df = data.frame(cgi = paste0(seqnames(cgi), ':', start(cgi), '_', end(cgi))[o2[,1]],
                   cpg = paste0(seqnames(cpg), ':', start(cpg), '_', end(cpg))[o2[,2]])

saveRDS(o.df, paste0(rdir,'o.df.rds'))
saveRDS(o2.df, paste0(rdir,'o2.df.rds'))

## load encode wgbs DNAm
## mean DNAm of these cgi
me = readRDS(paste0(ddir, '1.rds'))
res = t(sapply(unique(o2[,1]), function(i){
  tmp = o2[o2[,1]==i,2,drop=FALSE]
  colMeans(me[tmp, ,drop=F])
}))
rownames(res) = paste0(seqnames(cgi[unique(o2[,1])]), ':' , start(cgi[unique(o2[,1])]), '_', end(cgi[unique(o2[,1])]))
saveRDS(res, paste0(rdir, 'mean_DNAm_of_cgi_that_overlap_with_promoter2e3Up.rds'))

## sd of these cgi
source('/home/whou10/scratch16/whou10/resource/myfunc/01_function.R')
res2 = sapply(unique(o2[,1]), function(i){
  print(i)
  tmp = o2[o2[,1]==i,2]
  rowsds(me[tmp, ,drop=F])
})
names(res2) = paste0(as.character(seqnames(cgi[unique(o2[,1])])), ':' , start(cgi[unique(o2[,1])]), '_', end(cgi[unique(o2[,1])]))
saveRDS(res2, paste0(rdir, 'DNAm_sd_of_cgi_that_overlap_with_promoter2e3Up.rds'))

#### calculate the PCC  among the DNAm of the CpG within a same CGI
cgi.cor = sapply(unique(o2[,1]), function(i){
  print(i)
  tmp = o2[o2[,1]==i,2]
  c = corfunc(t(tmp), t(tmp))
  lower_tri_matrix <- c
  lower_tri_matrix[!lower.tri(lower_tri_matrix, diag=TRUE)] <- 0
  mean(as.vector(lower_tri_matrix))
})
names(cgi.cor) = paste0(seqnames(cgi[unique(o2[,1])]), ':' , start(cgi[unique(o2[,1])]), '_', end(cgi[unique(o2[,1])]))
saveRDS(cgi.cor, paste0(rdir, 'DNAm_cpg_cor_of_cgi_that_overlap_with_promoter2e3Up.rds'))


