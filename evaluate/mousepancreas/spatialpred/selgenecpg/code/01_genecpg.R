setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
dir.g = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/consistency/res/'
dir.r = '/home/whou10/data/whou/metpred/evaluate/mousepancreas/spatialpred/selgenecpg/res/'

expr <- readRDS(file="/home/whou10/data/whou/metpred/data/mousepancreas/spatial/procimpute/res/all/imputednormexpr.rds")
alls <- sub(':.*', '', colnames(expr))

for (s in unique(alls)){
  print(s)
  rs = rowsds(expr[, alls == s])
  prome <-  readRDS('/home/whou10/data/whou/FeinbergLab/mousePancreas/wgbs/promoter_methylation/promoter_tssup1kb_averaged_dnam.rds')
  int = intersect(rownames(prome), names(rs))
  rs = rs[int]
  prome = prome[int, ]
  gene.sel <- readRDS(paste0(dir.g, 'leaveout_', s, '_selgene1e3.rds'))
  
  ## select test sample 
  test.s <- s
  train.s <- setdiff(unique(alls), test.s)
  
  ## load cell type specific expression and wgbs data
  r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/imputed_nnls/res/all/cse.rds')
  w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')
  
  ## select the CpGs in 1k genes' promoters
  suppressMessages(library(bsseq))
  library(data.table)
  gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
  gtf <- gtf[gtf[,3]=='gene',]
  gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
  names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
  gene <- gene[gene.sel]
  pro <- promoters(gene,upstream=1000,downstream=0)
  wgbs.cpg = GRanges(seqnames = sub(':.*', '', rownames(w)), 
                     IRanges(start = sub('.*:', '', rownames(w))))
  o <- as.matrix(findOverlaps(wgbs.cpg,pro))
  cpg.id <- rownames(w)[o[,1]]
  saveRDS(data.frame(cpg = rownames(w)[o[,1]], gene = names(pro)[o[,2]]), 
          paste0(dir.r, s, '.rds'))

}

  
