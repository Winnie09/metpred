## read in CpG island (CGI) location information
tb = read.table('/home/whou10/data/whou/metpred/data/cpg_island/hgTables.txt', stringsAsFactors = FALSE)
# v = tb[,4] - tb[,3]
# > summary(v)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201.0   325.0   568.0   776.9   959.0 45712.0 
shelve <- upshore1 <- upshore2 <- upshore3 <- upshore4 <- tb[, 1:4]
colnames(upshore1) <- colnames(upshore2) <- colnames(upshore3) <- colnames(upshore4) <- c('#bin',  'chrom', 'chromStart',  'chromEnd')

upshore1[,3] <- tb[,3]-500
upshore1[,4] <- tb[,3]-1

upshore2[,3] <- tb[,3]-1000
upshore2[,4] <- tb[,3]-501

upshore3[,3] <- tb[,3]-1500
upshore3[,4] <- tb[,3]-1001

upshore4[,3] <- tb[,3]-2000
upshore4[,4] <- tb[,3]-1501

shelve[,3] <- tb[,3]-4000
shelve[,4] <- tb[,3]-2001

write.table(upshore1, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore1.txt', quote = F)
write.table(upshore2, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore2.txt', quote = F)
write.table(upshore3, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore3.txt', quote = F)
write.table(upshore4, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore4.txt', quote = F)
write.table(shelve, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_shelve.txt', quote = F)

library(GenomicRanges)
cgi <- GRanges(seqnames=tb[,2],IRanges(start=tb[,3],end=tb[,4]))
upshore1 <- GRanges(seqnames=upshore1[,2],IRanges(start=upshore1[,3],end=upshore1[,4]))
upshore2 <- GRanges(seqnames=upshore2[,2],IRanges(start=upshore2[,3],end=upshore2[,4]))
upshore3 <- GRanges(seqnames=upshore3[,2],IRanges(start=upshore3[,3],end=upshore3[,4]))
upshore4 <- GRanges(seqnames=upshore4[,2],IRanges(start=upshore4[,3],end=upshore4[,4]))
shelve <- GRanges(seqnames=shelve[,2],IRanges(start=shelve[,3],end=shelve[,4]))

## =========================================================
## identify CpG's location. Prepare a table of two columns:
## =========================================================
## cpgname, cgi/upshore1/upshore2/upshore3/shore4/shelve/sea
cpglocation <- function(cpg){
  o = as.matrix(findOverlaps(cgi, cpg))
  d = data.frame(cpg = cpgname[o[,2]], location = 'cgi', cgiIndex = o[,1])
  
  o = as.matrix(findOverlaps(upshore1, cpg))
  d = rbind(d,
            data.frame(cpg = cpgname[o[,2]], location = 'upshore1', cgiIndex = o[,1]))
  
  
  o = as.matrix(findOverlaps(upshore2, cpg))
  d = rbind(d,
            data.frame(cpg = cpgname[o[,2]], location = 'upshore2', cgiIndex = o[,1]))
  
  o = as.matrix(findOverlaps(upshore3, cpg))
  d = rbind(d,
            data.frame(cpg = cpgname[o[,2]], location = 'upshore3', cgiIndex = o[,1]))
  
  o = as.matrix(findOverlaps(upshore4, cpg))
  d = rbind(d,
            data.frame(cpg = cpgname[o[,2]], location = 'upshore4', cgiIndex = o[,1]))
  
  o = as.matrix(findOverlaps(shelve, cpg))
  d = rbind(d,
            data.frame(cpg = cpgname[o[,2]], location = 'shelve', cgiIndex = o[,1]))
  
  cpgsea = setdiff(cpgname, unique(d[,1]))
  d = rbind(d,
            data.frame(cpg = cpgsea, location = 'sea', cgiIndex = NA))
  return(d)
}

## TCGA WGBS
cpgname = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/cpg_names.rds')
seq = sub('_.*','',cpgname)
pos = as.numeric(sub('.*_','',cpgname))
cpg = GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
res = cpglocation(cpg)
write.table(res, '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/cpg_names_cgi_location.txt', quote = F)

## TCGA array
cpgname = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/me_cpgnames.rds')
seq = sub('_.*','',cpgname)
pos = as.numeric(sub('.*_','',cpgname))
cpg = GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
res = cpglocation(cpg)
write.table(res, '/home/whou10/data/whou/metpred/data/tcga/hg38/450k_2021/combine/Methylformer_data/cpg_names_cpg_location.txt', quote = F)

## ENCODE WGBS
cpgname = readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/Dreamland_data/me_cpgnames.rds')
seq = sub('_.*','',cpgname)
pos = as.numeric(sub('.*_','',cpgname))
cpg = GRanges(seqnames=seq,IRanges(start=pos,end=pos+1))
res = cpglocation(cpg)
write.table(res, '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/Dreamland_data/cpg_names_cpg_location.txt', quote = F)


