## read in CpG island (CGI) location information
tb = read.table('/home/whou10/data/whou/metpred/data/cpg_island/hgTables.txt', stringsAsFactors = FALSE)
# v = tb[,4] - tb[,3]
# > summary(v)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201.0   325.0   568.0   776.9   959.0 45712.0 
upshore1 <- upshore2 <- upshore3 <- upshore4 <- tb[, 1:4]
colnames(upshore1) <- colnames(upshore2) <- colnames(upshore3) <- colnames(upshore4) <- c('#bin',  'chrom', 'chromStart',  'chromEnd')

upshore1[,3] <- tb[,3]-500
upshore1[,4] <- tb[,3]-1

upshore2[,3] <- tb[,3]-1000
upshore2[,4] <- tb[,3]-501

upshore3[,3] <- tb[,3]-1500
upshore3[,4] <- tb[,3]-1001

upshore4[,3] <- tb[,3]-2000
upshore4[,4] <- tb[,3]-1501

write.table(upshore1, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore1.txt', quote = F)
write.table(upshore2, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore2.txt', quote = F)
write.table(upshore3, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore3.txt', quote = F)
write.table(upshore4, '/home/whou10/data/whou/metpred/data/cpg_island/hgTables_upshore4.txt', quote = F)
