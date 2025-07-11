i = as.numeric(commandArgs(trailingOnly = T)[1])
print(i)
# Load the required package
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
# Load the mouse genome
genome <- BSgenome.Hsapiens.UCSC.hg38
# Extract the sequence for chromosome 000, positions aaa to bbb
rdir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/sub/'
ddir <- '/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/combine/'
rn = readRDS(paste0(ddir, 'cpg_names.rds'))
str(rn)

runid = seq(1e4*(i-1)+1, min(1e4*i, length(rn)))
allseq <- sapply(rn[runid], function(rg){
  print(rg)
  chr_sequence <- genome[[sub('_.*','',rg)]]
  region_sequence <- subseq(chr_sequence, start= as.numeric(sub('.*_','',rg))-49, end=as.numeric(sub('.*_','',rg))+50)
  # Print the sequence
  # print(region_sequence)
  as.character(region_sequence)
})
print(length(allseq))
df = cbind(d[runid, ,drop=FALSE], DNA_sequence = allseq)
write.csv(df, paste0(rdir, 'CpG_name_location_DNAsequence_', i, '.csv'), row.names = F)

