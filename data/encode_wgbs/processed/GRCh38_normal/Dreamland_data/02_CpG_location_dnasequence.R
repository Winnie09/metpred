ddir <- '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/Dreamland_data/'
rdir <- '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/Dreamland_data/sub/'
i = as.numeric(commandArgs(trailingOnly = T)[1])
print(i)
# Load the required package
library(BSgenome.Hsapiens.UCSC.hg38)

# Load the mouse genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Extract the sequence for chromosome 000, positions aaa to bbb
# rdir <- ddir <- '/data/hji7/whou/metpred/tcga/hg38/combine/'

library(data.table)
# d = fread(paste0(ddir, 'CpG_name_location.csv'), sep = ',', data.table = F)
# print(dim(d))
rn = readRDS(paste0(ddir, 'wgbs_CpGnames.rds'))
str(rn)
rg = rn[1]

runid = seq(1e4*(i-1)+1, min(1e4*i, length(rn)))

allseq <- sapply(rn[runid], function(rg){
  print(rg)
  chr_sequence <- genome[[sub('_.*','',rg)]]
  region_sequence <- subseq(chr_sequence, start= as.numeric(sub('.*_','',rg))-49, end=as.numeric(sub('.*_','',rg))+50)
  # Print the sequence
  # print(region_sequence)
  as.character(region_sequence)
})
length(allseq)
# allseq = do.call(rbind, allseq)
dim(allseq)
df = cbind(CpG_location = rn[runid], DNA_sequence = allseq)
write.csv(df, paste0(rdir, 'CpG_name_location_DNAsequence_', i, '.csv'), row.names = F)


