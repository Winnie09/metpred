library(GenomicRanges)
suppressMessages(library(bsseq))
library(data.table)
for (technology in c('450k', 'EPIC')){
  # d: data.frame with columns CpG_name, CpG_location
  d = read.csv(paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/pos/', technology, '.csv'))
  # genecpglist: named list; names = gene, values = character vector of CpG_location strings
  genecpglist = readRDS(paste0('/home/whou10/scratch16/whou10/resource/genecpg/grch38_cpg_in_allgene_promotersTSS2k3Up_in_', technology, '.rds'))
  
  # 1) Invert genecpglist into long table: one row per (CpG_location, gene)
  inv <- data.table(
    gene         = rep(names(genecpglist), lengths(genecpglist)),
    CpG_location = unlist(genecpglist, use.names = FALSE)
  )
  
  # 2) Collapse to one row per CpG_location with all genes concatenated
  lookup <- inv[, .(Gene_promoter_TSS2k3Up = paste(unique(gene), collapse = ";")), by = CpG_location]
  
  # 3) Join back to d to annotate each row
  d_annot <- merge(as.data.table(d), lookup, by = "CpG_location", all.x = TRUE)
  
  # d_annot now has columns: CpG_location, CpG_name, gene_names
  # (rows with no match will have NA in gene_names)
  
  # Write to CSV
  write.csv(d_annot, file = paste0('/home/whou10/data/whou/metpred/data/tcga/hg38/2024/array/pos/', technology, '_with_promoter.csv'), row.names = FALSE, quote = FALSE)
  
}






