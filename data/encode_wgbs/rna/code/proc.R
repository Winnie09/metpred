library(data.table)
meta <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/rna.rds')
conv <- c('GRCh38'='Homo sapiens','mm10'='Mus musculus')
aaf <- sub('.tsv','',list.files('/home/whou10/data/whou/metpred/data/encode_wgbs/rna/tsv'))

for (org in c('GRCh38','mm10')) {
  af <- intersect(aaf,meta[which(meta$Organism==conv[org]),'file'])
  d <- sapply(af,function(f) {
    d <- fread(paste0('/home/whou10/data/whou/metpred/data/encode_wgbs/rna/tsv/',f,'.tsv'),data.table=F)
    v <- d$TPM
    names(v) <- d$gene_id
    v
  },simplify = F)
  print(mean(sapply(2:length(d),function(i) identical(names(d[[1]]),names(d[[i]])))))
  d <- do.call(cbind,d)
  d <- log2(d+1)
  d <- d[grep('^ENS',rownames(d)),]
  rownames(d) <- sub('\\..*','',rownames(d))
  saveRDS(d,file=paste0('/home/whou10/data/whou/metpred/data/encode_wgbs/rna/proc/',org,'.rds'))
}


