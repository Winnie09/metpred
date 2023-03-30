d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/encode_wgbs/data/wgbs.rds')
load('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/RNAseq/meta/processed/processed.rda')
hg38_experiment[,'Biosample summary'] <- gsub(' ','_',hg38_experiment[,'Biosample summary'])

e <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/RNAseq/matrix/TPM.rds')
colnames(e) <- sub('.tsv','',colnames(e))
e <- log2(e + 1)

e <- sapply(colnames(d),function(s) {
  aexpe <- hg38_experiment[hg38_experiment[,'Biosample summary']==s,2]
  af <- hg38_file[hg38_file[,'Experiment accession'] %in% aexpe,1]
  rowMeans(e[,af,drop=F])
})
row.names(e) <- sub('\\..*','',row.names(e))
saveRDS(e,file='/home-4/zji4@jhu.edu/scratch/metpred/data/encode_wgbs/data/expr.rds')
