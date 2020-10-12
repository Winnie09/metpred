load('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/metadata/RNAseq/processed/processed.rda')
hg19_experiment[,'Biosample summary'] <- gsub(' ','_',hg19_experiment[,'Biosample summary'])

e <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/hg19/RNAseq/matrix/FPKM.rds')
colnames(e) <- sub('.tsv','',colnames(e))
e <- log2(e + 1)
met <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
me <- sapply(colnames(met),function(ub) {
  af <- hg19_file[hg19_file[,'Experiment accession'] %in% hg19_experiment[hg19_experiment[,'Biosample summary']==ub,'Accession'],1]
  rowMeans(e[,af,drop=F])
})
suppressMessages(library(data.table))
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/hg19.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gid <- gsub('\"','',sub('gene_id ','',sapply(gtf[,9],function(i) strsplit(i,'; ')[[1]][1],USE.NAMES = F)))
me <- me[row.names(me) %in% gid,]
#me <- me[rowSums(me >= 1) >= 1,]
#cv <- apply(me,1,sd)/rowMeans(me)
#library(preprocessCore)
#dn <- dimnames(me)
#me <- normalize.quantiles(me)
#dimnames(me) <- dn
saveRDS(me,file='/home-4/zji4@jhu.edu/scratch/metpred/data/data/expr/mat.rds')


