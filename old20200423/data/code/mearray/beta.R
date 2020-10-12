expr <- read.table('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/meta/experiment.tsv',as.is=T,sep='\t',quote='',header=T)
expr[,'Biosample.summary'] <- gsub(' ','_',expr[,'Biosample.summary'])
m <- read.table('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/meta/metadata.tsv',as.is=T,sep='\t',quote='',header=T)
m <- m[m[,'Platform']=='Illumina Infinium Methylation EPIC BeadChip',]
expr <- expr[expr[,'Accession'] %in% m[,'Experiment.accession'],]
load('/home-4/zji4@jhu.edu/scratch/encode_compiled/May19/metadata/RNAseq/processed/processed.rda')
hg19_experiment[,'Biosample summary'] <- gsub(' ','_',hg19_experiment[,'Biosample summary'])
sum <- intersect(expr[expr[,'Accession'] %in% m[,'Experiment.accession'],'Biosample.summary'],hg19_experiment[,'Biosample summary'])

beta <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/beta.rds')
betaexpr <- sub('_.*','',colnames(beta))

m <- sapply(sum,function(i) {
  tar <- intersect(betaexpr,expr[expr[,'Biosample.summary']==i,'Accession'])
  rowMeans(beta[,betaexpr %in% tar,drop=F])
})
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/mearray/beta.rds')
