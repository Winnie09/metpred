mat = readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38mat/mat.rds')
m <- read.table('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/meta/metadata.tsv',as.is=T,sep='\t',header=T)

em <- read.table('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/meta/experiment.tsv',as.is=T,sep='\t',header=T)

v = gsub('_',' ',colnames(mat))
ems <- em[em[,'Biosample.summary']%in% v,]
# ems <- ems[,c('Accession','Biosample.summary','Organism','Replicates','Biosample.term.name')]
ems <- ems[order(ems[,'Biosample.summary']),]
write.csv(ems,'/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/meta/dnam/encode_wgbs/match_experiment.csv')

