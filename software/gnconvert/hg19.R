library(data.table)
d <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/hg19.gtf',data.table=F)
d <- d[d[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',d[,9]))
gid <- sub('\\..*','',sub('".*','',sub('gene_id "','',d[,9])))
m <- cbind(gn,gid)
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/software/gnconvert/hg19.rds')
