library(RCurl)
uf <- function(u) {
  strsplit(getURL(u,dirlistonly=T,ftp.use.epsv=F),'\n')[[1]]
}

res <- NULL
f1 <- uf('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/')
for (f1s in f1) {
  print(f1s)
  f2 <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',f1s,'/'))
  for (f2s in f2) {
    f3 <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',f1s,'/',f2s,'/'))
    for (f3s in f3) {
      f4 <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',f1s,'/',f2s,'/',f3s,'/'))
      res <- rbind(res,data.frame(f1s,f2s,f3s,f4))
    }
  }
}

type <- apply(res[,-4],1,paste0,collapse=':')
t1 <- type[res[,4]=='Bisulfite-Seq']
t2 <- type[res[,4]=='RNA-Seq']
int <- intersect(t1,t2)

burl <- sapply(int,function(i) {
  print(i)
  i <- which(res[,4]=='Bisulfite-Seq' & type == i)
  samp <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/'))
  print(length(samp))
  af <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/',samp,'/'))
  paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/',samp,'/',af[grepl('CPG_methylation_calls.bs_call',af)])
}) 


rurl <- sapply(int,function(i) {
  print(i)
  i <- which(res[,4]=='RNA-Seq' & type == i)
  samp <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/'))
  print(length(samp))
  af <- uf(paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/',samp,'/'))
  paste0('ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/',res[i,1],'/',res[i,2],'/',res[i,3],'/',res[i,4],'/',samp,'/',af[grepl('gene_quantification',af) & grepl('results',af)])
}) 

saveRDS(burl,file='/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/burl.rds')
saveRDS(rurl,file='/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/rurl.rds')
