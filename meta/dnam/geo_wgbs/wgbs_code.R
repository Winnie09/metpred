af <- list.files('Downloads/wgbs/')
ar <- NULL
for (f in af) {
  res <- readRDS(paste0('Downloads/wgbs/',f))
  res <- res[res[,'Type']=='Bisulfite-Seq',]
  ar <- rbind(ar,res)
}
arv <- apply(ar,1,paste0,collapse=';')
ar <- ar[grep('wgbs',arv,ignore.case = T),]
ar <- ar[ar[,'Organism'] %in% c('Homo sapiens','Mus musculus'),]
table(ar[,'Organism'])
write.csv(ar,file='Downloads/wgbstable.csv')
