af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/rnasep')
d <- sapply(af,function(f) {
  readRDS(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/rnasep/',f))
})
gn <- unique(unlist(sapply(d,names)))
m <- matrix(0,nrow=length(gn),ncol=length(af),dimnames = list(gn,af))
for (i in names(d)) {
  m[names(d[[i]]),i] <- d[[i]]
}
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/ge.rds')
