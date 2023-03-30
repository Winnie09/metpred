suppressMessages(library(rtracklayer))
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/wgbs')
n <- sapply(af,function(f) {
  print(f)
  d <- import(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/wgbs/',f))
  m <- unlist(mcols(d))
  names(m) <- paste0(as.character(d),'_',start(d))
  m
})

r <- NULL
for (i in names(n)) {
  r <- union(r,names(n[[i]]))
}
m <- matrix(NA,nrow=length(r),ncol=length(n),dimnames=list(r,af))
for (i in names(n)) {
  print(i)
  m[names(n[[i]]),i] <- n[[i]]
}
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/me.rds')
