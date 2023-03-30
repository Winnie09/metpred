suppressMessages(library(rtracklayer))
library(parallel)
af <- list.files('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/b')
d <- mclapply(af,function(f) {
  tmp <- import(paste0('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/b/',f))
  res <- unlist(mcols(tmp))
  names(res) <- paste0(as.character(seqnames(tmp)),'_',start(tmp))
  res
},mc.cores=detectCores())
names(d) <- af

s <- NULL
for (i in names(d)) s <- union(s,names(d[[i]]))
m <- matrix(NA,nrow=length(s),ncol=length(af),dimnames=list(s,af))
for (i in 1:length(d)) {
print(i)
  m[names(d[[i]]),i] <- d[[i]]
}
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/me.rds')
