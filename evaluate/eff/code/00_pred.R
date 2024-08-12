ddir <- '/home/whou10/data/whou/metpred/evaluate/eff/data/processed/'
afd <- list.files(ddir)
fd <- afd[as.numeric(commandArgs(trailingOnly = T)[1])]
print(fd)
rdir <- paste0('/home/whou10/data/whou/metpred/evaluate/eff/pred/', fd, '/')
dir.create(rdir, recursive = T)
print(rdir)

if (file.exists(paste0(rdir, '1.rds'))) {
  stop('The prediction has been completed, so the script has to be ended.')
}

ge <- readRDS(paste0(ddir, fd, '/rna.rds'))
me <- readRDS(paste0(ddir, fd, '/wgbs.rds'))  

# new program
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
set.seed(12345)
cid <- as.numeric(cut(1:ncol(me), 10))
names(cid) <- sample(colnames(me))

# i <- as.numeric(commandArgs(trailingOnly = T)) 
i <- 1
testid <- which(cid == i)
trainid <- which(cid != i)
pred <-
  trainpredict(trainexpr = ge[, trainid],
               testexpr = ge[, testid],
               trainmeth = me[, trainid])
saveRDS(pred, file = paste0(rdir, i, '.rds'))

