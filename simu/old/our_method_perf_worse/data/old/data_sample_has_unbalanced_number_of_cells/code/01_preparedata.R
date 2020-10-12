numcell1 = as.numeric(commandArgs(trailingOnly = T)[1])
numcell2 = as.numeric(commandArgs(trailingOnly = T)[2])
propCpG = as.numeric(commandArgs(trailingOnly = T)[3])
probRead= as.numeric(commandArgs(trailingOnly = T)[4])

fullmat <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/data/data/4To4_allcell/pred_ery.rds')
ap <- as.numeric(sub('BM','',sub(':.*', '', colnames(fullmat))))
as <- sub('_.*','', colnames(fullmat))
# BM4:29:male
# BM7:36:female
# 
# BM3:39:male
# BM8:32:female

if (numcell1 == numcell2 & numcell1 == 4){
  set.seed(12345)
  id1 <- c(sample(which(ap %in% c(4)), 2), sample(which(ap %in% c(7)), 2))
  id2 <- c(sample(which(ap %in% c(3)), 2), sample(which(ap %in% c(8)), 2))
} else if (numcell1 == 4){
  set.seed(12345)
  id1 <- c(sample(which(ap %in% c(4)), 2), sample(which(ap %in% c(7)), 2))
  id2 <- sample(which(!ap %in% c(3,8)), numcell2)
} else if (numcell2 == 4){
  set.seed(12345)
  id1 <- sample(which(ap %in% c(3,8)), numcell1)
  id2 <- c(sample(which(ap %in% c(4)), 2), sample(which(ap %in% c(7)), 2))
} else {
  set.seed(12345)
  id1 <- sample(which(ap %in% c(3,4,7,8)), numcell1)
  id2 <- sample(which(!ap %in% c(3,4,7,8)), numcell2)
  while (sum(table(ap[id1]) > 2) != length(unique(ap[id1]))){
    id1 <- sample(which(ap %in% c(3,4,7,8)), numcell1)
    id2 <- sample(which(!ap %in% c(3,4,7,8)), numcell2)
  }
}
id <- as.vector(c(id1,id2))

  
mat <- fullmat[,colnames(fullmat)[id]]
colnames(mat)[1:numcell1] <- paste0(colnames(mat)[1:numcell1],':1')
colnames(mat)[(numcell1+1):(numcell1+numcell2)] <- paste0(colnames(mat)[(numcell1+1):(numcell1+numcell2)],':2')

add <- fullmat[,setdiff(colnames(fullmat),colnames(mat))]
gn <- which(rowSums(mat) > 0)
mat <- mat[gn,]
oriadd <- add[gn,]

sampid <- sample(1:ncol(oriadd))
sampc1id <- sampid[1:length(id1)]
sampc2id <- sampid[(1+length(id1)):(length(id1)+length(id2))]

add1 <- sapply(sampc1id,function(i) {
  if (probRead< 1) {
    ssv <- oriadd[,i]  ## the counts form a cell i
    tarlib <- round(sum(ssv) * probRead) # library
    read <- rep(names(ssv),ssv*1e3) ## each site has a weight as its signal
    read <- sample(read,tarlib) ## sample the sites as many times as the library

    v <- rep(0,length(ssv))  # 0 vector for all sites
    names(v) <- names(ssv)
    tab <- table(read) # the signals for each site
    v[names(tab)] <- as.vector(tab) # give the signal back to the sites
    v
  } else if (probRead== 1) {
    oriadd[,i]
  } else {
    rowSums(oriadd[,sample(colnames(oriadd),probRead)])
  }
})
# add1 <- add1[rowSums(add1) > 0,] ## for probRead = 0

add2 <- sapply(sampc2id,function(i) {
  if (probRead< 1) {
    ssv <- oriadd[,i]
    tarlib <- round(sum(ssv) * probRead)
    read <- rep(names(ssv),ssv * 1e3)
    read <- sample(read,tarlib)
    v <- rep(0,length(ssv))
    names(v) <- names(ssv)
    tab <- table(read)
    v[names(tab)] <- as.vector(tab)
    v
  } else if (probRead== 1) {
    oriadd[,i]
  } else {
    rowSums(oriadd[,sample(colnames(oriadd),probRead)])
  }
})
## add2 <- add2[rowSums(add2) > 0,] ## for probRead = 0
## remove low wxpr genes
sampgn <- sample(intersect(row.names(add1),row.names(add2)),nrow(mat)*propCpG)     
simumat <- mat
g <- sampgn[1:round(length(sampgn)/2)]
dg1 <- g
simumat[g,1:numcell1] <- simumat[g,1:numcell1] + add1[g,]      
g <- sampgn[(round(length(sampgn)/2)+1):length(sampgn)]
dg2 <- g
simumat[g,((numcell1+1):(numcell1+numcell2))] <- simumat[g,((numcell1+1):(numcell1+numcell2))] + add2[g,]
simumat <- simumat[rowMeans(simumat > 0) >= 0.01,]
## save files
setwd('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/')
system(paste0('mkdir ',numcell1,'_',numcell2,'_',propCpG,"_",probRead))
dir.create(paste0(numcell1,'_',numcell2,'_',propCpG,"_",probRead), showWarnings = FALSE)
saveRDS(dg1,file=paste0(numcell1,'_',numcell2,'_',propCpG,"_",probRead,'/diffgn1.rds'))
saveRDS(dg2,file=paste0(numcell1,'_',numcell2,'_',propCpG,"_",probRead,'/diffgn2.rds'))
saveRDS(simumat,file=paste0(numcell1,'_',numcell2,'_',propCpG,"_",probRead,'/me.rds'))
saveRDS(sampgn,file=paste0(numcell1,'_',numcell2,'_',propCpG,"_",probRead,'/diffgn.rds'))



