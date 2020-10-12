# nSample = as.numeric(commandArgs(trailingOnly = T)[1]) ## number of sample in each group, 2, 4
# nCell = as.numeric(commandArgs(trailingOnly = T)[2]) ## number of cells in each sample: 2, 20, 100, 200, 400
# propCpG = as.numeric(commandArgs(trailingOnly = T)[3]) ## 0.2, 0.3, 0.4
# probRead= as.numeric(commandArgs(trailingOnly = T)[4]) ## 0.6, 0.8, 1
# nSample <- 4
# nCell <- 200
# propCpG <- 0.2
# probRead <- 0.6
for (nSample in c(2, 4)){
  for (nCell in c(2, 20, 100, 200, 400)){
    for (propCpG in c(0.2, 0.3, 0.4)){
      for (probRead in c(0.6, 0.8, 1)){
        fullmat <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/data/data/', nSample,'To', nSample, '_allcell/pred_ery.rds'))
        ap <- as.numeric(sub('BM','',sub(':.*', '', colnames(fullmat))))
        as <- sub('_.*','', colnames(fullmat))
        
        mat <- lapply(unique(as), function(i){
          set.seed(12345)
          id <- sample(which(as == i), nCell, replace = TRUE)
          fullmat[, id]
        })
        mat <- do.call(cbind, mat)
        colnames(mat) <- paste0(colnames(mat),1:ncol(mat))
        
        logit <- function(x) log(x/(1-x))
        logistic <- function(x) exp(x)/(1+exp(x))
        set.seed(12345)
        addgene <- sample(rownames(mat),propCpG*nrow(mat))
        cellid <- which(as.numeric(sub('BM','',sub(':.*','',colnames(mat)))) %in% c(1,2,5,6))
        addmat <- mat[sample(setdiff(rownames(mat),rownames(addgene)),length(addgene)),sample(colnames(mat),length(cellid))]
        mat[addgene, cellid] <- logistic(logit(mat[addgene, cellid]) + logit(addmat * probRead))
        
        t1 <- rowMeans(mat[addgene, cellid])
        t2 <- rowMeans(mat[addgene, setdiff(1:ncol(mat),cellid)])
        summary(abs(t1-t2))
        t1 <- rowMeans(mat[setdiff(rownames(mat),addgene), cellid])
        t2 <- rowMeans(mat[setdiff(rownames(mat),addgene), setdiff(1:ncol(mat),cellid)])
        summary(abs(t1-t2))
        rownames(mat)[rownames(mat) %in% addgene] <- paste0('add:',rownames(mat)[rownames(mat) %in% addgene])
        
        ## save files
        setwd('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/data/')
        dir <- paste0(nSample,'_',nCell,'_',propCpG,"_",probRead, '/')
        dir.create(dir, showWarnings = FALSE)
        saveRDS(mat,file=paste0(dir,'me.rds'))
      }
    }
  }
}

        
