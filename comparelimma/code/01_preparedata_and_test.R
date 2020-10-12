nSample = as.numeric(commandArgs(trailingOnly = T)[1]) ## number of sample in each group, 2, 4
nCell = as.numeric(commandArgs(trailingOnly = T)[2]) ## number of cells in each sample: 2, 20, 100, 200, 400
propCpG = as.numeric(commandArgs(trailingOnly = T)[3]) ## 0.2, 0.3, 0.4
probRead= as.numeric(commandArgs(trailingOnly = T)[4]) ## 0.6, 0.8, 1
nSample = 2
nCell = 2
propCpG = 0.2
probRead= 0.3
logistic <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))
source('/home-4/whou10@jhu.edu/work-zfs/whou10/trajectory_variability/function/01_function.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/software/raisin/raisin.R')
source('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/software/findDMR.R')
fun <- func[['limma_saver']]


fullmat <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/data/data/', nSample,'To', nSample, '_allcell/pred_ery.rds'))
ap <- as.numeric(sub('BM','',sub(':.*', '', colnames(fullmat))))
as <- sub('_.*','', colnames(fullmat))
ind <- NULL

Res <- NULL
for (seed in 1:3){
  mat <- lapply(unique(as), function(i){
    # set.seed(9884)
    set.seed(seed)
    #set.seed(6944)
    id <- sample(which(as == i), nCell, replace = TRUE)
    fullmat[, id]
  })
  mat <- do.call(cbind, mat)
  colnames(mat) <- paste0(colnames(mat), '_', 1:ncol(mat))
  
  ## add signals
  set.seed(12345)
  addgene <- sample(rownames(mat),propCpG*nrow(mat))
  cellid <- which(as.numeric(sub('BM','',sub(':.*','',colnames(mat)))) %in% c(1,2,5,6))
  addmat <- mat[sample(setdiff(rownames(mat),rownames(addgene)),length(addgene)),sample(colnames(mat),length(cellid))]
  mat[addgene, cellid] <- logistic(logit(mat[addgene, cellid]) + logit(addmat * probRead))
  rownames(mat)[rownames(mat) %in% addgene] <- paste0('add:',rownames(mat)[rownames(mat) %in% addgene])
  
  ## prepare for test
  sample <- sub(':.*','',colnames(mat))
  samplename <- unique(sample)
  sum(as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 0 == 1) < length(samplename)
  group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
  
  expr <- logit(mat) ## transform to -inf, +inf
  
  ## test our method
  design <- data.frame(sample=samplename,feature=group)
  
  # res1 <- RAISINtest(RAISINfit(mat,sample,design,testtype='unpaired',filtergene=F)) ## [0,1]
  res2 <- RAISINtest(RAISINfit(expr,sample,design,testtype='unpaired',filtergene=F)) ## -inf, +inf
  
  ## test other methods
  
  # res_limma01 <- fun()
  
  ind <- NULL
  expr <- mat
  res_limmainf <- fun()
  
  ## auc
  # res <- res_limma01
  # sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
  # auc <- AreaUnderSensFdr(sensfdr)
  # auc_limma01 <- auc
  # 
  res <- res_limmainf
  sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
  auc <- AreaUnderSensFdr(sensfdr)
  auc_limmainf <- auc
  
  # res <- res1
  # sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
  # auc <- AreaUnderSensFdr(sensfdr)
  # auc_our01 <- auc
  
  res <- res2
  sensfdr <- SensFdr(grep('add:',rownames(res),value=T), res)
  auc <- AreaUnderSensFdr(sensfdr)
  auc_ourinf <- auc
  
  # c(auc_limma01, auc_limmainf, auc_our01, auc_ourinf)
  
  Res <- cbind(Res, c(auc_limmainf, auc_ourinf))  
}

rownames(Res) <- c(paste0(c('auc_limmainf'), c(';fdr.diff', ';auc')), paste0(c('auc_ourinf'), c(';fdr.diff',';auc')))
saveRDS(t(Res), paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/comparelimma/res/', nSample, '_', nCell, '_', propCpG, '_', probRead, '.rds'))

