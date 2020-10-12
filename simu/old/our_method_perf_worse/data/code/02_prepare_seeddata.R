nSample <- 4
nCell <- 2
propCpG <- 0.2
probRead <- 0.6

fullmat <- readRDS(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/nullsimu/data/data/', nSample,'To', nSample, '_allcell/pred_ery.rds'))
ap <- as.numeric(sub('BM','',sub(':.*', '', colnames(fullmat))))
as <- sub('_.*','', colnames(fullmat))
ind <- NULL

logit <- function(x) log(x/(1-x))

res <- sapply(1:1e4, function(se){
  print(se)
  mat <- lapply(unique(as), function(i){
    set.seed(se)
    id <- sample(which(as == i), nCell, replace = TRUE)
    fullmat[, id]
  })
  mat <- do.call(cbind, mat)
  colnames(mat) <- paste0(colnames(mat), '_', 1:ncol(mat))
  
  sample <- sub(':.*','',colnames(mat))
  samplename <- unique(sample)
  sum(as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 0 == 1) < length(samplename)
  group <- as.numeric(sub('BM','',samplename)) %in% c(1,2,5,6) + 1
  
  expr <- logit(mat)
  
  
  suppressMessages(library(limma))
  if (is.null(ind) || sum(duplicated(ind))==0) {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(group==group[1]))
  } else {
    design <- cbind(Grp1=1,Grp2vs1=as.numeric(group==group[1]),ind=model.matrix(~x-1,data.frame(x=ind)))
  }
  sampsum <- sapply(samplename,function(us) rowMeans(expr[,sample==us,drop=F]))
  options(digits=3)
  fit <- lmFit(sampsum, design=design)
  fit <- eBayes(fit)
  res <- topTable(fit,coef=2,number=nrow(expr))
  colnames(res)[colnames(res)=='P.Value'] <- "pvalue"
  colnames(res)[colnames(res)=='t'] <- "stat"
  colnames(res)[which(colnames(res) == 'adj.P.Val')] <- 'FDR'
  res <- res[order(res[,'pvalue'],-abs(res[,'stat'])),]
  res[,5]
})

rownames(res) <- rownames(fullmat)
saveRDS(res, '/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/simu/data/seeddata/res.rds')

