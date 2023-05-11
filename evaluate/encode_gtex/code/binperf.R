sm <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/wgbs_bin.rds')

rowsds <- function(data) {
  cm <- rowMeans(data)
  sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
}

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}

type <- commandArgs(trailingOnly = T)
d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/binpred/',type,'.rds'))
int <- intersect(rownames(sm),rownames(d))
sm <- sm[int,colnames(d)]
d <- d[int,]
sample_cor <- corfunc(t(d),t(sm))
cpg_cor <- corfunc(d,sm)
cpg_sd <- rowsds(sm)
dif <- (d-sm)^2
sample_rmse <- sqrt(colMeans(dif))
cpg_rmse <- sqrt(rowMeans(dif))

sample <- data.frame(cor=sample_cor,rmse=sample_rmse,stringsAsFactors = F)
cpg <- data.frame(cor=cpg_cor,rmse=cpg_rmse,sd=cpg_sd,stringsAsFactors = F)


saveRDS(sample,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/binperf/',type,'/sample.rds'))
saveRDS(cpg,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/binperf/',type,'/cpg.rds'))


