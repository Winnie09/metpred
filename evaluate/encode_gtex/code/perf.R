sm <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/wgbs.rds')

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
d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/pred/',type,'.rds'))
sm <- sm[rownames(d),colnames(d)]
sample_cor <- corfunc(t(d),t(sm))
cpg_cor <- corfunc(d,sm)
cpg_sd <- rowsds(sm)
dif <- (d-sm)^2
sample_rmse <- sqrt(colMeans(dif))
cpg_rmse <- sqrt(rowMeans(dif))

sample <- data.frame(cor=sample_cor,rmse=sample_rmse,stringsAsFactors = F)
cpg <- data.frame(cor=cpg_cor,rmse=cpg_rmse,sd=cpg_sd,stringsAsFactors = F)


saveRDS(sample,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/perf/',type,'/sample.rds'))
saveRDS(cpg,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/perf/',type,'/cpg.rds'))

cb <- combn(ncol(d),2)
difcor <- do.call(rbind,sapply(1:ncol(cb),function(k) {
  smd <- sm[,cb[1,k]]-sm[,cb[2,k]]
  dd <- d[,cb[1,k]]-d[,cb[2,k]]
  data.frame(cor=cor(dd,smd),absdif=sqrt(mean((smd-dd)^2)),sd=sd(smd),comp=paste0(colnames(sm)[cb[1,k]],':',colnames(sm)[cb[2,k]]),stringsAsFactors = F)
},simplify = F))

saveRDS(difcor,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_gtex/perf/',type,'/difcor.rds'))

