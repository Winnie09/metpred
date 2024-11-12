me <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/wgbs.rds')

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

sid <- 10

d <- sapply(1:sid,function(i) {
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/pred/',type,'/',i,'.rds'))
  sm <- me[rownames(d),colnames(d)]
  sample_cor <- corfunc(t(d),t(sm))
  cpg_cor <- corfunc(d,sm)
  cpg_sd <- rowsds(sm)
  dif <- (d-sm)^2
  sample_rmse <- sqrt(colMeans(dif))
  cpg_rmse <- sqrt(rowMeans(dif))
  list(sample_cor=sample_cor,cpg_cor=cpg_cor,cpg_sd=cpg_sd,sample_rmse=sample_rmse,cpg_rmse=cpg_rmse)
},simplify = F)

sample <- cpg <- NULL
for (i in 1:sid) {
  sample <- rbind(sample,data.frame(cor=d[[i]]$sample_cor,rmse=d[[i]]$sample_rmse,cvid=i,stringsAsFactors = F))
  cpg <- rbind(cpg,data.frame(cor=d[[i]]$cpg_cor,rmse=d[[i]]$cpg_rmse,sd=d[[i]]$cpg_sd,cvid=i,stringsAsFactors = F))
}

saveRDS(sample,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/',type,'/sample.rds'))
saveRDS(cpg,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/',type,'/cpg.rds'))

difcor <- do.call(rbind,sapply(1:sid,function(i) {
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/pred/',type,'/',i,'.rds'))
  sm <- me[rownames(d),colnames(d)]
  cb <- combn(ncol(d),2)
  do.call(rbind,sapply(1:ncol(cb),function(k) {
    smd <- sm[,cb[1,k]]-sm[,cb[2,k]]
    dd <- d[,cb[1,k]]-d[,cb[2,k]]
    data.frame(cor=cor(dd,smd),absdif=sqrt(mean((smd-dd)^2)),sd=sd(smd),comp=paste0(colnames(sm)[cb[1,k]],':',colnames(sm)[cb[2,k]]),cvid=i,stringsAsFactors = F)
  },simplify = F))
},simplify = F))
saveRDS(difcor,file=paste0('/home/whou10/data/whou/metpred/evaluate/encode_cv/perf/',type,'/difcor.rds'))


