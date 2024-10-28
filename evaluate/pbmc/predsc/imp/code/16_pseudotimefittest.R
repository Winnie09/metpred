library(PreTSA)
ord <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/ord.rds')
af <- list.files('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred.bak.clu1e3.lambda10e-1/allcpg/')
ad <- NULL
for (f in af) {
  d <- readRDS(paste0('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pred.bak.clu1e3.lambda10e-1/allcpg/',f))
  ct <- sub(':.*','',colnames(d))
  d <- d[,ord]
  
  mysd <- function(mat){
    ## calculate the sd for mat rows.
    m <- rowMeans(mat)
    s <- sqrt((rowMeans(mat*mat) - m^2) * ncol(mat)/(ncol(mat)-1))
  }
  sdv <- mysd(d)
  
  d <- d[sdv > 0.05,]
  ad <- rbind(ad,d)
  rm('d')
}
ad <- ad[complete.cases(ad),]
pt <- 1:ncol(ad)
names(pt) <- colnames(ad)
r <- temporalFit(expr = ad, pseudotime = pt, knot = F)
saveRDS(r,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/mefit.rds')
r <- temporalTest(expr = ad, pseudotime = pt, knot = F)
saveRDS(r,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/metest.rds')

ddir <- '/home/whou10/data/whou/metpred/data/pbmc/imp/res/'
ad = readRDS(paste0(ddir, 'rna_test_sc_sub300.rds'))
ad = ad[,ord]
pt <- 1:ncol(ad)
names(pt) <- colnames(ad)
r <- temporalFit(expr = ad, pseudotime = pt, knot = F)
saveRDS(r,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/gefit.rds')
r <- temporalTest(expr = ad, pseudotime = pt, knot = F)
saveRDS(r,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/getest.rds')

