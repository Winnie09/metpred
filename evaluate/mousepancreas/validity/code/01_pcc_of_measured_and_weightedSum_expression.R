library(Matrix)
setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
source('/home/whou10/data/whou/metpred/software/trainpredict.R')
r <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/cse.rds')
w <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/wgbs/bs.rds')
int <- intersect(colnames(r),colnames(w))
samp <- sub(':.*','',int)
alls <- sub(':.*', '', colnames(r))
alls.w <- sub(':.*', '', colnames(w))
allct.w <- sub('.*:', '', colnames(w))
pdlist <- list()

for (i in sub('.rds', '',list.files('goldstandard/res/splotList_weightedSum_0.95_cpg_sd_top1e4'))[12]){
  print(i)
  trainid <- int[which(samp != i)]
  testid <- setdiff(int, trainid)
  cse.s <- r[, alls == i]
  allct <- sub('.*:', '', colnames(cse.s))
  ct.s <- allct.w[alls.w == i]
  
  e <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/normexpr.rds')
  se <- sub(':.*','',colnames(e))
  e <- e[,se == i]
  
  ## group the CpGs from its standard deviation in testing samples smallest to largets slots
  ## run the prediction for each CpG group
  rs <- rowsds(w[,testid])
  rs2 <- cut(rs,seq(0,1,0.05))
  names(rs2) <- names(rs)
  
  dir.gs <- 'goldstandard/res/splotList_weightedSum_0.95_cpg_sd_top1e4/'
  gs <- readRDS(paste0(dir.gs, i, '.rds'))
  
  ##  compare the weighted average of cell-type specific gene expression 
  ## and the observed spot-level gene expression for first sample
  weight <- readRDS('/home/whou10/data/whou/metpred/data/mousepancreas/spatial/weight.rds')
  weight <- weight[sub(':.*', '', rownames(weight)) == i, ]
  
  ### weighted averaged expression
  wae <- cse.s %*% t(weight[,allct])
  str(wae)
  e <- e[rownames(wae),]
  
  ### all spots correlation
  cc <- corfunc(t(e), t(wae))
  saveRDS(cc, paste0('validity/res/', i, '_spot_measure_wae_cor.rds'))
  
  # Bin size control + color palette
  p1 <- ggplot(data.frame(x = e[,3], y = wae[,3]), aes(x=x, y=y) ) +
    geom_bin2d(bins = 50) +
    # geom_point(alpha = 0.3)+
    scale_fill_continuous(type = "viridis") +
    xlab('Measured expression') + 
    ylab('WeightedSum expression') + 
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.3)
  p2 <- ggplot(data.frame(x = e[,colnames(gs)[10]], y = wae[,colnames(gs)[10]]), aes(x=x, y=y) ) +
    geom_bin2d(bins = 50) +
    # geom_point(alpha = 0.3)+
    scale_fill_continuous(type = "viridis") +
    xlab('Measured expression') + 
    ylab('WeightedSum expression')+
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.3)
  p3 <- ggplot(data.frame(x = e[,colnames(gs)[30]], y = wae[,colnames(gs)[30]]), aes(x=x, y=y) ) +
    geom_bin2d(bins = 50) +
    # geom_point(alpha = 0.3)+
    scale_fill_continuous(type = "viridis") +
    xlab('Measured expression') + 
    ylab('WeightedSum expression')+
    geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.3)
  plist <- list(p1,p2,p3)
  library(gridExtra)
  pdir <- 'validity/plot/'
  dir.create(pdir)
  pdf(paste0(pdir, i, '_example_spot.pdf'), width = 5.5, height = 1.5)
  grid.arrange(grobs = plist, nrow = 1)
  dev.off()
  
  pdlist[[i]] <- data.frame(weightSum = rowSums(weight[, ct.s]), 
                            pcc = cc, 
                            sample = i, 
                            stringsAsFactors = F)
  
}

pd <- do.call(rbind, pdlist)
pdf(paste0(pdir, 'pcc_measured_and_weightedSum_expr.pdf'), width = 5, height = 4)
ggplot(pd, aes(x=weightSum, y=pcc) ) +
  geom_bin2d(bins = 80) +
  scale_fill_continuous(type = "viridis") +
  xlab('Sum of deconvoluted weights') + 
  ylab('PCC of measured and weighted sum expression') +
  facet_wrap(~sample, ncol = 4 ) 
dev.off()

