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


setwd('/home/whou10/data/whou/metpred/evaluate/mousepancreas/')
dir.pf <- 'perf/res/top_1e4_cpg/pcc/'
dir.pf2 <- 'perf/res/top_1e4_cpg/mse/'
dir.create(dir.pf, showWarnings = F, recursive = T)
dir.create(dir.pf2, showWarnings = F, recursive = T)
pdir <- 'perf/plot/top_1e4_cpg/'
dir.create(pdir, showWarnings = F, recursive = T)

dir.gs <- 'goldstandard/res/splotList_weightedSum_0.95_cpg_sd_top1e4/'
dir.pd <- 'pred/ramp_subset_top_1e4_sd_cpg/'

alls <- sub('.rds','',intersect(list.files(dir.gs), list.files(dir.pd)))
cclist <- mselist <- rslist <- list()

# s <- 'A46_F_P_K4_D4'
for (s in alls){
  print(s)
  gs <- readRDS(paste0(dir.gs, s, '.rds'))
  pred <- readRDS(paste0(dir.pd, s, '.rds'))
  pred <- pred[, colnames(gs)]
  length(intersect(rownames(pred), rownames(gs)))
  cc <- corfunc(pred, gs)
  dir.create(paste0(dir.pf, s), recursive = T)
  saveRDS(cc, paste0(dir.pf, s, '.rds'))
  cclist[[s]] <- cc
  
  mse <- sapply(1:nrow(pred), function(i){
    tmpv <- pred[i, ] - gs[i, ]
    mean(tmpv * tmpv)
  })
  names(mse) <- rownames(pred)
  dir.create(paste0(dir.pf2, s), recursive = T)
  saveRDS(mse, paste0(dir.pf2, s, '.rds'))
  mselist[[s]] <- mse
  
  ## rethink:save a vector for testing set CpG sd, or otherwise need to read w all the time
  
  ### plot
  rs <- rowsds(gs)
  rs2 <- cut(rs, seq(0, 1, 0.05))
  names(rs2) <- names(rs)
  cc <- cc[names(rs2)]
  rslist[[s]] <- unique(rs2)
  
}


str(mse)
str(mselist)
str(rslist)

sapply(rslist, unique)

mse <- do.call(cbind, mselist)
mse <- reshape2::melt(mse)
str(mse)
apply(mse, 2, median)

cc <- do.call(cbind, cclist)
cc <- reshape2::melt(cc)
str(cc)

pdf(paste0(pdir, 'allsamples_pcc.pdf'),width=5,height=3)
print(ggplot(cc,aes(y=value,x=Var2,fill=Var2)) + 
        geom_violin(alpha=0.2, scale = 'area') + 
        geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.1) + 
        theme_classic() + 
        scale_fill_brewer(palette = 'Set2') + 
        xlab('Samples') + 
        ylab('PCC') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

pdf(paste0(pdir, 'allsamples_mse.pdf'),width=5,height=3)
print(ggplot(mse,aes(y=value,x=Var2,fill=Var2)) + 
        geom_violin(alpha=0.2, scale = 'area') + 
        geom_boxplot(alpha=0.3, width = 0.2, outlier.size = 0.1) + 
        theme_classic() + 
        scale_fill_brewer(palette = 'Set2') + 
        xlab('Samples') + 
        ylab('MSE') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

