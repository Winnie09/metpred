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


sfunc <- function(tmpm) {
  cm <- rowMeans(tmpm) ## CpG site mean
  csd <- sqrt((rowMeans(tmpm*tmpm) - cm^2) / (ncol(tmpm) - 1) * ncol(tmpm)) ## CpG site sd
  names(head(sort(csd,decreasing = T),100000))  
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
s <- 'A46_F_P_K4_D4'

gs <- readRDS(paste0(dir.gs, s, '.rds'))
str(gs)
pred <- readRDS(paste0(dir.pd, s, '.rds'))
str(pred)
pred <- pred[, colnames(gs)]
length(intersect(rownames(pred), rownames(gs)))
cc <- corfunc(pred, gs)

saveRDS(cc, paste0(dir.pf, s, '.rds'))


mse <- sapply(1:nrow(pred), function(i){
  tmpv <- pred[i, ] - gs[i, ]
  mean(tmpv * tmpv)
})
names(mse) <- rownames(pred)

saveRDS(mse, paste0(dir.pf2, s, '.rds'))


## rethink:save a vector for testing set CpG sd, or otherwise need to read w all the time

### plot
rs <- rowsds(gs)
rs2 <- cut(rs, seq(0, 1, 0.05))
names(rs2) <- names(rs)
cc <- cc[names(rs2)]



pdf(paste0(pdir, s ,'_pcc.pdf'),width=2.5,height=2.5)
print(ggplot(data.frame(sd=rs2,cor=cc),aes(y=cor,x=sd,fill=rs2)) + 
        geom_violin(alpha=0.2, scale = 'area') + 
        geom_boxplot(alpha=0.3, width = 0.2) + 
        theme_classic() + 
        scale_fill_manual(values=rainbow(length(unique(cc)))) + 
        xlab('Testing measured value sd') + 
        ylab('PCC') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

pdf(paste0(pdir, s ,'_mse.pdf'),width=2.5,height=2.5)
print(ggplot(data.frame(sd=rs2,mse=mse),aes(y=mse,x=sd,fill=rs2)) + 
        geom_violin(alpha=0.2, scale = 'area') + 
        geom_boxplot(alpha=0.3, width = 0.2) + 
        theme_classic() + 
        scale_fill_manual(values=rainbow(length(unique(cc)))) + 
        xlab('Testing measured value sd') + 
        ylab('MSE') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()


