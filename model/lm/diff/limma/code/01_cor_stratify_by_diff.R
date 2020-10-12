id = as.numeric(commandArgs(trailingOnly = T)[[1]])
allf = list.files('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geonly/')
f = allf[id]
load(paste0('/home-4/zji4@jhu.edu/scratch/metpred/model/lm/cv/res/lm_genecluster1000_geonly/',f))
proj <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/tcga/proc/combine/project.rds')
proj = gsub('TCGA-','',proj)
sproj = proj[colnames(testm)]
cc <- lapply(1:(length(unique(sproj))-1), function(ii){
  i = unique(sproj)[ii]
  t(sapply((ii+1):length(unique(sproj)), function(jj){
    j = unique(sproj)[jj]  
    u <- rowMeans(pred[, sproj==i]) - rowMeans(pred[, sproj==j])
    v <- rowMeans(testm[, sproj==i]) - rowMeans(testm[, sproj==j])
    a <- cor(u,v)
    print(paste0(i,'_',j))
    expr = testm[, sproj %in% c(i,j)]
    ssproj = proj[colnames(expr)]
    library(limma)
    des <- cbind(1,ifelse(ssproj==i,1,0))
    fit <- eBayes(lmFit(voom(expr,des),design=des))
    res <- topTable(fit,coef=2,number=nrow(expr))
    res <- res[res[,'adj.P.Val']<0.05,]
    gs <- rownames(res)
    b <- cor(u[gs],v[gs])
    c <- cor(u[setdiff(names(u), gs)], v[setdiff(names(v), gs)])
    c(a, b, c)
  }))
})
cc = do.call(rbind, cc)
colnames(cc) <- c('all_loci', 'diff_loci', 'nondiff_loci')
saveRDS(cc,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/diff/limma/tmpres/',f,'.rds'))

# smoothScatter(u, v, xlab='predicted difference', ylab='measured difference', main=paste0(i, ',', j, ', PCC=', round(cor(u,v),2)))

