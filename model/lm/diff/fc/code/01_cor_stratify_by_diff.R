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
    a <- cor(u, v)
    gs1 <- which(v <= quantile(v)[2])
    gs2 <- which(v > quantile(v)[2] & v <= quantile(v)[3])
    gs3 <- which(v > quantile(v)[3] & v <= quantile(v)[4])
    gs4 <- which(v > quantile(v)[4] & v <= quantile(v)[5])
    c(cor(u,v), cor(u[gs1],v[gs1]), cor(u[gs2],v[gs2]), cor(u[gs3],v[gs3]), cor(u[gs4],v[gs4]))
  }))
})
cc = do.call(rbind, cc)
colnames(cc) <- c('all_loci', 'loci_quantile1', 'loci_quantile2', 'loci_quantile3', 'loci_quantile4')
saveRDS(cc,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/lm/diff/fc/tmpres/',f,'.rds'))
