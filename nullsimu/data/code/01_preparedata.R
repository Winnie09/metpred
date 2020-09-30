setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/')
# read data: WGBS_hg38 model, HCA imputed gene expression
model <- readRDS('./model/wgbs_hg38/model.rds')
data <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/matrix/saver.rds')

# subset cluster 2 (all erythroids) 
clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/hca/data/proc/cluster/cluster.rds')
clu <- clu[colnames(data)]
selectcell <- names(clu[clu==2])  # all erythroid
table(clu[selectcell])
table(ct[selectcell])
expr <- data[, selectcell]

# use WGBS prediction model but predict only on 1e4 CpG
source('./software/predict.R')
set.seed(12345)
pred <- predict(expr = expr, model = model)
saveRDS(pred, './nullsimu/data/data/4To4_allcell/pred_ery.rds')

# 2 to 2 samples
ap <- paste0(sapply(colnames(expr), function(i) strsplit(i, ':')[[1]][1]),
             ':',
             sapply(colnames(expr), function(i) strsplit(i, ':')[[1]][2]))

ap <- sub('_.*', '', colnames(pred))
group1 <- c('BM5:29:male', 'BM6:26:female' )
group2 <- c('BM4:29:male', 'BM8:32:female') 
pred2 <- pred[, ap %in% c(group1, group2)]
saveRDS(pred2, './nullsimu/data/data/2To2_allcell/pred_ery.rds')

# subsample 2 cells each sample
set.seed(12345)
id <- c(sample(which(ap == group1[1]), 2), 
sample(which(ap == group1[2]), 2),
sample(which(ap == group2[1]), 2),
sample(which(ap == group2[2]), 2))
pred3 <- pred[,id]
saveRDS(pred3, './nullsimu/data/data/2To2_2cell/pred_ery.rds')



