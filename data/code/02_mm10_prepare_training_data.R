meth <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/mm10filtermat/mat.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/TPM.rds')
colnames(expr) <- sub('.tsv', '', colnames(expr))
mm10experiment <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/meta/mm10_experiment.rds')
mm10file <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/meta/mm10_file.rds')
mm10experiment[, 'Biosample summary'] <- gsub(' ', '_', mm10experiment[, 'Biosample summary'])
acc <- mm10experiment[match(colnames(meth), mm10experiment[, 'Biosample summary']), 'Accession']
facc <- mm10file[match(acc, mm10file[,'Experiment accession']), 'File accession']
expr <- expr[, facc]
colnames(expr) <- colnames(meth)
saveRDS(meth, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')
saveRDS(expr, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')
d <- meth
cm <- rowMeans(d)
csv <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
d <- d[csv >= 0.2,]
set.seed(12345)
d <- d[sample(1:nrow(d),10000),]
saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/sampleme.rds')


