library(TSCAN)
d <- readRDS('/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/perf/ramp/pr.rds')
d <- d$x
d <- d[,1:10]
d <- d[sub(':.*','',rownames(d)) %in% c('naive_t','cytotoxic_t'),]
ct <- sub(':.*','',rownames(d))
names(ct) <- rownames(d)
ord <- guided_tscan(d,ct,c('naive_t','cytotoxic_t'))
saveRDS(ord,file='/home/whou10/data/whou/metpred/evaluate/pbmc/predsc/imp/res/pseudotime/ord.rds')
