m <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38filtermat/mat.rds')
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/wgbs/hg38filtermat/gr.rds')
end(gr)=start(gr)
row.names(m)=sub(':','_',as.character(gr))
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/encode_wgbs/data/wgbs.rds')
