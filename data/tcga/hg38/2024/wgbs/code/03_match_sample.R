ge_meta = readRDS('/home/whou10/data/whou/metpred/data/tcga/hg38/450k/filelist/ge.rds') ## 10499    33
w_meta = read.table('/home/whou10/data/whou/metpred/data/tcga/hg38/wgbs/meta/tcga_wgbs_filelist.txt') ## 47  4

ge_pt = sapply(ge_meta$patient, function(i) strsplit(i, '-')[[1]][3]) 
w_pt = sapply(sub('.bed.gz','',w_meta[,4]), function(i){
  strsplit(i, '_')[[1]][3]
})
unname(w_pt)
int = intersect(w_pt, ge_pt)
str(int)

diff = setdiff(w_pt, ge_pt)
str(diff)

diff_sp = sapply(diff, function(i){
  sub('.bed.gz','',w_meta[which(w_pt == i), 4])
})

