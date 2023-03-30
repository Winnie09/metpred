d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/burl.rds')
writeLines(c('cd /home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/b',paste0('wget ',d)),'/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/code/downb.sh')
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/rurl.rds')
writeLines(c('cd /home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/down/r',paste0('wget ',d)),'/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/code/downr.sh')
