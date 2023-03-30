d <- readLines('https://epigenomesportal.ca/api/datahub/download?build=2019-11&assembly=4&format=text')
rna <- grep('rpkm',d,value=T)
rnasamp <- unique(sub('.*\\.','',sub('\\.mRNA-Seq.*','',rna)))
wgbs <- grep('CEEHRC',grep('WGB-Seq.methylation_profile',d,value=T),value=T)
wgbssamp <- unique(sub('.*\\.','',sub('\\.WGB-Seq.*','',wgbs)))
writeLines(c('cd /home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/rna',paste0('wget ',rna)),'/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/code/down/rna.sh')
wgbs <- wgbs[wgbssamp %in% rnasamp]
writeLines(c('cd /home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/down/wgbs',paste0('wget ',wgbs)),'/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/code/down/wgbs.sh')
