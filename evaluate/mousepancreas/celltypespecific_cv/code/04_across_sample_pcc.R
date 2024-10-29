method = commandArgs(trailingOnly = T)[[1]] # method = 'permu'
source('/home/whou10/scratch16/whou10/resource/startup.R')
rdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/perf/all/'
cvid = 1
cc3 <- lapply(seq(1, 5), function(cvid){
  print(cvid)
  dir <- paste0('evaluate/mousepancreas/celltypespecific_cv/perf/', method, '/pcc/cv', cvid, '/')
  allf <- list.files(dir)
  cc <- lapply(allf, function(f){
    tmp <- readRDS(paste0(dir, f))
    chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
    tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y'))]
    set.seed(12345)
    tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
    tmp2 <- data.frame(cor = tmp,  
                       sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
  })
  cc2 <- do.call(rbind, cc)
})
cc4 <- do.call(rbind, cc3)  
saveRDS(cc4, paste0(rdir, method, '_across_sample_pcc.rds'))
