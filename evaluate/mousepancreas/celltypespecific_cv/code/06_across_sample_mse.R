method = commandArgs(trailingOnly = T)[[1]] # method = 'permu'
print(method)
rdir <- '/home/whou10/data/whou/metpred/evaluate/mousepancreas/celltypespecific_cv/perf/all/'
setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
mse <- lapply(seq(1, 5), function(cvid){
  print(cvid)
  dir <- paste0('evaluate/mousepancreas/celltypespecific_cv/perf/', method, '/mse/cv', cvid, '/')
  allf <- list.files(dir)
  m <- lapply(allf, function(f){
    tmp <- readRDS(paste0(dir, f))
    print(str(tmp))
    Sys.sleep(2)
    chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
    tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y'))]
    set.seed(12345)
    tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
    tmp2 <- data.frame(mse = tmp,  
                       sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
    print(str(tmp2))
    Sys.sleep(2)
    return(tmp2)
  })
  m2 <- do.call(rbind, m)
  return(m2)
})
mse2 <- do.call(rbind, mse)  
print(mse2)
saveRDS(mse2, paste0(rdir, method, '_across_sample_mse.rds'))
# acrossRowCor_plot(plotdata = mse2,
#                   ylab = 'across-sample MSE',
#                   savefile = TRUE,
#                   filename = 'evaluate/mousepancreas/celltypespecific_cv/plot/across_sample_mse.pdf')
# 
