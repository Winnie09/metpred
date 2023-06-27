source('/home/whou10/scratch16/whou10/resource/startup.R')
cc3 <- lapply(seq(1, 5), function(cvid){
  dir <- paste0('evaluate/mousepancreas/celltypespecific_cv/perf/pcc/cv', cvid, '/')
  allf <- list.files(dir)
  cc <- lapply(allf, function(f){
    tmp <- readRDS(paste0(dir, f))
    chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
    tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y')), ]
    set.seed(12345)
    tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
    tmp2 <- data.frame(cor = tmp,  
                       sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
  })
  cc2 <- do.call(rbind, cc)
})
cc4 <- do.call(rbind, cc3)  
acrossRowCor_plot(plotdata = cc4,
                  ylab = 'across-sample PCC',
                  savefile = TRUE,
                  filename = 'evaluate/mousepancreas/celltypespecific_cv/plot/across_sample_pcc.pdf')
