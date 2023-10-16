setwd('/home/whou10/data/whou/metpred/')
source('/home/whou10/scratch16/whou10/resource/startup.R')
dir <- 'evaluate/mousepancreas/trainDuctAcinar_predADM/perf/pcc/'
allf <- list.files(dir)
cc <- lapply(allf, function(f){
  print(f)
  tmp <- readRDS(paste0(dir, f))
  chr <- sub(':.*', '', sub('_.*', '', names(tmp)))
  tmp <- tmp[chr %in% paste0('chr', c(seq(1,22), 'X', 'Y'))]
  set.seed(12345)
  tmp <- sample(tmp, min(1e2, length(tmp)), replace = F)
  tmp2 <- data.frame(cor = tmp,  
                     sd = paste0('[',sub('.rds', '', sub('testsd_cpg_', '', f)),')'))
})
cc2 <- do.call(rbind, cc)
acrossRowCor_plot(plotdata = cc2,
                  ylab = 'across-sample PCC',
                  savefile = TRUE,
                  filename = 'evaluate/mousepancreas/trainDuctAcinar_predADM/plot/across_sample_pcc.pdf')


