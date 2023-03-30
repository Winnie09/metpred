# pred <- readRDS('/home/whou10/scratch16/whou10/metpred/predict/res/predicted_DNAm_on_testset_cancer.rds')
pred <- readRDS('/home/whou10/scratch16/whou10/metpred/predict/res/predicted_DNAm_on_testset_normal.rds')
# testid <- readRDS('/home/whou10/scratch16/whou10/metpred/predict/res/testid_cancer.rds')
testid <- readRDS('/home/whou10/scratch16/whou10/metpred/predict/res/testid_normal.rds')
m.bak <- readRDS('/home/whou10/scratch16/whou10/metpred/final/nonawgbs_hg38.rds')
true = m.bak[, testid]

ct = readRDS('/home/whou10/scratch16/whou10/metpred/quality/ct.rds')
table(ct[testid])
ct.uni = unique(ct[testid])
ct.uni = ct.uni[!is.na(ct.uni)]
ct1 = ct.uni[1]

for (ct1id in 1:(length(ct.uni)-1)){
  ct1 = ct.uni[ct1id]
  print(ct1)
  print('\n')
  for (ct2id in ct1id:length(ct.uni)){
    ct2 = ct.uni[ct2id]
    print(ct2)
    s1 <- names(which(ct[testid]==ct1))
    s2 <- names(which(ct[testid]==ct2))
    
    pred1 <- pred[,s1]
    pred2 <- pred[,s2]
    true1 <- true[,s1]
    true2 <- true[,s2]
    
    set.seed(12345)
    id = sample(1:nrow(pred), 1e4)
    
    predpval <- sapply(id,function(i) {
      wilcox.test(pred1[i,],pred2[i,])$p.value
    })
    names(predpval) <- rownames(pred1)[id]
    
    truepval <- sapply(id,function(i) {
      wilcox.test(true1[i,],true[i,])$p.value
    })
    names(truepval) <- rownames(true1)[id]
    predpval.sort = sort(predpval)
    truepval.sort = sort(truepval)
    o = sapply(1:(length(id)/10),function(j) {
      length(intersect(names(predpval.sort)[1:(10*j)], names(truepval.sort)[1:(10*j)]))/(10*j)
    })
    # saveRDS(o, '/home/whou10/scratch16/whou10/metpred/diffDNAm/plot/cancer_randomly1e4CpG_overlap.rds')
    saveRDS(o, paste0('/home/whou10/scratch16/whou10/metpred/diffDNAm/plot/normal_randomly1e4CpG_overlap_', ct1, '_', ct2,'.rds'))
  }
}

