setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/meta/dnam')
a = read.csv('./encode_wgbs/match_experiment.csv',as.is=T,header=T)
a = a[,-1]
ct1 = a$Biosample.term.name
b = read.csv('./geo_wgbs/wgbstable.csv',as.is=T,header=T)
b = b[,-1]
b = b[b$Organism == 'Homo sapiens',]
b$Source = sub('; WGBS.*','',b$Source)
b$Source = sub(', WGBS.*','',b$Source)
ct2 = b$Source
length(union(unique(ct1),unique(ct2)))
write.csv(unique(ct1), './encode_wgbs/unique_ct.csv')
write.csv(unique(ct2), './geo_wgbs/unique_ct.csv')
