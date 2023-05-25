d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
tis <- sub(' originated from.*','',sub('motor neuron .*','motor neuron',gsub(' tissue','',sub(' male.*| female.*','',sub('Homo sapiens ','',colnames(d))))))
tis[grep('esophagus|gastroesophageal',tis)] <- 'esophagus'
tis[grep('heart|atrium',tis)] <- 'heart'
tis[grep('lung',tis)] <- 'lung'
tis[grep('muscle',tis)] <- 'muscle'
tis[grep('skin',tis)] <- 'skin'
tis[grep('liver',tis)] <- 'liver'
tis[grep('colon',tis)] <- 'colon'
tab <- table(tis)
pd <- data.frame(n=names(tab),co=as.vector(tab),stringsAsFactors = F)
pd[,1] <- factor(pd[,1],levels=pd[order(pd[,2]),1])
library(ggplot2)
pdf('/home/whou10/data/whou/metpred/data/encode_wgbs/plot/tissue.pdf')
ggplot(pd,aes(x=co,y=n)) + geom_bar(stat='identity',fill='black',alpha=0.7) + theme_classic()
dev.off()
