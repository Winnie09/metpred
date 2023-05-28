## ENCODE
d <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38_normal/rna.rds')
tis <- sub(' originated from.*','',sub('motor neuron .*','motor neuron',gsub(' tissue','',sub(' male.*| female.*','',sub('Homo sapiens ','',colnames(d))))))
tis[grep('esophagus|gastroesophageal',tis)] <- 'esophagus'
tis[grep('heart|atrium',tis)] <- 'heart'
tis[grep('lung',tis)] <- 'lung'
tis[grep('muscle',tis)] <- 'muscle'
tis[grep('skin',tis)] <- 'skin'
tis[grep('liver',tis)] <- 'liver'
tis[grep('colon',tis)] <- 'colon'

sex <- rep(NA,ncol(d))
sex[grep(' male ',colnames(d))] <- 'male'
sex[grep(' female ',colnames(d))] <- 'female'
age <- rep(NA,ncol(d))
age[grep('\\(',colnames(d))] <- sub('\\)','',sub('.*\\(','',colnames(d)[grep('\\(',colnames(d))]))
lf <- rep(NA,ncol(d))
lf[grep('adult \\(',colnames(d))] <- 'adult'
lf[grep('child \\(',colnames(d))] <- 'child'
lf[grep('embryo \\(',colnames(d))] <- 'embryo'
df <- data.frame(DataSource = 'ENCODE', 
                 Species = 'Homo sapiens',
                 TissueName=colnames(d),
                 GeneralTissue=tis,
                 Sex=sex,
                 Age=age,LifeStage=lf,stringsAsFactors = F)
apply(df, 2, table)
rdir <- '/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/'
saveRDS(df, paste0(rdir, 'meta_included_95_normal_samples.rds'))
write.csv(df, paste0(rdir, 'meta_included_95_normal_samples.csv'), row.names = F)
write.csv(df, '/home/whou10/data/whou/metpred/data/processed_meta/meta_included_95_normal_samples.csv', row.names = F)


## ===============================
## ============================ matched ===================================
### human GRCh38 (normal 104)
## encode
rdir1 <- '/home/whou10/data/whou/metpred/data/encode_wgbs/processed/GRCh38/' #100
rna = readRDS(paste0(rdir1, 'rna.rds'))
str(rna)


[meta]
r.m = readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/rna.rds') #361
str(r.m)
encode.r.m <- r.m[r.m[,3] %in% df[,3], ]

w.m <- readRDS('/home/whou10/data/whou/metpred/data/encode_wgbs/meta/combine/wgbs.rds') #255
str(w.m)
encode.w.m <- w.m[w.m[,3] %in% df[,3], ]

write.csv(encode.r.m, paste0(rdir, 'meta_included_95_normal_samples_rna_encff.csv'), row.names = F)
write.csv(encode.w.m, paste0(rdir, 'meta_included_95_normal_samples_wgbs_encff.csv'), row.names = F)


write.csv(encode.r.m, '/home/whou10/data/whou/metpred/data/processed_meta/meta_included_95_normal_samples_rna_encff.csv', row.names = F)
write.csv(encode.w.m, '/home/whou10/data/whou/metpred/data/processed_meta/meta_included_95_normal_samples_wgbs_encff.csv', row.names = F)





## gtex
/home/whou10/data/whou/metpred/data/gtex_wgbs/processed 9

[meta]
gtex.m <- readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/rnaseq/sample.rds') 
str(gtex.m)

r = readRDS('/home/whou10/data/whou/metpred/data/gtex_wgbs/processed/rna.rds')
wgbs.rds
str(r)
gtex.m[gtex.m[,'SMTS'] %in% colnames(r), ]
intersect(colnames(r), gtex.m[,'SMTS']) ## check with Jason


## ====================================================    
##  not used version 2 (456 samples)
## ====================================================
/home/whou10/data/whou/metpred/data/v_20220624/combine 456 (Blueprint: 116, ENCODE: 96, GEO: 244) [exclude encode]

## =======
##  mouse
## =======
### mouse mm10
/home/whou10/data/whou/metpred/data/notused/mm10 72
expr.rds
meth.rds

## =========== not used version 1 (arxive) ===========
## =======
##  human 
## =======

## CEEHRC 
/home/whou10/data/whou/metpred/data/v_2021_notused/CEEHRC_wgbs/processed 92
rna <- readRDS('/home/whou10/data/whou/metpred/data/v_2021_notused/CEEHRC_wgbs/processed/rna.rds')
wgbs <- readRDS('/home/whou10/data/whou/metpred/data/v_2021_notused/CEEHRC_wgbs/processed/wgbs.rds')
str(rna)

m_ceehrc <- read.table('/home/whou10/data/whou/metpred/data/v_2021_notused/CEEHRC_wgbs/meta/meta.txt',
                       sep = '\t',
                       header = T) #107
str(m_ceehrc)

