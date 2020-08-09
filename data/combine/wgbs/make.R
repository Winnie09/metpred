e1 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/ge.rds')
e1 <- log2(e1 + 1)
ef <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/rurl.rds')
ef <- sub('.*/','',ef)
m1 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/proc/filterme.rds')
mf <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/blueprint_wgbs/data/filelist/burl.rds')
mf <- sub('.*/','',mf)
colnames(e1) <- paste0('blueprint:',gsub(':','-',names(ef)[match(colnames(e1),ef)]))
colnames(m1) <- paste0('blueprint:',gsub(':','-',names(mf)[match(colnames(m1),mf)]))

e2 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/ge.rds')
colnames(e2) <- paste0('CEEHRC:',sub('.rds','',colnames(e2)))
m2 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/CEEHRC_wgbs/data/proc/filterme.rds')
m2 <- m2/100
colnames(m2) <- paste0('CEEHRC:',sub('.*\\.','',sub('.WGB-Seq.methylation_profile.bigWig','',colnames(m2))))

e3 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/encode_wgbs/data/expr.rds')
m3 <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/encode_wgbs/data/wgbs.rds')
colnames(e3) <- paste0('encode:',colnames(e3))
colnames(m3) <- paste0('encode:',colnames(m3))

print(summary(apply(e1,2,max)))
print(summary(apply(e2,2,max)))
print(summary(apply(e3,2,max)))
print(summary(apply(m1,2,max,na.rm=T)))
print(summary(apply(m2,2,max,na.rm=T)))
print(summary(apply(m3,2,max,na.rm=T)))

g <- unique(c(row.names(e1),row.names(e2),row.names(e3)))
g2 <- unique(c(row.names(m1),row.names(m2),row.names(m3)))
n <- unique(c(colnames(e1),colnames(e2),colnames(e3)))
n2 <- unique(c(colnames(m1),colnames(m2),colnames(m3)))

m <- matrix(NA,nrow=length(g2),ncol=length(n),dimnames=list(g2,n))
e <- matrix(NA,nrow=length(g),ncol=length(n),dimnames=list(g,n))
m[rownames(m1),colnames(m1)] <- m1
m[rownames(m2),colnames(m2)] <- m2
m[rownames(m3),colnames(m3)] <- m3

e[rownames(e1),colnames(e1)] <- e1
e[rownames(e2),colnames(e2)] <- e2
e[rownames(e3),colnames(e3)] <- e3
colnames(m) <- gsub(',','',colnames(m))
colnames(e) <- gsub(',','',colnames(e))

saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/me.rds')
saveRDS(e,file='/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/ge.rds')

m <- m[rowMeans(is.na(m))==0,]
e <- e[rowMeans(is.na(e))==0,]
saveRDS(m,file='/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/filterme.rds')
saveRDS(e,file='/home-4/zji4@jhu.edu/scratch/metpred/data/combine/wgbs/filterge.rds')

