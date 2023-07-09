library(data.table)
d <- fread('/home/whou10/data/whou/metpred/data/hsc/raw/rna/GSE87195_rnaseq_ensT_all.csv.gz',data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
w <- readRDS('/home/whou10/data/whou/metpred/data/hsc/proc/wgbs.rds')

rnact <- paste0(sub('.*_','',sub('_100','',colnames(d))),'_',sub('_.*','',sub('RNA_','',colnames(d))))
wct <- sub('_R[12].*','',sub('_1000','',sub('GSM[0-9]*_','',colnames(w))))
colnames(d) <- rnact
colnames(w) <- wct

## ranct cannot match to w in ct_donor: "HSC_D1"  "MLP2_D1" "MPP_D1"  "MPP_D2" 
## w: 
# training: "HSC_D10"  "HSC_D7" --> 'HSC_cb'; 
# "MPP_D6"   "MPP_D8"  "MPP_D7"   "MPP_D9"   "MPP_D10" --> 'MPP_cb'
# testing: "HSC_D5"   "HSC_D6" --> 'HSC_cb';
#  "MPP_D4"   "MPP_D5" --> 'MPP_cb'
## rna:
# training: "HSC_D1" --> "HSC_cb"; 
# "MPP_D1"  "MPP_D2"  --> "MPP_cb"

cw = colnames(w)
cr = colnames(d)
cw[cw %in% c("HSC_D10", "HSC_D7")] <- 'HSC_cb'
cw[cw %in% c("MPP_D6", "MPP_D8", "MPP_D7", "MPP_D9", "MPP_D10")] <- 'MPP_cb'
cw[cw %in% c("HSC_D5", "HSC_D6")] = 'HSC_cb'
cw[cw %in% c("MPP_D4", "MPP_D5")] = 'MPP_cb'
cr[cr %in% c("HSC_D1")] = "HSC_cb"
cr[cr %in% c("MPP_D1", "MPP_D2")] = "MPP_cb"
rnact <- colnames(d) <- cr
wct <- colnames(w) <- cw

##
int <- intersect(rnact,wct)
d.bak <- d <- d[,colnames(d) %in% int]
w.bak <- w <- w[,colnames(w) %in% int]

aggregatefunc <- function(d){
  ## aggregate the data frame according to colnames. 
  ## same colnames will only preserve average
  mat = matrix(0, nrow = ncol(d), ncol = length(unique(colnames(d))))
  for (id in 1:length(unique(colnames(d)))){
    mat[which(colnames(d)==unique(colnames(d))[id]), id] = 1
  }
  d.combine = d %*% mat
  colnames(d.combine) = unique(colnames(d))
  return(d.combine)
}

aggregatefunc2 <- function(d, fun = 'mean'){
  ## using na.rm = T, ignore the NA case when averaging
  d2 = sapply(unique(colnames(d)), function(i){
    idv = which(colnames(d) == i)
    if (length(idv) > 1) {
      d.tmp = d[, colnames(d) %in% i, drop = FALSE]
      if (fun == 'mean') {
        return(rowMeans(d.tmp, na.rm = T))
      } else if (fun == 'sum'){
        return(rowSums(d.tmp, na.rm = T))
      }
    } else {
      return(d[, colnames(d) %in% i])
    }
  })
}

d2 = aggregatefunc2(d, fun = 'sum')  
w2 = aggregatefunc2(w, fun = 'mean')  
d = d2
w = w2

gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='transcript',]
trid <- sub('\\..*','',sub('.*transcript_id "','',gtf[,9]))
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
d <- d[rownames(d) %in% trid,]
d <- rowsum(d,gn[match(rownames(d),trid)])

gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grch38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='exon',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
gl <- tapply(gtf[,5]-gtf[,4],list(gn),sum)
glv <- as.vector(gl)
names(glv) <- names(gl)
d <- d/glv[rownames(d)]*1e3
d <- t(t(d)/colSums(d)*1e6)
d <- log2(d + 1)
saveRDS(d,file='/home/whou10/data/whou/metpred/data/hsc/proc/rna_matchedtrain.rds')
saveRDS(w,file='/home/whou10/data/whou/metpred/data/hsc/proc/wgbs_matchedtrain.rds')

