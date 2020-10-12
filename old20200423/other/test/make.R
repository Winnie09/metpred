# ml python/3.7
# ml tensorflow/1.10.1-gpu-py3-keras
# ml cuda
# ml R/3.6.1

library(keras)
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
d <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta.rds')
set.seed(12345)
trainid <- grep('37',colnames(d),value=T,invert=T)
testid <- setdiff(colnames(d),trainid)
sdv <- apply(d[,trainid],1,sd)
id <- which(sdv > 0.25)
#id <- sample(1:nrow(d),10000)
#id <- 1:nrow(d)
d <- d[id,]
ecl <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/cellline/data/data/database/clumean.rds')
d <- d[,colnames(ecl)]
gr <- readRDS('/home-4/zji4@jhu.edu/scratch/encode_compiled/Nov19/mearray/data/filtermat/gr.rds')
gr <- gr[id,]
seq <- toupper(as.data.frame(getSeq(Hsapiens, as.character(seqnames(gr)), start = start(gr)-100, end = end(gr) + 101))[,1])
nchr <- nchar(seq[1])
chrmat <- do.call(rbind,strsplit(seq,''))
chrmat <- matrix(as.numeric(cbind(chrmat=='A',chrmat=='T',chrmat=='G',chrmat=='C')),nrow=nrow(chrmat))
#bm <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta_mean.rds')[id]
#bsd <- readRDS('/home-4/zji4@jhu.edu/scratch/metpred/data/data/mearray/beta_sd.rds')[id]
nsite <- nrow(d)
nexp <- nrow(ecl)
xchr <- array_reshape(chrmat,dim = c(nsite,4,nchr,1),order = 'C')

main_input <- layer_input(shape = c(4,nchr,1), dtype = 'float32', name = 'main_input')

lstm_out <- layer_conv_2d(main_input,filters = 100,strides = c(1,1), kernel_size = c(5,5), padding = 'same',activation = 'relu') %>%
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size = c(1,4)) %>%
  layer_dropout(0.2) %>%
  layer_conv_2d(filters = 100, kernel_size = c(5,5), padding = 'same',strides = c(1,1),activation = 'relu') %>%
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size = c(1,4)) %>%
  layer_dropout(0.2) %>%
  # layer_conv_2d(filters = 150, kernel_size = c(1,12), padding = 'same',strides = c(1,1),activation = 'relu') %>%
  # layer_batch_normalization() %>%
  # layer_max_pooling_2d(pool_size = c(1,4)) %>%
  # layer_dropout(0.2) %>%
  layer_flatten() 

# library(densenet)
# model <- application_densenet(include_top = TRUE, input_tensor = main_input, dropout_rate = 0.2) %>%
#   layer_flatten()

auxiliary_input <- layer_input(shape = c(nexp), name = 'aux_input')

main_output <- layer_concatenate(c(lstm_out, auxiliary_input)) %>%  
  layer_dense(units = 32, activation = 'relu') %>% 
  layer_dense(units = 32, activation = 'relu') %>% 
  #layer_dense(units = 32, activation = 'relu') %>% 
  layer_dense(units = 2, activation = 'sigmoid', name = 'main_output')

model <- keras_model(
  inputs = c(main_input, auxiliary_input), 
  outputs = c(main_output)
)

#parallel_model <- multi_gpu_model(model, gpus=2)

model %>% compile(
  optimizer = 'adam',
  loss = 'categorical_crossentropy',
  metrics = c('accuracy')
)


for (cellid in trainid) {
  print(cellid)
  #xexp <- cbind(t(matrix(rep(ecl[,cellid],nsite),nrow=nexp)),bm,bsd)
  xexp <- t(matrix(rep(ecl[,cellid],nsite),nrow=nexp))
  y <- to_categorical(as.numeric(d[,cellid] > 0.5))
  model %>% fit(
    x = list(xchr,xexp),
    y = list(y),
    epochs = 10,
    batch_size = 64
  )
  save_model_hdf5(model,file=paste0('/home-4/zji4@jhu.edu/scratch/metpred/test/model/',cellid,'.h5'))
}

pred <- sapply(testid,function(cellid) {
  xexp <- t(matrix(rep(ecl[,cellid],nsite),nrow=nexp))
  predict(model,list(xchr,xexp))[,2]
})

true <- d[,testid]
print(summary(sapply(1:ncol(true),function(i) {
  cor(true[,i],pred[,i])
})))

print(summary(sapply(1:nrow(true),function(i) {
  cor(true[i,],pred[i,])
})))


