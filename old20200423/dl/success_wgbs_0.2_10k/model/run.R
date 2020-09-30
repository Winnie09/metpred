library(keras) ## platform for DL
trainx <- readRDS('/home/ec2-user/run/data/trainx.rds') ## traiconvning gene expr
trainy <- readRDS('/home/ec2-user/run/data/trainy.rds') ##  training DNAm
trainym <- rowMeans(trainy)  ## for each CpG, DNAm mean
trainysd <- apply(trainy,1,sd) ##  for each CpG, DNAm sd
xchr <- readRDS('/home/ec2-user/run/data/seqmat.rds') ## DNA sequence info, used for training and prediction, 4 * 500
nsite <- nrow(trainy) 
nexp <- nrow(trainx)
nchr <- ncol(xchr)/4
xchr <- array_reshape(xchr,dim = c(nsite,4,nchr,1),order = 'C')

dnaseq_input <- layer_input(shape = c(4,nchr,1), dtype = 'float32', name = 'main_input')

dnaseq_out <- layer_conv_2d(dnaseq_input,filters = 100,strides = c(1,1), kernel_size = c(5,5), padding = 'same',activation = 'relu') %>% ##  CNN 
  ## first layer CNN
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size = c(1,4)) %>%
  layer_dropout(0.2) %>%
  layer_conv_2d(filters = 100, kernel_size = c(5,5), padding = 'same',strides = c(1,1),activation = 'relu') %>%
  ## second layer CNN
  layer_batch_normalization() %>%
  layer_max_pooling_2d(pool_size = c(1,4)) %>%
  layer_dropout(0.2) %>%
  # layer_conv_2d(filters = 150, kernel_size = c(1,12), padding = 'same',strides = c(1,1),activation = 'relu') %>%
  # layer_batch_normalization() %>%
  # layer_max_pooling_2d(pool_size = c(1,4)) %>%
  # layer_dropout(0.2) %>%
  layer_flatten()  ## flatten matrix to a vector 

expr_input <- layer_input(shape = c(nexp), name = 'aux_input') ## prior mean
meansd_input <- layer_input(shape = 2) ## prior sd

main_output <- layer_concatenate(c(dnaseq_out, expr_input,meansd_input)) %>%  ## 3 hidden layers, 1 sigmoid layer
  layer_dense(units = 32, activation = 'relu') %>% 
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dense(units = 32, activation = 'relu') %>% 
  layer_dense(units = 1,activation = 'sigmoid') ## use sigmoid becasue beta values in [0,1]

model <- keras_model(
  inputs = c(dnaseq_input, expr_input, meansd_input),
  outputs = c(main_output)
)

#parallel_model <- multi_gpu_model(model, gpus=2)



model %>% compile(
  loss = "binary_crossentropy", ## loss function, minimize, use gradient decent to optimize
  #optimizer = 'adam'
  optimizer = optimizer_adam(lr = 0.01)
)

testx <- readRDS('/home/ec2-user/run/data/testx.rds')
testy <- readRDS('/home/ec2-user/run/data/testy.rds')


## divide the matrix into 20 blocks and train each one, otherwise the full matrix is too large to train
set.seed(12345)
for (looptime in c(1:20)) {
  # sampcellid <- sample(1:ncol(trainy),50000,replace = T)
  # sampsiteid <- sample(1:nrow(trainy),50000,replace = T)
  # samptrainy <- trainy[cbind(sampsiteid,sampcellid)]
  # samptrainxchr <- xchr[sampsiteid,,,,drop=F]
  # samptrainxexp <- t(trainx[,sampcellid])
  # sampmeansd <- cbind(trainym,trainysd)[sampsiteid,]
  gridid <- expand.grid(((looptime-1)*500+1):(looptime*500),1:ncol(trainx))
  sampcellid <- gridid[,2]
  sampsiteid <- gridid[,1]
  samptrainy <- trainy[cbind(sampsiteid,sampcellid)]
  samptrainxchr <- xchr[sampsiteid,,,,drop=F]
  samptrainxexp <- t(trainx[,sampcellid])
  sampmeansd <- cbind(trainym,trainysd)[sampsiteid,]
  ## training
  model %>% fit(
    x = list(samptrainxchr,samptrainxexp,sampmeansd),
    y = list(samptrainy),
    epochs = 50,
    batch_size = 200
  )
  
  xmeansd <- cbind(trainym,trainysd)
  ## prediction
  pred <- sapply(colnames(testx),function(cellid) {
    print(cellid)
    xexp <- t(matrix(rep(testx[,cellid],nsite),nrow=nexp))
    tmp <- predict(model,list(xchr,xexp,xmeansd))
  })
  ## calculate cor using test sets, equivalent to cross validation
  samplecor <- sapply(colnames(testy),function(i) {
    cor(testy[,i],pred[,i])
  })
  
  sitecor <- sapply(1:nrow(testy),function(i) {
    cor(testy[i,],pred[i,])
  })
  print(summary(samplecor))
  print(summary(sitecor))
}

meansamplecor <- apply(testy,2,cor,trainym)
saveRDS(list(samplecor=samplecor,sitecor=sitecor,meansamplecor=meansamplecor),file='cor.rds')
