library(superpc)
library(randomForestSRC)
library(mboost)
source("Cindex.R")
library("e1071")
library("pROC")
library(Matrix)
library(Coxnet)
library(survcomp)
library(fastcox)
library(glmnet)
library(survival)
source('E:/RWorkspace/TCGA/data_585_mklsvm/min_max_nor.R')
########################################
geneExp_1020 <- read.table('geneExp_1020.txt')
CNA_1020 <- read.table('CNA_1020.txt')
Methyl_1020 <- read.table('Methyl_1020.txt')
CP_1020 <- read.table('CP_1020.txt')
protein_1020 <- read.table('protein_1020.txt')
label <- as.matrix(read.table('label_585.txt'))
time_os <- as.matrix(read.table('os_585.txt'))
status <- as.matrix(read.table('status_585.txt'))

time_os_nor <- min_max_nor(time_os)
####

time_os[which(time_os<=0)] <- 0.00001#(0) ALL
xx <-  as.matrix(cbind(geneExp_1020,CNA_1020,Methyl_1020,protein_1020,CP_1020))

#(1) exp
xx <-  as.matrix(geneExp_1020)
#(2) CP
xx <-  as.matrix(CP_1020)
#(3) linear_CNA
xx <-  as.matrix(CNA_1020)
#(4) methy
xx <-  as.matrix(Methyl_1020)
#(5) protein
xx <- as.matrix(protein_1020)
#(6)-CP
xx <-  as.matrix(cbind(geneExp_1020,CNA_1020,Methyl_1020,protein_1020))
#(7) 2
xx <-  as.matrix(cbind(Methyl_1020,CNA_1020))
#(7) 3
xx <-  as.matrix(cbind(geneExp_1020,CNA_1020,Methyl_1020))
#(7) 4
xx <-  as.matrix(cbind(geneExp_1020,CNA_1020,Methyl_1020,protein_1020))

i <- 0
#numepochs_20 <- c(200,100,50,150,150,50,100,200,150,100,100,150,50,150,150,50,150,100,100,100)
auc_m = matrix(0,20,1)
for (epos in (1:20)){
  i <- i + 1
  set.seed(epos)
  foldid=sample(rep(seq(585),length=585))
  rate = 0.8
  train_ind = foldid[1:floor(rate*585)]
  test_ind = foldid[(floor(rate*585)+1):585]
  
  train_data <- xx[train_ind,]
  test_data <- xx[test_ind,]
  train_label <- label[train_ind]
  test_label <- label[test_ind]
  train_time_os <- time_os[train_ind]
  test_time_os <- time_os[test_ind]
  test_status <- status[test_ind]
  

  
  #(1) svm
  #model <- svm(train_data, train_label)
  #pred_train <- predict(model, train_data)
  #test
  #pred_test <- predict(model, test_data)
  #modelroc1 <- roc(test_label,pred_test)
  #print(modelroc1$auc)
 # cindex <- concordance.index(1-pred_test,surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  
  #(2) superpc
  
  #data_train<-list(x=t(train_data),y=train_label)
  #data_test<-list(x=t(test_data),y=test_label)
  
  #superpc_m<- superpc.train(data_train, type="regression")
  #fit<- superpc.predict(superpc_m,data_train, data_test, threshold=1.0, n.components=0,prediction.type="continuous")
  
  #cindex <- concordance.index(1-as.numeric(fit$v.pred),surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  
  
  #modelroc1 <- roc(test_label,as.numeric(fit$v.pred))
  #print(modelroc1$auc)
  
  
  #(3)
  #fit <- glmnet(train_data, Surv(as.double(time_os[train_ind]),as.double(status[train_ind])), family = "cox",alpha = 0.5)
  #pred_glm = predict(fit, test_data,type="response",s=0.01)
  #modelroc1 <- roc(test_label,pred_glm)
  #print(modelroc1$auc)
  
  
  #cindex <- concordance.index(pred_glm,surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  
  #(4)
  #mytime <- as.double(time_os[train_ind])
  #mystatus <- as.double(status[train_ind])
  #my_data <- data.frame(mytime,mystatus,train_data)

  #fit <- survreg(Surv(mytime, mystatus) ~., my_data, dist="exponential")
  #pred <- predict(fit, newdata=data.frame(test_data), type='response')
  #modelroc1 <- roc(test_label,pred)
  #print(modelroc1$auc)
  
  
  #cindex <- concordance.index(1-pred,surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  
  #(5) rfsrc
  #mytime <- as.double(time_os_nor[train_ind])
  #mystatus <- as.double(status[train_ind])
  #my_data <- data.frame(mytime,mystatus,train_data)
  
  #fit_rfsrc <- rfsrc(Surv(mytime, mystatus) ~., my_data, ntree = 100)
  #pred_rfsrc <- predict(fit_rfsrc, newdata=data.frame(test_data))
  
  #modelroc1 <- roc(test_label,pred_rfsrc$predicted)
  #print(modelroc1$auc)
  
  #cindex <- concordance.index(pred_rfsrc$predicted,surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  #(6)  boostCI
  mytime <- as.double(time_os_nor[train_ind])
  mystatus <- as.double(status[train_ind])
  my_data <- data.frame(mytime,mystatus,train_data)
  mod_boostCI <- glmboost(Surv(mytime, mystatus) ~ ., family = Cindex(sigma = 0.1),
                   control = boost_control(mstop = 500, trace = TRUE, nu = 0.1),
                   data = my_data)
  mod_boostCI <- mod_boostCI[5000]
  pred_boostCI <- predict(mod_boostCI, newdata = data.frame(test_data))
  #modelroc1 <- roc(test_label,as.numeric(pred_boostCI))
  #auc_m[i] = modelroc1$auc
  #print(modelroc1$auc)
  
  cindex <- concordance.index(1-pred_boostCI,surv.time = test_time_os, surv.event = test_status,method = "noether")
  #print(cindex$c.index)
  auc_m[i] =cindex$c.index
}


