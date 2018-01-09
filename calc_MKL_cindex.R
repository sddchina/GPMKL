library(survival)
source('E:/RWorkspace/TCGA/data_585_mklsvm/min_max_nor.R')


time_os <- as.matrix(read.table('os_585.txt'))
status <- as.matrix(read.table('status_585.txt'))
label <- as.matrix(read.table('label_585.txt'))
time_os_nor <- min_max_nor(time_os)
MKL_prediction_score <- read.table('MKL_prediction_score.txt')


for (epos in (1:20)){
  set.seed(epos)
  foldid=sample(rep(seq(585),length=585))
  rate = 0.8
  train_ind = foldid[1:floor(rate*585)]
  test_ind = foldid[(floor(rate*585)+1):585]
  #test_label <- label[test_ind]
  
  train_time_os <- time_os[train_ind]
  test_time_os <- time_os[test_ind]
  test_status <- status[test_ind]
  start = 1+117*(epos-1)
  end = 117+117*(epos-1)
  mkl_label <- MKL_prediction_score[start:end,1]
  mkl_pre <- MKL_prediction_score[start:end,2]
  cindex <- concordance.index(1-mkl_pre,surv.time = test_time_os, surv.event = test_status,method = "noether")
  print(cindex$c.index)
}