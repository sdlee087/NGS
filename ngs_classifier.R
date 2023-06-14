library(MASS)
library(tidyverse)
library(glmnet)
library(class)
library(rpart)
library(sparsepca, rpca)


logistic.kcv=function(trainX, trainY, testX, testY){
  k=5
  n = nrow(trainX)
  p = ncol(trainX)
  grid=10^seq(-4, 0, length=100)
  cv.err.vec = rep(0, length(grid))
  folds = rep_len(1:k, n)
  
  for(j in 1:k){
    tr=which(folds!=j)
    val=which(folds==j)
    x_train=trainX[tr, ]
    y_train=trainY[tr]
    x_val = trainX[val,]
    y_val = trainY[val]
    
    fit_train=glmnet(x_train, y_train, lambda=grid, family = "binomial")
    
    for (i in 1:length(grid)){
      pred=predict(fit_train, newx=x_val, s=grid[i], type = "class")
      cv.err.vec[i]=cv.err.vec[i]+sum((as.numeric(pred)-y_val)^2)
    }
  }
  best.index=max(which(cv.err.vec==min(cv.err.vec)))
  best.lam=grid[best.index]
  best.cv.err=cv.err.vec[best.index]/n
  
  best.model=glmnet(trainX, trainY, alpha=1, lambda=best.lam, family = "binomial")
  pred_test=predict(best.model, newx=testX, s=best.lam, type = "class")
  test.err=sum((as.numeric(pred_test)-testY)^2)/nrow(testX)
  
  return(list("best.model" = best.model, "best.lam" = best.lam, "cv.acc" = 1 - best.cv.err, "test.acc" = 1 - test.err))
}

pred_err = function(trainX, trainY, testX, testY, model, par_tuning, extended="0"){
  trainY = as.matrix(trainY)
  testY = as.matrix(testY)
  if(model=="logistic"){
    trainX = as.matrix(trainX)
    testX = as.matrix(testX)
    fit_train = glmnet(trainX, trainY, lambda=par_tuning[1], family = "binomial")
    pred = as.numeric(predict(fit_train, newx=testX, s=par_tuning[1], type = "class")) # par_tuning[1] or par_tuning[p+1] ?
  } else if(model=="knn"){
    trainX = as.matrix(trainX)
    testX = as.matrix(testX)
    pred = as.numeric(as.vector(knn(trainX, testX, trainY, k=par_tuning[1])))
  } else if(model=="tree"){
    train.df = as.data.frame(cbind(as.factor(trainY), trainX))
    train.df = setNames(train.df, c("Y", paste0("Var", 1:ncol(trainX))))
    testX = setNames(as.data.frame(testX), paste0("Var", 1:ncol(trainX)))
    
    fit_train = rpart(Y~., data=train.df,
                      control=rpart.control(minsplit = par_tuning[1], # 1,3,5,7,9
                                            minbucket = par_tuning[2], # 1,3,5,7,9
                                            xval = 0, 
                                            maxdepth = par_tuning[3])) # 1,3,5,7,9
    pred <- as.numeric(as.vector(predict(fit_train, as.data.frame(testX), type="class")))
    test.err = mean((pred-testY)^2);test.err
  }
  
  test.err = mean((pred-testY)^2)
  
  if(extended=="0"){
    return(list("acc" = 1-test.err))
  } else if(extended=="1"){
    return(list("acc" = 1-test.err, "prediction" = pred))
  } else if(extended=="2"){
    TP = sum(pred*testY)
    precision = TP/sum(pred) # Here, NaN means 0/0
    recall = TP/sum(testY)
    return(list("acc" = 1-test.err, "prediction" = pred, "precision" = precision, "recall" = recall,
                "F1.score" = 2*precision*recall/(precision+recall)))
  }
}


# load and select data
genes = c("abl1", "bcr", "mll", "runx1", "etv6")
data = list()
for(i in 1:5){
  data[[i]] = read_csv(paste0("data/", genes[i], ".csv"))
}
panel <- read_csv("data/panel_info.csv")
panel[ , 4:8] = abs(panel[ , 4:8]-1) # 0-1 exchanged


# CV-test split

val_X1=val_X2=val_X3=val_X4=test_X=list(0) # 4-fold used
val_y1=val_y2=val_y3=val_y4=test_y=y=test_ind=list(0)

set.seed(123)
for(i in 1:5){
  n_row = nrow(data[[i]])
  y = panel[ , genes[i]][1:n_row, ]
  ind_0 = which(y==0); ind_1 = which(y==1)
  ind_0 = sample(ind_0, length(ind_0), replace = F) # shuffle
  ind_1 = sample(ind_1, length(ind_1), replace = F)
  split_0 = seq(1, length(ind_0), length.out=6)
  split_1 = seq(1, length(ind_1), length.out=6)
  
  test_ind[[i]] = c(ind_0[split_0[1]:split_0[2]],
                    ind_1[split_1[1]:split_1[2]])
  test_X[[i]] = data[[i]][test_ind[[i]], -1] # remove an index column
  test_y[[i]] = y[test_ind[[i]], ]
  
  for(j in 1:4){
    .GlobalEnv[[paste0("val_X",j)]][[i]] = data[[i]][c(ind_0[split_0[j+1]:split_0[j+2]],
                                                       ind_1[split_1[j+1]:split_1[j+2]]), -1]
    .GlobalEnv[[paste0("val_y",j)]][[i]] = y[c(ind_0[split_0[j+1]:split_0[j+2]],
                                               ind_1[split_1[j+1]:split_1[j+2]]), ]
  }
}


# Sample graph
for_tib = c()
for(i in 1:10){
  for_tib = bind_rows(for_tib, tibble(pos = 1:p, cov = unlist(train_X[,-1] %>% slice(i),use.names=F), index=i, normal = y[i]))
}
for_tib = for_tib %>% mutate(index = as.factor(index), normal = as.factor(normal))
for_tib %>% mutate(cov = sqrt(cov), index %in% 1:5) %>% ggplot(aes(x = pos, y = cov, colour = index, linetype = normal)) + geom_line() + labs(title = 'sqrt')

# preprocess
train_mat = sqrt(as.matrix(train_X[,-1]))
test_mat = sqrt(as.matrix(test_X[,-1]))

# prediction example

i = 5

pred_err(trainX=rbind(val_X1[[i]],val_X2[[i]],val_X3[[i]],val_X4[[i]]), 
         trainY=rbind(val_y1[[i]],val_y2[[i]],val_y3[[i]],val_y4[[i]]), 
         testX=test_X[[i]], testY=test_y[[i]], model="logistic", par_tuning=0.01, extended="2")

pred_err(trainX=rbind(val_X1[[i]],val_X2[[i]],val_X3[[i]],val_X4[[i]]), 
         trainY=rbind(val_y1[[i]],val_y2[[i]],val_y3[[i]],val_y4[[i]]), 
         testX=test_X[[i]], testY=test_y[[i]], model="knn", par_tuning=13, extended="2")

pred_err(trainX=rbind(val_X1[[i]],val_X2[[i]],val_X3[[i]],val_X4[[i]]), 
         trainY=rbind(val_y1[[i]],val_y2[[i]],val_y3[[i]],val_y4[[i]]), 
         testX=test_X[[i]], testY=test_y[[i]], model="tree", par_tuning=c(1,1,1), extended="2")
# {Error: protect(): protection stack overflow} for bcr, etv6


# Dimension Reduction

# Robust PCA
?rpca::rpca
rb = rpca::rpca(train_mat)

# Robust Sparse PCA
?robspca
rbs = robspca(train_mat, k = 10)

# Sparse PCA

?sparsepca::spca
sp <- sparsepca::spca(train_mat, k = 10, verbose = F)

