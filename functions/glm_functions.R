

#making model from glm
elastic_net_model_two_class<-function(data,lab,seed=1){ 
  set.seed(seed)
  require(caret)
  require(glmnet)
  trcntrl <- trainControl(method = "repeatedcv",
                                number = 5,
                                repeats =2,
                                search = "random",
                                verboseIter = FALSE,
                                classProbs=TRUE,  
                                summaryFunction=twoClassSummary)
  
  
  modl<-train(data, lab,
              method = "glmnet",
              family = "binomial",
              metric = "ROC",
              tuneLength = 20,
              trControl = trcntrl)
  
  mfit.ctcf.2.coef = as.matrix(coef(modl$finalModel, modl$bestTune$lambda))
  sel<-rownames(mfit.ctcf.2.coef)[which(abs(mfit.ctcf.2.coef)>0.0)] 
  return(list(modl=modl,sel=sel,cof=mfit.ctcf.2.coef))

}   #two_class

splitted_CV<-function(cpu_index,caret_index,x,y){
  require(caret);require(glmnet)
  all_F1<-rep(NA,length(cpu_index))
  for (i in 1:length(cpu_index)){
    to_pick<-cpu_index[i]
    train_index<-caret_index[[to_pick]];x_train<-x[train_index,];y_train<-y[train_index]
    test_index<-which(!1:length(y) %in% train_index);x_test<-x[test_index,];y_test<-y[test_index]
    
    modl<-elastic_net_model_two_class(x_train,y_train)
    predicted_labels<-predict(modl$modl,x_test)
    conf_mat<-confusion_matrix(y_test,predicted_labels)
    acuracy<-mean(cm_based_evaluation(conf_mat)$f1)
    all_F1[i]<-acuracy
  }
  return(all_F1)
}



confusion_matrix<-function(real,pred){#two vectors of real observations and predicted observations it can be two or multiclass labeled
  all_conditions<-unique(union(real,pred))
  cm<-matrix(NA,ncol = length(all_conditions),nrow = length(all_conditions))
  rownames(cm)<-colnames(cm)<-all_conditions
  for (i in all_conditions){
    for (j in all_conditions){
      cm[i,j]<-length(which(real %in% i ==pred %in% j & real %in% i==TRUE & pred %in% j==TRUE))
    } #end for pred
  }   # end for real
  return(cm)
}# end function

cm_based_evaluation <- function(cm, type = "basic"){ ## http://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class
  p = rowsums / n # distribution of instances over the actual classes
  q = colsums / n # distribution of instances over the predicted classes
  # metrics
  accuracy = sum(diag) / n 
  precision = diag / colsums 
  recall = diag / rowsums 
  # macro-averaged Metrics
  f1 = 2 * precision * recall / (precision + recall) 
  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1,na.rm = T)
  # one vs all
  oneVsAll = lapply(1 : nc,function(i){
    v = c(cm[i,i],
          rowsums[i] - cm[i,i],
          colsums[i] - cm[i,i],
          n-rowsums[i] - colsums[i] + cm[i,i]);
    return(matrix(v, nrow = 2, byrow = T))})
  s = matrix(0, nrow = 2, ncol = 2)
  for(i in 1 : nc){s = s + oneVsAll[[i]]}
  #print(s)
  # - avearge accuracy
  avgAccuracy = sum(diag(s)) / sum(s)
  # - micro average metrics
  micro_prf = (diag(s) / apply(s,1, sum))[1]
  # Majority-class Metrics
  mcIndex = which(rowsums==max(rowsums))[1] # majority-class index
  mcAccuracy = as.numeric(p[mcIndex]) 
  mcRecall = 0*p;  mcRecall[mcIndex] = 1
  mcPrecision = 0*p; mcPrecision[mcIndex] = p[mcIndex]
  mcF1 = 0*p; mcF1[mcIndex] = 2 * mcPrecision[mcIndex] / (mcPrecision[mcIndex] + 1)
  if(type == "basic")
    return(list(accuracy=accuracy, precision=precision, recall=recall, f1=f1))
  else if(type == "macro")
    return(list(totAcc = accuracy, allPrec = precision, allRec = recall, allF1 = f1,
                macAcc = avgAccuracy, macPrec = macroPrecision, macRec = macroRecall, macF1 = macroF1,
                micAcc = micro_prf))
  else if(type == "majority")
    return(list(mcAccuracy=mcAccuracy, mcPrecision=mcPrecision, mcRecall=mcRecall, mcF1=mcF1))
}


