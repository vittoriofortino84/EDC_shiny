
#FUNCTION jaccard_dissimilarity returns  n*n matrix of pairwise jaccard dissimilarity from n*m binary matrix 
#         the columns must have names
jaccrd_dissimilarity<-function(pro){ #row wise calculation of jaccard index
  jac<-matrix(NA, nrow = dim(pro)[1], ncol = dim(pro)[1])
  nu<-dim(pro)[1]
  nu0<-nu-1
  for (i in 1:nu0){
    er<-i+1
    for (j in er:nu){
      g1<-colnames(pro)[which(pro[i,] %in% 1)]
      g2<-colnames(pro)[which(pro[j,] %in% 1)]
      tmp<-1-(length(intersect(g1,g2))/length(union(g1,g2)))
      jac[i,j]<-tmp
      jac[j,i]<-tmp
    } #end for j
  }#end for i
  diag(jac)<-0
  rownames(jac)<-colnames(jac)<-rownames(pro)
  return(jac)
}#end function   

#retunrs n*n from a list with gene sets
list_jaccrd_similarity<-function(pro){ #row wise calculation of jaccard index
  jac<-matrix(NA, nrow = length(pro), ncol = length(pro))
  nu<-length(pro)
  nu0<-nu-1
  for (i in 1:nu0){
    er<-i+1
    for (j in er:nu){
      tmp<-length(intersect(pro[[i]],pro[[j]]))/length(union(pro[[i]],pro[[j]]))
      jac[i,j]<-tmp
      jac[j,i]<-tmp
    } #end for j
  }#end for i
  diag(jac)<-1
  rownames(jac)<-colnames(jac)<-names(pro)
  return(jac)
}#end function





name_order<-function(inp){
  nams<-factor(inp[c(grep(x=inp,pattern = 'Drug_Matrix_Rat_invitro'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Single_Dose_1_day'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_3_days'),
                                grep(x=inp,pattern = 'Drug_Matrix_Rat_invivo_Repeated_Dose_5_days'),
                                grep(x=inp,pattern = 'TG_GATEs_Human_invitro'),
                                grep(x=inp,pattern = 'TG_GATEs_Rat_invitro'),
                                grep(x=inp,pattern = 'Low'),grep(x=inp,pattern = 'Middle'),
                                grep(x=inp,pattern = 'High'),
                                grep(x=inp,pattern = '8_days'),grep(x=inp,pattern = '15_days'),
                                grep(x=inp,pattern = '29_days'),
                                grep(x=inp,pattern = 'Consensus'),
                                grep(x=inp,pattern = 'PPI'))])
  return(nams)
} #end function

#calculates harmonic sum data can be a matrix or a data frame, col.scores are the selected columns
#the rows must have names
harmonic_sum<-function(data,col.scores){
  HS <- rep(NA, nrow(data))
  for(i in 1:nrow(data)){
    dtemp<-as.numeric(data[i,col.scores])
    dtemp<-dtemp[!is.na(dtemp)]
    HS[i] <- sum(sort(dtemp,decreasing=T)/(1:length(dtemp))^2)
  } # end for
  names(HS)<-rownames(data)
  return(HS)
} # end function


mean_function<-function(data,col.scores){
  means<-rep(NA,nrow(data))
  for (i in 1:length(means)){
    dtemp<-as.numeric(data[i,col.scores]) 
    means[i]<-  mean(dtemp,na.rm = T)
  }
  names(means)<-rownames(data)
  return(means)
}
