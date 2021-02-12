edc_score<-function(all_vam_library,inp,all_data_levels,col.scores=1:24){
  require(tidyr)
  
  one_compound_harmonic_sum<-function(data,col.scores){
    dtemp<-as.numeric(data[col.scores])
    dtemp<-dtemp[!is.na(dtemp)]
    HS <- sum(sort(dtemp,decreasing=T)/(1:length(dtemp))^2)
    return(HS)
  } # end function
  
  
  one_compound_mean_function<-function(data,col.scores){
    dtemp<-as.numeric(data[col.scores]) 
    means<-  mean(dtemp,na.rm = T)
    return(means)
  }
  
  
  res<-matrix(NA, nrow = length(inp), ncol = ncol(all_vam_library),dimnames = list(inp,colnames(all_vam_library)))
  for (i in 1:length(inp)){
    if (any(tolower(inp[i]) %in% tolower(all_vam_library$comp_names)))res[i,]<-as.matrix(all_vam_library[tolower(all_vam_library$comp_names)==tolower(inp[i]),])
    if (any(inp[i] %in% all_vam_library$mesh)) res[i,]<-as.matrix(all_vam_library[all_vam_library$mesh==inp[i],])
    if (any(inp[i] %in% all_vam_library$cas)) res[i,]<-as.matrix(all_vam_library[all_vam_library$cas==inp[i],])
    res[i,"average_edc_score"]<-one_compound_mean_function(res[i,1:24],col.scores)
    res[i,"harmonic_sum_edc_score"]<-one_compound_harmonic_sum(res[i,1:24],col.scores)
  }
  # 
  res<-as.data.frame(res)
  res$nnames<-paste(res$comp_names,res$mesh,res$is_in_training,sep = '_')
  # colnames(res)<-gsub(x=colnames(res),
  #                     pattern ="Consensus_Rat_invitro_Drug.Matrix_TG_GATEs",
  #                     replacement = "Consensus_Rat_invitro_DM_TG_GATEs" )
  res<-res[,c(all_data_levels,'nnames')]

  #
 
  res<-tidyr::gather(res,key = 'network',value = 'score',-c(nnames))
  res$score<-as.numeric(res$score)
  
  return(res)
} # end function
