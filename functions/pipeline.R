
# ptw is pathways
# perc is the percentage of the edges from the network to be used and 
# ng is the number of sorted genes after random walk to be used for fgsea
# com_seed is the seeds or MIEs for the compounds
# p is parralellization  option
pipeline<-function(network,ptw,perc,ng,comp_seeds,p){
  options(warn=-1)
  
  require(fgsea)
  require(dnet)
  require(igraph)
  
  #graph optimization

  gr<-graph_from_data_frame(network[order(network[,3],decreasing = T),1:2][1:round(nrow(network)*perc),],directed=F)  #getting the   percent of weighted edge values
  
  
    #gettting the   percent of weighted edge values
  d_sign<-lapply(comp_seeds,function(x)intersect(x,V(gr)$name))    #intersection of the genes in the graph and the genes related to the compounds
  d_sign<-d_sign[sapply(d_sign,length)>0]                     #omitting chemicals with no intersected genes inside the graph 
  mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)           #matrix of the genes in the row and chemicals in the columns
  rownames(mat)<-V(gr)$name 
  
  # Randomwalk With Restart (RWR)
  probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
  probm<-as.matrix(probm)
  colnames(probm)<-names(d_sign)                        #Columns are chemicals
  rownames(probm)<-V(gr)$name
  for (i in 1:dim(probm)[2]){
    cutof_q<-sort(probm[,i],decreasing = T)[ng]
    probm[probm[,i]<cutof_q,i]<-0
  }
  
  
  # FGSEA
  q<-apply(probm,2,function(y)sum(y!=0)) # Counts the number of genes for each chemical which are at least once seen in random walk (probability >0)
  probm<-probm[,q>=ng]                     # We take the chemicals with at least more than ng observed genes after random walk on the GCN
  probm<-as.matrix(probm)
  colnames(probm)<-names(q)[which(q>=ng)]
  ES<-NES<-pval<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),dimnames=list(names(ptw),colnames(probm))) #making a matrix with null values  the rows will be the pathways and the columns will be the compounds
 
   for(i in 1:ncol(probm)){
    #for each chemical and its sets of genes we do this procedure
    res<-fgsea(pathways=ptw,stats=probm[probm[,i]!=0,i],nperm=10000,minSize=1,maxSize=200,scoreType='pos',BPPARAM=p)
    #the pathways more than 200 genes are not considered and those with at least one gene will be considered
    
    ES[res$pathway,i]<-res$ES
    NES[res$pathway,i]<-res$NES
    pval[res$pathway,i]<-res$pval

    if(i%%100==0) gc()           # at the end of the loop clear the memory
  }
  
  return(list(ES=ES,NES=NES,pval=pval))
}# end function


# based on the precompiled networks

onecompound_mies2classprob<-function(networks,patways,models,compound_mies){
  options(warn=-1)
  
  require(fgsea)
  require(dnet)
  require(igraph)
  grp<-1:length(networks$networks)
  cl_probs<-rep(NA,length(grp))
  withProgress(message = "Wait for RWR-FGSEA-GLM",value = 0,{
  for (i in 1:length(grp)) {

    incProgress(1/length(grp),detail = paste('network',i,' out of ',length(grp)))
    indd<-grp[i]
    gr<-networks$networks[[indd]]
    ng<-networks$genes_number[indd]
    ptw<-patways$all_patways_list[names(patways$all_patways_list) %in% patways$used_patways_in_models[[indd]]]
    
    d_sign<-lapply(compound_mies,function(x)intersect(x,V(gr)$name))  #intersection of the genes in the graph and the genes related to the compounds
    
    if (length(unlist(d_sign))>0){
      mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)             #matrix of the genes in the row and chemicals in the columns
      rownames(mat)<-V(gr)$name 
      
      #random walk
      probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=F,verbose = F) 
      probm<-as.matrix(probm)
      colnames(probm)<-names(d_sign)                        #Columns are chemicals
      rownames(probm)<-V(gr)$name 
      cutof_q<-sort(probm[,1],decreasing = T)[ng]
      probm[probm<cutof_q]<-0
      
      NES<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),dimnames=list(names(ptw),colnames(probm))) 
      # fgsea
      res<-fgsea(pathways=ptw,stats=probm[probm[,1]!=0,1],nperm=10000,minSize=1,maxSize=200,scoreType='pos')
      NES[res$pathway,1]<-res$NES
      NES[is.na(NES)]<-0
      NES<-t(NES)
      cl_probs[i]<-predict(models[[indd]]$modl,NES,type='prob')$edc
    }
  }
  })
  names(cl_probs)<-names(models)[grp]
  return(cl_probs)
}
# Example
# all_pipeline_precompiled_data<-readRDS('inputData/all_precompiled_pipeline.rds')
# selected_layers<-c("TG_GATEs_Rat_invivo_Repeated_Dose_8_days" ,"TG_GATEs_Rat_invivo_Single_Dose_High_1_day" )
# patways<-all_pipeline_precompiled_data$pathways
# patways$used_patways_in_models<-patways$used_patways_in_models[selected_layers]
# networks<-lapply(all_pipeline_precompiled_data$networks,function(x)x[selected_layers])
# models<-all_pipeline_precompiled_data$models[selected_layers]
#newcomp<-list(c("4617","4654" ,"4656" ,"5077","25937")) #genes as MIEs 
# example3<-onecompound_mies2classprob(networks,patways,models,newcomp)


