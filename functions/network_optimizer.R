
network2silhouette<-function(network,perc,nglist,edcs_list,decoys_list){
  require(dnet)
  require(igraph)
  require(cluster)
  
  edc_decoy<-as.array(c(edcs_list,decoys_list))
  
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
  
  silm<-matrix(NA, nrow=length(perc),ncol=length(nglist))
  colnames(silm)<-nglist
  rownames(silm)<-perc
  ind_perc<-0
  for (p in 1:length(perc)){
    ind_perc=ind_perc+1
    gr<-graph_from_data_frame(network[order(network[,3],decreasing = T),1:2][1:round(nrow(network)*perc[p]),],directed=F)  #gettting the   percent of weighted edge values
    d_sign<-lapply(edc_decoy,function(x)intersect(x,V(gr)$name))    #intersection of the genes in the graph and the genes related to the compunds
    
    d_sign<-d_sign[sapply(d_sign,length)>0]                     #omitting chemicals with no intersected genes inside the graph 
    mat<-sapply(d_sign,function(x)(V(gr)$name%in%x)*1)           #matrix of the genes in the row and chemicals in the columns
    rownames(mat)<-V(gr)$name 
    
    #random walk with restart
    probm<-dRWR(gr,setSeeds=mat,normalise.affinity.matrix='quantile',parallel=T,verbose = F) 
    probm<-as.matrix(probm)
    colnames(probm)<-names(d_sign)                        #Columns are chemicals
    rownames(probm)<-V(gr)$name 
    
    ind_ng<-0  
    for (ng in nglist){                       # for each cutoff on number of the top visited genes
      ind_ng<-ind_ng+1
      print(ng)
      pro<-probm
      
      # cutoff based on  number of top visited genes 
      for (i in 1:dim(pro)[2]){
        cutof_q<-sort(pro[,i],decreasing = T)[ng]
        pro[pro[,i]<cutof_q,i]<-0
        pro[pro[,i]>=cutof_q,i]<-1
      }
      
      #jaccard  matrix
      pro<-t(pro)
      jac<-jaccrd_dissimilarity(pro)
      ndecoy<-length(which(names(decoys_list) %in% rownames(pro)))
      nedc<-length(which(names(edcs_list) %in% rownames(pro)))
      resp<-c(rep(1,nedc),rep(0,ndecoy)) 
      km <- kmeans(jac, centers = 2, nstart=25)
      km$cluster<-resp
      ss<-silhouette(km$cluster, jac)
      #print(ss)
      silm[ind_perc,ind_ng]<-mean(ss[1:nedc,3])
      
    } # end for number genesets
  }# end of perc
  return(silm)   
}



cpu_splitter<-function(n_cpu=2,objects){
  plist<-rep(NA,n_cpu)
  plist<-as.list(plist)
  ind<-1
  for (i in 1:length(objects)){
    plist[[ind]]<-c(plist[[ind]],objects[i])
    ind<-ind+1
    if (ind>n_cpu)ind=1
    
  }
  
  flist<-lapply(plist, function(x)na.omit(x))
  flist[which(sapply(flist, length)==0)]<-c()
  return(flist)
}

