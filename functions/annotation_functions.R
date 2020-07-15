

mesh2name<-function(inp,ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(ixns$mesh==inp[i])])>0){
    nn[i]=ixns$name[which(ixns$mesh==inp[i])]
    }
    else {
      nn[i]=''
    }
  }
return(nn)
} #end function

name2mesh<-function(inp,ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(tolower(ixns$name)==tolower(inp[i]))])>0){
      nn[i]=ixns$mesh[which(tolower(ixns$name)==tolower(inp[i]))]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

name2cas<-function(inp,ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(tolower(ixns$name)==tolower(inp[i]))])>0){
      nn[i]=ixns$cas[which(tolower(ixns$name)==tolower(inp[i]))]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

cas2name<-function(inp,
                   ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F),
                   return_back_NA='TRUE'){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$name[which(ixns$cas==inp[i])])>0){
    nn[i]=ixns$name[which(ixns$cas==inp[i])]
    }
    else {
      if(return_back_NA=='TRUE')nn[i]=''
    }
  }
  return(nn)
} #end function

cas2mesh<-function(inp,ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F),
                   return_back_NA='TRUE'){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$mesh[which(ixns$cas==inp[i])])>0){
    nn[i]=ixns$mesh[which(ixns$cas==inp[i])]
    }
    else {
      if(return_back_NA=='TRUE')nn[i]=''
    }
  }
  return(nn)
} #end function

mesh2cas<-function(inp,
                   ixns=read.csv('inputData/annotaion/chem_mesh_cas_dictionary.csv',stringsAsFactors = F)){
  nn<-inp
  for (i in 1:length(inp)){
    if (length(ixns$cas[which(ixns$mesh==inp[i])])>0){
    nn[i]=ixns$cas[which(ixns$mesh==inp[i])]
    }
    else {
      nn[i]=''
    }
  }
  return(nn)
} #end function

pathway_annotate<-function(inp,patway_dic=readRDS('inputData/annotaion/pathway_annot_dictionary.rds')){
  patway_dic$pathway<-as.character(patway_dic$pathway)
  patway_dic$pathway_annotation<-as.character(patway_dic$pathway_annotation)
  nn<-inp
  
  for (i in 1:length(inp)){
    if (any(inp[i] %in% patway_dic$pathway)){
    nn[i]<-patway_dic$pathway_annotation[which(patway_dic$pathway %in% inp[i])]
    } #end if
  }
  return(nn)
} #end function


