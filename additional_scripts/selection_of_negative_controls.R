# Returns negative controls and related set of MIEs for a given set of EDCs


# Pool_compounds: a named list of compounds and theirs MIES
# edcNames      : a character vector with names of EDCs which are also inside the pool_compounds list
# mean_similarity_cutoff: the compounds with average Jaccard similarity bellow this cut-off percent will
# be selected as negative controls

EDC_negativeControl_from<-function(Pool_compounds,edcNames,mean_similarity_cutoff=10){ #row wise calculation of Jaccard index
  jac<-matrix(NA, nrow = length(Pool_compounds), ncol = length(Pool_compounds))
  nu<-length(Pool_compounds)
  nu0<-nu-1
  for (i in 1:nu0){
    er<-i+1
    for (j in er:nu){
      tmp<-length(intersect(Pool_compounds[[i]],Pool_compounds[[j]]))/length(union(Pool_compounds[[i]],Pool_compounds[[j]]))
      jac[i,j]<-tmp
      jac[j,i]<-tmp
    } #end for j
  }#end for i
  diag(jac)<-1
  rownames(jac)<-colnames(jac)<-names(Pool_compounds)
  negativeContorls=Pool_compounds[which(apply(jac[edcNames,], 2,mean)<=mean_similarity_cutoff*0.01)]
  return(negativeContorls)
}


# Example:
# We have a pool list of 131 compounds and their related MIEs, and  10 compounds are known as EDCs in this list
# we get the compounds below average Jaccard similarity of 3% with respect to EDCs as negative controls
load('example_negativeControl.RData')
negative_controls=EDC_negativeControl_from(all,edcNames,mean_similarity_cutoff = 3)