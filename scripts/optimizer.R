
# 1.prep ------------------------------------------------------------------


load('inputData/network/wTO_DM.RData')
x<-wTO_DM$DM_hepatocytes.1
x$Node.1<-as.character(x$Node.1) #pairwise calculation so length(node1)=length(node2)=length(wTO)=length(pval)=length(pval.adjusted)      
x$Node.2<-as.character(x$Node.2)
x$wTO[x$pval.adj>5*10^-3]<-0 
x$wTO<-abs(x$wTO)
DM_hep_data_frame<-as.data.frame(list(x$Node.1,x$Node.2,x$wTO))
colnames(DM_hep_data_frame)<-c('Node.1','Node.2','wTO')
save(DM_hep_data_frame,file = 'inputData/network/DM_hepatocytes.RData')
###

# Optimization ------------------------------------------------------------

rm(list = ls())
source('functions/network_optimizer.R')
load('inputData/network/DM_hepatocytes.RData')

load('outputData/edc_decoys.RData')
files<-ls()
network<-files[which(sapply(files, function(x)class(get(x))) %in% 'data.frame')]
load('outputData/edc_decoys.RData')

perclist<-c(0.02,0.03,0.05,.06,.07,0.1)        # which percent of the edges to extract for randomwalking
geneset<-c(200)            # The number of most top visited genes to consider for jaccard dissimilarity

n_cpu<-6
edgelist<-cpu_splitter(n_cpu,perclist)


library(doParallel)
cl<-makeCluster(length(edgelist))
registerDoParallel(cl)
modls<-foreach(i=1:length(perclist),.packages=c('dnet','igraph','cluster')) %dopar% network2silhouette(DM_hep_data_frame,
                                                                                                       edgelist[[i]],
                                                                                                       geneset,
                                                                                                       edcs_list,decoys_list)
stopCluster(cl)
names(modls)<-perclist
silhouette_mat<-as.data.frame(do.call(rbind,modls))
silhouette_mat$edge<-rownames(silhouette_mat)
paret_input<-tidyr::gather(silhouette_mat,key='genes',value='silhouette',-edge)
paret_input$genes<-as.numeric(paret_input$genes)
paret_input$edge<-as.numeric(paret_input$edge)
ps<-rPref::psel(paret_input,high(silhouette)*low(genes)*low(edge))

n_genes<-ps$genes[which(ps$silhouette==max(ps$silhouette))]
n_edges<-ps$edge[which(ps$silhouette==max(ps$silhouette))]




