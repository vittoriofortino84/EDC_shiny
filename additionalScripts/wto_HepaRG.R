#data source:   https://www.ebi.ac.uk/biostudies/studies/S-DIXA-AN-012
# Generation of wTO_HepaRG_1_day Gene co-expression network-----------------------------------
Hepa<-read.csv('inputData/log2ratio/HepaRG/HepaRG_log2_ratios.txt',sep='\t',stringsAsFactors=F,check.names=F)
Hepa_annotation<-read.csv('inputData/log2ratio/HepaRG/HepaRG_Exp_annotations.txt',sep='\t',stringsAsFactors=F)
  rownames(Hepa)<-Hepa$`Row Names`
  Hepa_annotation<-Hepa_annotation[Hepa_annotation$Control.s_ %in% F,] # removing rows containing "control" 
  Hepa_annotation<-Hepa_annotation[Hepa_annotation$Time.point..C. %in% '24h',]  #24h of exposure with the compounds
 
 library(magrittr)
  all_compounds<-Hepa_annotation$Hybridization %>% strsplit('-') %>% lapply('[',3) %>% unlist()

 # mapping of probe to genes
require(xlsx)
 #https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=info:doi/10.1371/journal.pcbi.1004847.s026
pr<-read.xlsx2('inputData/journal.pcbi.1004847.s026.XLS',1,stringsAsFactors=F)   #list of probe ids and their corresponding gene ids   in_9071<-pr$in_9071_rat_liver_expressed_set=='Yes'
  probe_id<-pr$HPH.HGU133.2 
  has_orth<-pr$HUMAN_ORTHOLOG_ENTREZ!=''
  logi<-has_orth&in_9071&(probe_id!='')
  probe_set<-probe_id[logi]
  human_orth<-pr$HUMAN_ORTHOLOG_ENTREZ[logi]
  probe_gene_lib<-as.data.frame(list(human_orth,probe_set),stringsAsFactors = F,col.names=c('gene','probe'))
  Final_Hepa<-Hepa[probe_set,colnames(Hepa) %in% Hepa_annotation[,1]]
 # function used for orthology mapping  
  prob2gene<-function(inp){
    res<-rep(NA,length(inp))
    for (i in 1:length(inp)){
      res[i]<-probe_gene_lib$gene[probe_gene_lib$probe==inp[i]] }
    return(res) }
   rownames(Final_Hepa)<-prob2gene(probe_set)
  Final_Hepa<-as.data.frame(Final_Hepa)
    library(caret)
  # removing near to zero genes over all samples
  near_zero_genes<-nearZeroVar(t(Final_Hepa))
  Final_Hepa<-Final_Hepa[-near_zero_genes,]
    
  # generation of the network 
 library(wTO)
 wTO_HepaRG<-wTO.fast(Data=Final_Hepa,method = 'p', sign = 'sign',method_resampling = "Bootstrap",delta = 0.2) 
 wTO_HepaRG$wTO[wTO_HepaRG$pval.adj>5*10^-3]<-0  
 wTO_HepaRG$wTO<-abs(wTO_HepaRG$wTO) 
 wTO_HepaRG<-as.data.frame(list(wTO_HepaRG$Node.1,wTO_HepaRG$Node.2,wTO_HepaRG$wTO))
 colnames(wTO_HepaRG)<-c('gene1','gene2','wTO')
 saveRDS(wTO_HepaRG,file='example/wTO_HepaRG_network.rds')