st<-readRDS('inputData/GLM_coeffs_edcs_decoys.rds')
all_pathways<-unique(st$annotated_pathways)
pat_types<-c('KEGG','REACTOME','GO','WIKI')
pat_sts<-rep(NA,length(pat_types))
for (i in 1:length(pat_types))  pat_sts[i]<-length(grep(x=all_pathways,pattern = pat_types[i]))
patwa_statistics<-as.data.frame(list(pat_sts,pat_types))
colnames(patwa_statistics)<-c('number','type')

ggplot(patwa_statistics, aes(x=type,
                            y=number,
                                    fill=type)) +
  geom_bar(stat="identity",show.legend = F)+
  ylab('Number of Pathways')+xlab('Type')+
  theme_minimal()+
  theme(legend.position = 'bottom',
        axis.text.x=element_text(angle = 0 ,hjust = 1 ))
saveRDS(patwa_statistics,file = 'inputData/statistics/patway_st.rds')
