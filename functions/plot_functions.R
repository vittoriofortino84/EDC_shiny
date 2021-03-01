# plots
bar_stat<-function(patwa_statistics){ 
require(ggplot2)
p<-ggplot(patwa_statistics, aes(x=type,
                             y=number,
                             fill=type)) +
  geom_bar(stat="identity",show.legend = F)+
  ggtitle('Pathways Used in the pipeline')+
  ylab('Number of Pathways')+xlab('')+
  theme_minimal()+
  theme(legend.position = 'bottom',
        axis.title.y =  element_text(size = 14),
        axis.text.x=element_text(size = 14,angle = 45 ,hjust = 1 ))
return(p)
}# used in summary tab for barplot of pathways


pie_chart<-function(pie_data,title){
  require(ggplot2)
  ggplot(pie_data, aes(x = 2, y = prop, fill=factor(Category,levels=rev(Category)))) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y")+labs(fill=title)+
    geom_text(aes(y = lab.ypos, label = Number), color = "Black",size=4)+
    scale_fill_manual(values = rev(pie_data$colr))+
    facet_wrap(.~pie_data$type)+
    theme_void()+theme(legend.title = element_text(size = 14),
                       strip.text = element_text(size = 16),
                       legend.text = element_text(size = 14))+
    xlim(0.5, 2.5)
} # ggplot



box_plot<-function(dat,nams,color_values,min_point=0.7){ 
  require(ggplot2)
  p<-ggplot(dat, aes(x=factor(networks,level=nams), y=values,fill=factor(networks,level=nams))) +
    geom_boxplot(show.legend = F) +theme_minimal()+ ggtitle('Accuracy of Data layers')+
    xlab('Data Layers')+
    ylab('F1-score')+scale_fill_manual(values=color_values)+
    labs(x='',fill='Data Layers')+
    scale_y_continuous(breaks = seq(min_point,1,0.05), limits = c(min_point,1)) +
    theme(axis.text.x=element_text(angle = 45,hjust = 1,size = 14),
      axis.title.y = element_text(size = 16),
      axis.ticks.x=element_blank(),legend.position = 'none',
      plot.margin = margin(0, 0, 0, 5, "cm"),
      panel.grid.minor = element_blank())
  return(p)
}#used in summary tab for k-fold-cv results

bubble_plot<-function(data_x,net_names,ranges,Title=''){
  require(ggplot2);require(RColorBrewer)
  ggplot(data_x,aes(x = reorder(annotated_pathways,-Average_NES_EDC),
                    factor(network,level=net_names),size =glm_coefs ))+
    geom_point(shape=21,aes(fill=Average_NES_EDC))  +
    scale_size(range=c(4,10),breaks=c(0.1,0.5,1,1.5),limits=c(0.0,1.5),labels = c('<= 0.1','0.5','1','1.5'))+
    scale_fill_gradientn(colours = brewer.pal(n = 9, name ='YlOrRd'),
                         breaks=seq(0,max(data_x$Average_NES_EDC),0.4),
                         limits=c(0,max(data_x$Average_NES_EDC)))+
    labs(size= 'GLM-Coef',fill='Activation Score',alpha='Pathway CV-Stability')+
    ylab('')+ xlab('')+ guides(fill= guide_colorbar(label.theme = element_text(angle = 0,size=10)  ),
                               size=guide_legend(nrow = 4,order = 2))+ theme_minimal() + 
    coord_cartesian(ranges$x, ylim = ranges$y, expand = T)+
    theme(strip.text.y = element_text(angle=0),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.position  = 'right',
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 25, hjust = 1,size = 10))+
    ggtitle(Title)
}# used in putatuve pathways tab

edc_multi_compare_bar<-function(edc_data,net_lvls,ranges,ylbl,y_scale,angle_t=60,show_leg=T){
  require(ggplot2)
  ggplot(edc_data, aes(x=factor(network,levels = net_lvls),
                                      y=score,
                                      fill=factor(nnames))) +
    geom_bar(show.legend =show_leg ,stat="identity",position = position_dodge2(preserve = 'single'))+ 
    labs(fill='')+
    scale_y_continuous(limits = y_scale)+
    ylab(ylbl)+xlab('')+
    theme_minimal()+
    guides(fill=guide_legend(nrow = 5))+
    coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = T)+
    theme(legend.position = 'bottom',
          plot.margin = margin(0, 0, 0, 5, "cm"),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x=element_text(size = 12,angle = angle_t ,hjust = 1 ))
  
  
}# used in edc score tab

toxpi_plot<-function(all_mat,ranges,min_toxpi=0.1,min_edc_score=0.8){ 
  require(ggplot2)
  all_mat$colr<-'grey50'
  all_mat$colr[all_mat$unk_VAM_scores >=min_edc_score & all_mat$toxpi >=min_toxpi]<-'red'
  ggplot(all_mat,aes(x=unk_VAM_scores,y=toxpi))+
    geom_point(size=3,color=all_mat$colr)+
    geom_text(label=all_mat$X,size=5,segment.size = 0.2, segment.color = "grey",show.legend = F) +
    ylab('TOXPI score')+ xlab('Average EDC score')+
    theme_minimal()+ 
    coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = T)+
    ggtitle('Average EDC score VS TOXPI scores')+xlab('EDC score')+ylab('TOXPI score')
}#used in evaluation with toxpi tab

brush_adjust<-function(brush,ranges){
if (!is.null(brush)) {
  ranges$x <- c(brush$xmin, brush$xmax)
  ranges$y <- c(brush$ymin, brush$ymax)
} else {
  ranges$x <- NULL
  ranges$y <- NULL
}
} # used in zoom for the putative pathways, toxpi and edc score tabs







