# server.R

library(shiny)
library(shinyjs)
options(shiny.maxRequestSize=100*1024^2)   # maximum 100 mb upload limit for the bigger networks

function(input,output,session){
  
  source('functions/general_functions.R')        # harmonic sum, jaccard distance, etc
  source('functions/network_optimizer.R')        # funcitons for network optimization pareto method
  source('functions/pipeline.R')                 # functions for RWR-FGSEA 
  source('functions/plot_functions.R')           # plot functions
  source('functions/vam_functions.R')            # function edc_score  
  source('functions/glm_functions.R')            # elastic net GLM functions
  source('functions/annotation_functions.R')     # annotation functions
  
### Tab 1: Summary -------------------------------------------------------------
  
  output$plot_st1 <- renderPlot({   # pie chart of average edc scores in the summary tab
    data_edcs_count_average <- readRDS("inputData/statistics/average_distribution_compounds_EDCs.rds")
    data_edcs_count_average$type = "Average"
    data_edcs_count_harmonic <- readRDS("inputData/statistics/harmonic_distribution_compounds_EDCs.rds")
    data_edcs_count_harmonic$type = "Harmonic Sum"
    p1 <- pie_chart(data_edcs_count_average, "EDC score")
    p2 <- pie_chart(data_edcs_count_harmonic, "EDC score")
    ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "bottom")
    })
  
  output$plot_st2 <- renderPlot({   # barplot of the number of pathways in summary tab
    data_pathways <- readRDS("inputData/statistics/patway_st.rds")
    print(bar_stat(data_pathways))
    })

  rv_all_f1 <- reactiveValues(min_p = .7, f1_names_colrs = readRDS("inputData/statistics/f1_scores.rds"))
  
  output$plot_st4 <- renderPlot({   #boxplot of k-fold-cross-validation in the summary tab
    print(box_plot(rv_all_f1$f1_names_colrs$f1_scores,
                   rv_all_f1$f1_names_colrs$net_names,
                   rv_all_f1$f1_names_colrs$color_vals,
                   rv_all_f1$min_p))
    })

  observeEvent(input$F1_scores_input,{ # adding a new data F1 layer for comparison
    temp_f1_scores <- rv_all_f1$f1_names_colrs
    f1 <- readRDS(input$F1_scores_input$datapath)
    temp_f1_scores$f1_scores <- rbind(temp_f1_scores$f1_scores, f1)
    rv_all_f1$min_p <- min(temp_f1_scores$f1_scores$values, na.rm = T)
    temp_f1_scores$color_vals <- c(temp_f1_scores$color_vals, input$new_F1_Colr_input)
    temp_f1_scores$net_names <- c(temp_f1_scores$net_names, unique(f1$networks))
    rv_all_f1$f1_names_colrs <- temp_f1_scores
  })

   
# 2.  Tab2: toxicogenomics pipeline----------------------------------------------------------------------
  library(parallel)
  det_cpus<-detectCores()
  
  updateNumericInput(session,'number_cpu_input',
                     value = det_cpus)
  
  rv_pipeline<-reactiveValues(f1_scores=data.frame(values=c(0,0,0,0,0)),
                              status='Toxicogenomics pipeline')
  
  dataset_export<-reactive({
    switch(input$export_input,
         'pathway activation score'=list(NES_scores=rv_pipeline$NES,labels=rv_pipeline$class),
         'All_model_parameters_for_prediction'=rv_pipeline$params,
         'F1_scores'=rv_pipeline$f1_scores[,c('values','networks')],
         'data_frame_GLM_coefs'=rv_pipeline$data_frame,
         'predicted_items'=rv_pipeline$predicted_data)})
   output$export_btn <- downloadHandler(
     filename = function() {
       paste(input$export_input, ".rds", sep = "")
     },
     content = function(file) {
       showNotification("Please wait until the complete message")
       saveRDS(dataset_export(), file)
       showNotification("Saving Finished")
     }) # saving  the output  of the pipeline
  
 
  
   observeEvent(input$pareto_btn,{ 
    if (!is.null(input$network_input$datapath) & !is.null(input$edc_input$datapath) & !is.null(input$decoy_input$datapath) ){
     showNotification("Please wait")
     network<-readRDS(input$network_input$datapath)
     print(input$network_input$name)
     edcs_list<-readRDS(input$edc_input$datapath)
     decoys_list<-readRDS(input$decoy_input$datapath)
     n_cpuc<-input$number_cpu_input
     
       print('Pareto is running')
       genes<-as.numeric(unlist(strsplit(input$gene_set_input,',')))
       edge<-as.numeric(unlist(strsplit(input$edge_set_input,',')))
       edgelist<-cpu_splitter(n_cpuc,edge)
       library(doParallel)
       cl<-makeCluster(length(edgelist))
       registerDoParallel(cl)
       modls<-foreach(i=1:length(edgelist),
                      .export = 'network2silhouette',
                      .packages=c('dnet','igraph','cluster')) %dopar% network2silhouette(network,
                                                                                         edgelist[[i]],
                                                                                         genes,
                                                                                         edcs_list,
                                                                                         decoys_list)
       stopCluster(cl)
       #pareto
       silhouette_mat<-as.data.frame(do.call(rbind,modls))
       silhouette_mat$edge<-rownames(silhouette_mat)
       paret_input<-tidyr::gather(silhouette_mat,key='genes',value='silhouette',-edge)
       paret_input$genes<-as.numeric(paret_input$genes)
       paret_input$edge<-as.numeric(paret_input$edge)
       library(rPref)
       ps<-psel(paret_input,high(silhouette)*low(genes)*low(edge))
       # winning gene and edge combination based on pareto
       n_genes<-ps$genes[which(ps$silhouette==max(ps$silhouette))]
       n_edges<-ps$edge[which(ps$silhouette==max(ps$silhouette))]
       #
       updateTextInput(session,'final_gene_input',value = n_genes)
       updateTextInput(session,'final_edge_input',value = n_edges)
       print(paste('The optimum number of genes is:',n_genes))
       print(paste('The optimum number of edges is:',n_edges))
       print('finished')
       rv_pipeline$status<-'Pareto is finished'
       showNotification("Finished")
       rv_pipeline$pareto<-ps}else{
         showNotification('The fields network,edc list, decoy list should not be empty',type = 'error')}
     })  # Pareto solution to edge and gene size
   
   
   observeEvent(input$network_btn,{ 
     if (!is.null(input$network_input$datapath) & !is.null(input$edc_input$datapath) & !is.null(input$decoy_input$datapath) ){
    showNotification("Please wait")
    net_pat<-input$network_input$datapath
    network<-readRDS(net_pat)
    print(input$network_input$name)
    edcs_list<-readRDS(input$edc_input$datapath)
    decoys_list<-readRDS(input$decoy_input$datapath)
    n_cpuc<-input$number_cpu_input
    n_genes<-as.numeric(input$final_gene_input)
    n_edges<-as.numeric(input$final_edge_input)
    withProgress(message = 'wait',value = 0,{ 
    print('Wait for randomwalk and gene set enrichment analysis')
  
     # RWR-FGSEA
     patways<-readRDS('inputData/pathways.rds')
     mies<-c(edcs_list,decoys_list)
     library(BiocParallel);library(stats)
     if (length(grep(x=as.character(Sys.info()),pattern = 'windows',ignore.case = T))){
       p<-bpstart(SnowParam(n_cpuc))}else{
       p<-bpstart(MulticoreParam(n_cpuc))
       }
     fgsea_res<-pipeline(network,patways,n_edges,n_genes,mies,p)
     bpstop(p)
     print('Gene set enrichment analysis is finished')
     incProgress(1/3,detail = 'RWR-Geneset enrichment analysis is finished')
     # pretreatment of NES values
     fgs<-t(fgsea_res$NES)
     fgs[is.na(fgs)]<-0
     class<-rep('decoy',nrow(fgs))
     names(class)<-rownames(fgs)
     class[which(names(class) %in% names(edcs_list))]<-'edc'
     rv_pipeline$NES<-fgs
     rv_pipeline$class<-class
     
     # ML part training model
     modl<-elastic_net_model_two_class(fgs,class)
     print('Machine Learning is finished')
     incProgress(1/3,detail = 'Machine learning is finished')
     
     # Cross validtion
     print('starting kfold CV validation')
     all_index<-caret::createMultiFolds(class,k=input$k_input,times=input$repeat_input)
     cpu_ind<-cpu_splitter(n_cpuc,1:length(all_index))
     library(doParallel)
     cl<-makeCluster(length(cpu_ind))
     registerDoParallel(cl)
     F1_scores<-foreach(j=1:length(cpu_ind),
                        .export = c('splitted_CV','elastic_net_model_two_class','confusion_matrix','cm_based_evaluation'),
                        .packages=c('caret','glmnet')) %dopar% splitted_CV(cpu_ind[[j]],all_index,fgs,class)
     stopCluster(cl)
     F1_scores<-do.call(c,F1_scores)
     rv_pipeline$f1_scores<-data.frame('values'=F1_scores)
     rv_pipeline$f1_scores$networks<-input$job_name_input
     print(F1_scores)
     
     # finalization of data frame and updating the putative pathways tab
     glm_cof<-modl$cof[-which(rownames(modl$cof) %in% '(Intercept)'),]
     glm_dta_frame<-as.data.frame(glm_cof)
     glm_dta_frame$pathways<-rownames(glm_dta_frame)
     nes<-fgs[which(class %in% 'edc'),]
     mean_nes<-apply(nes, 2, mean)
     mean_nes<-as.data.frame(mean_nes)
     mean_nes$pathways<-rownames(mean_nes)
     all_df<-merge(mean_nes,glm_dta_frame)
     annotations<-pathway_annotate(all_df$pathways)
     all_df$pathways<-annotations
     all_df$network<- input$job_name_input
     all_df<-all_df[,c("network","mean_nes","glm_cof","pathways")]
     colnames(all_df)<- c("network" ,"Average_NES_EDC", "glm_coefs" ,"annotated_pathways")
     rv_pipeline$data_frame<-all_df
     rv_bubble_plot$networks<-c(rv_bubble_plot$networks,input$job_name_input)
     updateSelectInput(session,'data_layer_input',choices = rv_bubble_plot$networks)
     rv_bubble_plot$data_glm<-rbind(rv_bubble_plot$data_glm,all_df)
     # 
     rv_pipeline$status<-'RWR-FGSEA-GLM is finished'
     incProgress(1/3,detail = 'Kfold_CV is finished')
     })
     showNotification("Finished")
     rv_pipeline$params<-list(n_genes=n_genes,n_edges=n_edges,network=network,model=modl)}else{
       showNotification('The fields network,edc list, decoy list should not be empty',type = 'error')}
      }) #   RWR-FGSEA-GLM
   
  
   
   observeEvent(input$predict_btn,{
     if (!is.null(input$all_model_parameters_input$name))rv_pipeline$params<-readRDS(input$all_model_parameters_input$datapath)
     # preparation
     showNotification("Please wait")
     withProgress(message = 'wait',value = 0,{ 
     network<-rv_pipeline$params$network
     model<-rv_pipeline$params$model
     n_cpuc<-input$number_cpu_input
     n_genes<-as.numeric(rv_pipeline$params$n_genes)
     n_edges<-as.numeric(rv_pipeline$params$n_edges)
     print(n_genes);print(n_edges)
     showNotification('wait for random walk and FGSEA')
    
     # RWR-FGSEA
     patways<-readRDS('inputData/pathways.rds')
     mies<-readRDS(input$test_compounds_input$datapath)
     library(BiocParallel);library(stats)
     if (length(grep(x=as.character(Sys.info()),pattern = 'windows',ignore.case = T))){
       p<-bpstart(SnowParam(n_cpuc))}else{
       p<-bpstart(MulticoreParam(n_cpuc))}
     fgsea_res<-pipeline(network,patways,n_edges,n_genes,mies,p)
     bpstop(p)
     print('Gene set enrichment was finished, wait')
     incProgress(1/2,detail = 'RWR-FGSEA is finished')
     
     #prediction 
     fgs<-t(fgsea_res$NES)
     fgs[is.na(fgs)]<-0
     predicted_class<-predict(model$modl,fgs,type='raw')
     predicted_prob<-predict(model$modl,fgs,type='prob')
     rv_pipeline$predicted_data<-as.data.frame(list(rownames(fgs),predicted_prob[,'edc'],predicted_class))
     colnames(rv_pipeline$predicted_data)<-c('names','EDC_probability','class')
     rv_pipeline$status<-'Prediction of newcompounds is finished'
     incProgress(1/2,detail = 'Prediction is finished')
     showNotification("Finished")
     })
   }) #predcition of new compounds
   
   
   # outputs
   output$prediction_tab1<-renderTable({rv_pipeline$predicted_data}) #predicted compounds table
   
   output$pareto_tabl<-renderTable({rv_pipeline$pareto}) #predicted compounds table
   
   output$status_lbl<-renderText({print(rv_pipeline$status)}) # status verbatim update
   
   output$f1_plt<-renderPlot({barplot(rv_pipeline$f1_scores$values,main = 'K-Fold-CV',ylab = 'F1 Scores')}) # k-fold-CV plot
  

   
   
  
# 3.  Tab3: putative pathways-------------------
   networks<-readRDS('inputData/network_names.rds')     # names of the networks
   glm_coefs<-readRDS('inputData/GLM_coeffs_edcs_decoys.rds')
   updateSelectInput(session,'data_layer_input',choices = networks)
   rv_bubble_plot<-reactiveValues(data_glm=glm_coefs,
                                  data_subset=glm_coefs,
                                  networks=networks,
                                  networks_selected=networks)
   ranges_ptw <- reactiveValues(x = NULL, y = NULL)
   
  
  observeEvent(input$ptw_dbl_click, {brush_adjust(input$ptw_brush,ranges_ptw)}) #plot fucntion call
 
   
  output$export_btn_pathways <- downloadHandler(
    filename = function() {
      'pat.csv'
    },
    content = function(file) {
      write.csv(rv_bubble_plot$data_subset, file)
    }
  )  # export pathways as bubble plot
  
  
  observeEvent(input$calc2,{                   # pressing the show button on mode of action Tab
  data_P<-rv_bubble_plot$data_glm
  selected_networks<-input$data_layer_input   # selected networks by user
  selected_categories<-input$pathway_category_input # selected cateogires by user
  data_P<-data_P[grep(x=data_P$annotated_pathways,
          pattern=paste(selected_categories,'_',sep='',collapse="|")),]  # pathway category subsetting
 
  data_P<-data_P[data_P$network %in% selected_networks,]# Network subsetting
  data_P<-data_P[which(data_P$glm_coefs>input$GLM_coef &
                      data_P$Average_NES_EDC>input$NES),] # cutoffs subsetting
  rv_bubble_plot$data_subset<-data_P
  rv_bubble_plot$networks_selected<-selected_networks
  rv_bubble_plot$categories_selected<-selected_categories }) #observe event
  
  
  observeEvent(input$new_data_glm_input,{ # adding a new data layer
    new_data_layer<-readRDS(input$new_data_glm_input$datapath)
    if(is.data.frame(new_data_layer)){
    if (ncol(new_data_layer)==4){
    colnames(new_data_layer)<-colnames(rv_bubble_plot$data_glm)
    rv_bubble_plot$networks<-c(rv_bubble_plot$networks,unique(new_data_layer$network))
    updateSelectInput(session,'data_layer_input',choices = rv_bubble_plot$networks)
    rv_bubble_plot$data_glm<-rbind(rv_bubble_plot$data_glm,new_data_layer)}}
    else{showNotification('The data frame does not meet the criteria')}
    })
  
  output$plt_ptway<-renderPlot({
  bubble_plot(rv_bubble_plot$data_subset,rv_bubble_plot$networks_selected,ranges_ptw)}) # plot_functions call
  


# 4.  Tab4: Edc scores and class probabilties  -----------------------------
  # params and variables
  all_vam<-readRDS('inputData/all_edc_scores.rds')     # edc scores and class probs for networks
  networks<-readRDS('inputData/network_names.rds')     # names of the networks
  network_score_levelss=c(networks,'average_edc_score','harmonic_sum_edc_score')
  
  all_possible_mesh_cas_names=c(unique(all_vam$mesh),
                                unique(all_vam$cas),
                                unique(all_vam$comp_names))
  
  rv_edc_library<-reactiveValues(all_vam=all_vam,
                                 all_possible_mesh_cas_names=all_possible_mesh_cas_names)
  
  rv_edc_score<-reactiveValues(table_scors=c(),
                               class_prob_scores=edc_score(all_vam,'1962-83-0',
                                                           network_score_levelss)[1:24,],
                               harmonic_average_scores=edc_score(all_vam,'1962-83-0',
                                                           network_score_levelss)[25:26,]) #Default compound at launch
  updateSelectizeInput(session, "cmpname", selected = '1962-83-0',choices  =all_possible_mesh_cas_names)
  updateSelectInput(session,"edc_score_layer_input",
                    choices = networks,selected = networks[1:24])
  ranges_class_prob <- reactiveValues(x = NULL, y = NULL)
  ranges_edc_score <- reactiveValues(x = NULL, y = NULL)
  path_2_network<-'large_file/all_precompiled_pipeline.RDSS'  # all networks and paramterers and models
  if(!file.exists(path_2_network))shinyjs::hide('mie2classprob_btn')
  # events 
  observeEvent(input$plot_class_prob_scores_dbl_click,{
               brush_adjust(input$class_prob_scores_brush, ranges_class_prob)}) #plot fucntion call
  
  # data before multi compate plot are calculated
  observeEvent(input$calc,{  
  if (!is.null(input$cmpname)){
  data_for_edc_score_plot<-edc_score(rv_edc_library$all_vam,
                                     input$cmpname,
                                     network_score_levelss,
                                     col.scores = input$edc_score_layer_input)
  rv_edc_score$class_prob_scores<-data_for_edc_score_plot[data_for_edc_score_plot$network %in% input$edc_score_layer_input,]
  rv_edc_score$harmonic_average_scores<-data_for_edc_score_plot[data_for_edc_score_plot$network %in% c('average_edc_score',
                                                                                                       'harmonic_sum_edc_score'),]

table_data<-data_for_edc_score_plot[data_for_edc_score_plot$network %in% 'average_edc_score'  ,c('nnames','score')]
table_data$harmonic<-data_for_edc_score_plot$score[data_for_edc_score_plot$network %in% 'harmonic_sum_edc_score']
colnames(table_data)<-c('Compounds','Average','Harmonic_Sum')
rv_edc_score$table_scors<-table_data
  }else{showNotification('Enter at least one compound for visualization',type = 'error')}
   })   # show on plot BUTTON
  # plot
  
  output$plot_class_prob_scores<-renderPlot({edc_multi_compare_bar(rv_edc_score$class_prob_scores,
                                                                    network_score_levelss[factor(networks,
                                                                                                 levels = networks)],
                                                                    ranges_class_prob,
                                                                    'Class Prob.',
                                                                    c(0,1),angle_t = 60)}) #plot_functions call
  output$table_edc_scores<-renderTable({rv_edc_score$table_scors}) # table
  
   # output$plot_edc_score<-renderPlot({edc_multi_compare_bar(rv_edc_score$harmonic_average_scores,
   #                                                          network_score_levelss[network_score_levelss %in% c('average_edc_score',
   #                                                                                                             'harmonic_sum_edc_score')],
   #                                                          ranges_edc_score,
   #                                                          'edc score',
   #                                                          c(0,2),angle_t = 0,show_leg = F)}) #plot_functions call

  
   # mie 2 class prob
  observeEvent(input$mie2classprob_btn,{
    if(any(tolower(trimws(input$txt_input_newcompound_name)) %in% tolower(rv_edc_library$all_possible_mesh_cas_names)) %in% F &
       (!trimws(input$txt_input_newcompound_name)=='')){ 
      showNotification('wait',duration = NULL)
      newcomp<-strsplit(input$txt_input_mies,',') #genes as MIEs
      if (!"all_pipeline_precompiled_data" %in% ls())all_pipeline_precompiled_data<-readRDS(path_2_network)
      results<-onecompound_mies2classprob(all_pipeline_precompiled_data$networks,
                                          all_pipeline_precompiled_data$pathways,
                                          all_pipeline_precompiled_data$models,
                                          newcomp)
      if(length(which(!is.na(results) %in% T))>0){
        results<-as.data.frame(t(results),col.names=names(results))
        results$harmonic_sum_edc_score<-results$average_edc_score<-0
        results$comp_names<-trimws(input$txt_input_newcompound_name)
        results$mesh<-results$cas<-''
        results$is_in_training<-'unk'
        results<-results[,colnames(rv_edc_library$all_vam)]
        rv_edc_library$all_vam<-rbind(rv_edc_library$all_vam,results)
        data_for_edc_score_plot<-edc_score(rv_edc_library$all_vam,
                                           results$comp_names,
                                           network_score_levelss,
                                           col.scores = input$edc_score_layer_input)
        rv_edc_score$class_prob_scores<-data_for_edc_score_plot[data_for_edc_score_plot$network %in% input$edc_score_layer_input,]
        rv_edc_score$harmonic_average_scores<-data_for_edc_score_plot[data_for_edc_score_plot$network %in% c('average_edc_score',
                                                                                                             'harmonic_sum_edc_score'),]
        print(rv_edc_score$harmonic_average_scores)
        rv_edc_library$all_possible_mesh_cas_names<-c(results$comp_names,rv_edc_library$all_possible_mesh_cas_names)
        updateSelectizeInput(session, "cmpname", selected = results$comp_names,choices  =rv_edc_library$all_possible_mesh_cas_names)
        
        showNotification('finished')}else{showNotification('Wrong MIEs')}
      
    }else{showNotification('the compounds already exists in the data base or the field is empty')}
    
  })  #MIE 2 edc score
  output$export_btn_edcscores <- downloadHandler(
    filename = function() {
      'pat.csv'
    },
    content = function(file) {
      write.csv(rbind(rv_edc_score$class_prob_scores,rv_edc_score$harmonic_average_scores), file)
    }
  )  
    
    
# observeEvent(input$chkbox_most_informatiave,{
#  if (input$chkbox_most_informatiave==T){
#   updateSelectInput(session,"edc_score_layer_input",choices = networks,selected = networks[c(1,2,4,13,14,15)])}})

# 5.  Tab5: Toxpi evaluation ------------------------------------------------------------------
rv_eval_toxpi<-reactiveValues(plot_data=readRDS('inputData/toxpi_scores.rds'),
                              selected_for_score=networks)
  
updateSelectInput(session,'toxpi_layer_input',choices = networks)
ranges_toxpi <- reactiveValues(x = NULL, y = NULL)

  
 # observeEvent(input$toxpi_plot_dbl_click, {brush_adjust(input$toxpi_brush,ranges_toxpi)}) #plot fucntion call
# refresh button of ToxPi
observeEvent(input$toxpiBtn_refresh,{
# if(!is.null(input$toxpi_layer_input)){
#   updateSelectInput(session,'toxpi_layer_input',selected =input$toxpi_layer_input )}
  rv_eval_toxpi$plot_data<-readRDS('inputData/toxpi_scores.rds')
  })
  observeEvent(input$toxpi_btn,{
    if(!length(input$toxpi_layer_input)==0){ # at least one layer should be selected by the user
    primary_data<-rv_eval_toxpi$plot_data
    rv_eval_toxpi$selected_for_score<-input$toxpi_layer_input
    primary_data$unk_VAM_scores<-mean_function(primary_data,rv_eval_toxpi$selected_for_score)
    primary_data<-primary_data[!is.na(primary_data$unk_VAM_scores),]
    rv_eval_toxpi$plot_data<-primary_data}})  #Show  button click
  
  output$plot_toxpi<-renderPlot({
    rv_eval_toxpi$mintoxpi=input$toxpi_plot_click$y-input$slider_click_toxpi
    rv_eval_toxpi$maxtoxpi=input$toxpi_plot_click$y+input$slider_click_toxpi
    rv_eval_toxpi$minedc=input$toxpi_plot_click$x-input$slider_click_toxpi
    rv_eval_toxpi$maxedc=input$toxpi_plot_click$x+input$slider_click_toxpi
    all_mat<-rv_eval_toxpi$plot_data
    res<-all_mat[which(all_mat$toxpi<=rv_eval_toxpi$maxtoxpi               &
                         all_mat$toxpi>rv_eval_toxpi$mintoxpi              &
                         all_mat$unk_VAM_scores<= rv_eval_toxpi$maxedc     &
                         all_mat$unk_VAM_scores>rv_eval_toxpi$minedc),c("compound_name" ,
                                                                        "cas","toxpi",
                                                                        "unk_VAM_scores" )]
    ind_c<-which(all_mat$compound_name %in% res$compound_name)
    rv_eval_toxpi$plot_data$X[ind_c]<-rv_eval_toxpi$plot_data$compound_name[ind_c]
 
    toxpi_plot(rv_eval_toxpi$plot_data,
              ranges_toxpi,
              min_toxpi=input$slider_toxpi_plt_toxpi_cutoff,
              min_edc_score=input$slider_toxpi_plt_edc_score_cutoff)}) #plot_functions call
  
 
 # data table bellow the plot
   output$table_toxpi<-renderDataTable({
     p_data<-rv_eval_toxpi$plot_data
     ind_comps<-unique(p_data$X)
     ind_comps[!ind_comps %in% '']
     res<-p_data[which(p_data$compound_name %in% ind_comps),
                                        c("compound_name" ,"cas","toxpi","unk_VAM_scores",'status' )]
     
     colnames(res)<-c("compound_name" ,"cas","toxpi","Average_EDC_score",'type' )
     res$link<-sprintf('<a href="https://comptox.epa.gov/dashboard/dsstoxdb/results?abbreviation=TOXCAST_PH2&search=%s" target="_blank" class="btn btn-primary">link_comptox</a>',res$cas)
     res<-res[,c("compound_name","toxpi","Average_EDC_score",'link','type')]
     res$toxpi<-round(res$toxpi,3)
     res$Average_EDC_score<-round(res$Average_EDC_score,3)
    
     return(res)
     },escape=F,options = list(
       pageLength = 5))
  
  
  output$export_btn_toxpi <- downloadHandler(
    filename = function() {
      'pat.csv'
    },
    content = function(file) {
      write.csv(rv_eval_toxpi$plot_data, file)
    }
  )  # export pathways as bubble plot
  
  
  # output$txt_toxpi_click<- renderText({
  # rv_eval_toxpi$mintoxpi=input$toxpi_plot_click$y-input$slider_click_toxpi
  # rv_eval_toxpi$maxtoxpi=input$toxpi_plot_click$y+input$slider_click_toxpi
  # rv_eval_toxpi$minedc=input$toxpi_plot_click$x-input$slider_click_toxpi
  # rv_eval_toxpi$maxedc=input$toxpi_plot_click$x+input$slider_click_toxpi
  #   res2<-paste('Toxpi= ',round(as.numeric(input$toxpi_plot_click$y),2),
  #              '  +/- ',input$slider_click_toxpi ,
  #              '\n','EDC_score= ',
  #              round(as.numeric(input$toxpi_plot_click$x),2),'  +/-  ',
  #              input$slider_click_toxpi ,sep='')
  # 
  #   return(res2)
  #  })
  
} # server function

#Updated by: Arindam Ghosh