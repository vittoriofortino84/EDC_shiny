#ui.R
# libraries and menues---------------------------------------------------------------


library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
library(DT)
shinyUI(
  dashboardPage(
    dashboardHeader(title="EDCmet"),
    dashboardSidebar(
      sidebarMenu(id="tabs",
                  menuItem("Summary", tabName = "dashboard", icon = icon("dashboard")),
                  menuItem('Toxicogenomics Pipeline',tabName = 'p1',icon = icon('th')),
                  menuItem("Pathway activation scores",tabName="p2",icon = icon('th')),
                  menuItem("Predcited EDC scores",tabName="p3",icon = icon('th')),
                  menuItem("Evaluation with ToxPi ",tabName="p4",icon = icon('th'))
                #  menuItem("Prediction from MIEs ",tabName="p5",icon = icon('th'))
      )
    ),
    dashboardBody(
      tabItems(

        #  # 1. tab summary-----------------------------------------------------------------------

 tabItem(tabName = "dashboard",
         fluidRow(
           box(plotOutput("plot_st1",height = 250),width = 6,height =300),
           box(plotOutput("plot_st2",height = 250),width = 6,height = 300)
         ),
         fluidRow(
        # box(plotOutput("plot_st3",height = 400),width = 12,height = 500),
           box(plotOutput("plot_st4",height = 400),
               bsButton("qtf1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
               bsPopover(id = "qtf1", title = "F1_scores",
                         content = paste0('a data frame as RDS file with two columns. The first column should be named values (the F1 scores)',
                                          '. The second column should be named networks (a similar name for all rows i.e example/outPut/F1_scores.rds).'),
                         placement = "right", 
                         trigger = "hover",
                         options = list(container = "body")),
               fluidRow(box(
                 fileInput('F1_scores_input',label = 'Add new F1-scores',accept = c("rds","A list",".rds")),width = 6,solidHeader = T),
               box(textInput("new_F1_Colr_input", "Color for the new layer", value = "orange"),width = 6,solidHeader = T)

               ),
               width = 12,height = 580)
       
           
         )
),        
        #  # 2. tab Toxicogenomics pipeline -----------------------------------------------
        
        tabItem(tabName = 'p1',
                headerPanel('Network Optimization, RWR, FGSEA,elastic net GLM'),
                mainPanel(
                fluidRow(
                 box(tags$script('
                             $(document).ready(function(){
                             var d = new Date();
                             var target = $("#clientTime");
                             target.val(d.toLocaleString());
                             target.trigger("change");
                             });
                             '),
                 textInput("clientTime", "Client Time", value = ""),
                 numericInput('number_cpu_input','number of cpus',4,min=1,step = 1),
                 bsButton("qtjob", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                 bsPopover(id = "qtjob", title = "Job name",
                           content = paste0('Please select a proper job name based on the gene-gene network.  ',
                                            'The pathways activation score tab and other genereated parameters use this name for comparison'),
                           placement = "right", 
                           trigger = "hover",
                           options = list(container = "body")),
                 textInput('job_name_input',label = 'Job name',value = 'my_job'),
                 verbatimTextOutput('status_lbl'), br()
                 ,height = 400,width = 6,solidHeader = T),
                
                 box(
                   bsButton("qt1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qt1", title = "Network",
                             content = paste0('A dataframe as rds file with three columns, the first two columns are gene entrez IDs ',
                                              ', the third column can be the weights from wTO package or weights from WGCNA. The maximum size should not exceed 100mb.'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('network_input',label = 'Gene-Gene network',accept = c(
                     "rds","A data frame with three columns gene1 gene2 weight/p_value",".rds")),
                   bsButton("qtedc", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtedc", title = "EDCs",
                             content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                              'to that EDC as ENTREZ gene IDs'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('edc_input',label = 'A list of EDCs and MIEs',accept = c("rds","A list", ".rds")),
                   bsButton("qtdec", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtdec", title = "Negative controls",
                             content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                              'to that negative control as ENTREZ gene IDs'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('decoy_input',label = 'A list of negative controls and MIEs',accept = c("rds","A list",".rds"))
                   ,height = 400,width = 6,solidHeader = T)
                 ), 
                box(
                  fluidRow(
                   box(textInput('gene_set_input',label = 'Number of top Genes after Random walk :',value = '200,500,1000'),width = 6,solidHeader = T),
                   box(textInput('edge_set_input',label = 'Network Edge proportion before random walk:',value = '0.02,0.03,0.05,0.1'),width = 6,solidHeader = T)
                  ),
                   actionButton('pareto_btn',label = 'Use Pareto solution for optimization'), 
                   tableOutput('pareto_tabl')
                   ,width = 12,collapsible = T,collapsed = T,title = 'Pareto based optimization of the pipeline'
                  ),
                box(
                  fluidRow(
                   box( bsButton("qfinal_gene", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                        bsPopover(id = "qfinal_gene", title = "Final Gene set for enrichment",
                        content = paste0('The maximum number of sorted genes based on visit probability to be selected for gene set enrichment analysis after randomwalk',
                        '. In case of running pareto solution this parameter will be automatically updated'),
                                                    placement = "right", 
                                                    trigger = "hover",
                                                    options = list(container = "body")),
                     textInput('final_gene_input',label = 'Final Number of genes after random walk:',value = '500'),width = 6,solidHeader = T),
                   box(bsButton("qfinal_edge", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                        bsPopover(id = "qfinal_edge", title = "Edge percent from the network",
                                  content = paste0('The portion of sorted edges to be extracted from the network ',
                                                   '. In case of running pareto solution this parameter will be automatically updated'),
                                  placement = "right", 
                                  trigger = "hover",
                                  options = list(container = "body")),
                     textInput('final_edge_input',label = 'Final Edge  proportion from the network:',value = '.05'),width = 6,solidHeader = T)
                  ),
                  fluidRow(
                   box(numericInput('k_input','K folds Cross Validation',5,min=2,max = 20,step = 1),width = 6,solidHeader = T),
                   box(numericInput('repeat_input','Cross validation repeats',1,min=1,max = 10,step = 1),width = 6,solidHeader = T)),
                   actionButton('network_btn',label = 'Run RWR-FGSEA-GLM'),
                   plotOutput('f1_plt')
                ,width = 12,collapsed = T,collapsible = T,title = 'RWR-Gene set enrichment-GLM'),
                box(
                   selectInput('export_input', 'Select a data set to export as rds file:',
                               choices = c('pathway activation score',
                                           'All_model_parameters_for_prediction',
                                           'F1_scores','data_frame_GLM_coefs',
                                           'predicted_items')),
                   downloadButton('export_btn','Export selected item'),width = 12,title = 'Export',collapsible = T,collapsed = T),
                box(
                  fluidRow(
                   box(fileInput('test_compounds_input',label = 'A list of MIEs for unknown compounds',accept = c("rds","A list",".rds")),width = 6,solidHeader = T),
                   box(fileInput('all_model_parameters_input',label = ' Model and paramters for prediction',accept = c("rds","A list",".rds")),width = 6,solidHeader = T)),
                   actionButton('predict_btn',label = 'Predict Unknown compounds'),
                   tableOutput('prediction_tab1'),width = 12,collapsible = T,title = 'prediction of new compounds',collapsed = T)
                                                      #end of side bar panel2
                 ,width = 12)),
        
        
        #  # 3. tab Putative pathways -----------------------------------------------
        
        
        tabItem(tabName = 'p2',
                 headerPanel('Putative pathways as the mode of action for EDCs'), 
                 sidebarPanel(
                   bsButton("q2", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q2", title = "GLM coefficients",
                             content = paste0('Set the minimum cutoff for GLM coefficient'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'GLM_coef',
                               label='GLM Coefficients',
                               min = 0,max=1,value=0,step = .01,round = F),
                   bsButton("q3", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q3", title = "Pathway activation score",
                             content = paste0('Set the minimum cutoff for the average of enrrichment pathway activation scores for EDCs'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'NES',
                               label=' pathway activation scores',
                               min = 0,max=1.7,value=0,step = .01,round = F),
                   hr(),
                   selectInput('pathway_category_input', 'Pathway category (Hold Ctrl to select multi categories)', 
                               c('KEGG','REACTOME','WIKI','GO'), multiple=TRUE, selectize=FALSE),
                  
                   hr(),
                   selectInput('data_layer_input', 'Data layers (Hold Ctrl to select multi data layers)', state.name, multiple=TRUE, selectize=FALSE),
                   actionButton(inputId = 'calc2', label = 'Show'),hr(),br(),
                   bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q1", title = "data structure",
                             content = paste0('A data frame as rds file with four columns is needed, the first column should be ', 
                                              'the name of data layer ','the second column should be the average of pathway activations scores ',
                                              'across all EDCs for that pathway , the third column is the GLM coefficients for the pathway ',
                                              'and the last column is the name of the pathway. You can select any ',
                                              'names for the columns but the order should be as above'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('new_data_glm_input',label = 'Add a new data layer',accept = c(
                     "rds", 
                     "A list",
                     ".rds")
                   )
                 
                   
                   
                 ),                                   #end of side bar panel2
                 mainPanel(  
                   plotOutput('plt_ptway',
                              height = 900,
                              dblclick = 'ptw_dbl_click',
                              brush = brushOpts(id = 'ptw_brush',
                                                resetOnNew = T)),
                   downloadButton('export_btn_pathways','Export Plot data as csv'),br(),hr()
                   
                 )                                   # end of main panel2
        ),                                           # end tabpanel 2 prediction form MIEs
        
        
        #  # 4. tab Class probabities for EDCs -----------------------------------------
       
        tabItem(tabName = 'p3',         
                 headerPanel('Compound EDC Probability'), 
                 sidebarPanel( 
                   bsButton("qtedc1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtedc1", title = "EDC_scores",
                             content = paste0('Selected data layers will be used for the prediction of ', 
                                              'average and harmonic sum scores. Hold Ctrl key to select ',
                                              'mutli data layers. The class probabilites of',
                                              ' maximum 5 compounds can be compared.'
                                              ),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   selectizeInput('cmpname','CAS, MESH ID, Compound name',
                                  choices=c('1962-83-0'),
                                  options = list(maxOptions = 10,maxItems=5)),
                  selectInput('edc_score_layer_input', 'Data layers', c('PPI_STRINGdb'),
                               multiple=TRUE, selectize=FALSE),
                   #checkboxInput('chkbox_most_informatiave',
                   #'Select the most correlated networks with ToxCast',F),   # most informative 

                   actionButton(inputId = 'calc',
                                label = 'Show on the plot'),br(),hr(),
                  downloadButton('export_btn_edcscores','Export Plot data as csv'),
                   # helpText("Note: You can compare the class probabilities ",
                   #          "as well as the EDC scores for maximum 5 compounds ")
                       hr(),
                   bsButton("qtfile", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtfile", title = "Precompiled networks",
                             content = paste0('The MIEs (Genes) of one compound will be used to ',
                                              'predict class probabilities ',
                                              'on the 24 data layers. ',
                                              'For this module to run the file all_precompiled_pipeline.RDSS  ',
                                              'should exist in the folder large_file'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   textInput("txt_input_newcompound_name", "Name of the compound",
                             value = "new_compound"),

                   textInput("txt_input_mies", "MIES of one compound", 
                             value = "4617,4654,4656,5077,25937"),
                   useShinyjs(),
                   actionButton(inputId = 'mie2classprob_btn',label = 'Calculate from MIEs'),hr()
                 ),                                         # End of side bar panel 1
                 mainPanel(
                  box(tableOutput('table_edc_scores'),collapsible = T,collapsed = T,width = 12),
                  box(
                   plotOutput('plot_class_prob_scores',
                              height = 500,
                              dblclick = 'plot_class_prob_scores_dbl_click',
                              brush = brushOpts(id = 'class_prob_scores_brush',
                                                resetOnNew = T)),width = 12,collapsible = T)
                 
                 
                  
                  #plotOutput('plot_edc_score',
                  #            height = 300),
                   
                  
                 )                                           # End of main panel 1
        ),


        #  # 5. tab Evaluation with toxpi ----------------------------------------------


       tabItem(tabName = 'p4',         
              headerPanel('Evaluation with ToxPi Scores'), 
              sidebarPanel(
                bsButton("qt4_1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                bsPopover(id = "qt4_1", title = "Data Layers",
                          content = paste0('Selected layers will be used to calculate the average EDC scores, ',
                                            'hold Ctrl to select multi layers',
                                           ', In case that one layer is selected the class probability of ', 
                                           'that layer will be used for comparison'),
                          placement = "right", 
                          trigger = "hover",
                          options = list(container = "body")),
                selectInput('toxpi_layer_input', 'Data layers', c('PPI_STRINGdb'), multiple=TRUE, selectize=FALSE),
                sliderInput(inputId = 'slider_toxpi_plt_toxpi_cutoff',
                            label='Minimum Toxpi score:',
                            min = 0.000,max=.4,value=.1,step = .01,round = F),
                sliderInput(inputId = 'slider_toxpi_plt_edc_score_cutoff',
                            label=' Minimum EDC score:',
                            min = 0.000,max=1,value=.8,step = .01,round = F),
                sliderInput(inputId = 'slider_click_toxpi',
                            label=' Plot Click prescision ',
                            min = 0.000,max=0.05,value=0.005,step = .0001,round = F),
                actionButton(inputId = 'toxpi_btn',label = 'Calculate'),br(),hr(),
                actionButton(inputId = 'toxpiBtn_refresh',label = 'Refresh'),br(),hr(),
                downloadButton('export_btn_toxpi','Export Plot data as csv')
               # verbatimTextOutput("txt_toxpi_selected_layers",placeholder = F)
                ),
            mainPanel(                    
          box(plotOutput('plot_toxpi',click = "toxpi_plot_click"                             
                   # dblclick = 'toxpi_plot_dbl_click',
                    #brush = brushOpts(id = 'toxpi_brush',
                    #                  resetOnNew = T)
                    ),collapsible = T,width = 12),
          verbatimTextOutput("txt_toxpi_click"),
          # verbatimTextOutput("txt_toxpi_compounds"),
          #tags$head(tags$style("#txt_toxpi_compounds{color:blue; font-size:14px; font-style:bold; overflow-y:scroll; max-height: 100px; background: ghostwhite;}")),
          dataTableOutput('table_toxpi')
        )                                           # End of main panel 1
       ), # end of tab item toxcast

        #  # 6. mie 2 class prob -----------------------------------------------------


     tabItem(tabName = 'p5',         
        headerPanel('Predcition of class probabilities from MIEs'), 
        sidebarPanel(

          
          
          
        ),
        mainPanel(
          
          
          
          
          
        ) 
      ) # end of prediction from mIEs
 
      #ANY new tabs here

        )                           # End of all tabsets panel
    )
  )
)







