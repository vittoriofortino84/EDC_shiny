#ui.R

# Libraries to import ----------------------------------------------------------

library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
library(DT)

# Sidebar menu -----------------------------------------------------------------

shinyUI(
  dashboardPage(
    dashboardHeader(title = "EDTox"),
    dashboardSidebar(
      sidebarMenu(id = "tabs",
                 #menuItem("Home", tabName = "home", icon = icon("dashboard")),
                  menuItem("Summary", tabName = "dashboard", icon = icon("dashboard")),
                  menuItem("Toxicogenomics Pipeline", tabName = "p1", icon = icon('th')),
                  menuItem("Pathway activation scores", tabName = "p2", icon = icon('th')),
                  menuItem("Predicted EDC scores", tabName = "p3", icon = icon('th')),
                  menuItem("Comparison with ToxPi scores", tabName = "p4", icon = icon('th'))
                 #menuItem("Prediction from MIEs", tabName = "p5", icon = icon('th'))
                )
      ),
      dashboardBody(
       tabItems(
        

### Tab 0: Home ----------------------------------------------------------------

#        tabItem(),
# To create home page with pipeline diagram and text explaining the pipeline

### Tab 1: Summary -------------------------------------------------------------

        tabItem(tabName = "dashboard",
                fluidRow(
                  box(plotOutput("plot_st1", height = 250), width = 6, height = 300),
                  box(
                    title = "Pathways used in the pipeline",
                    bsButton("qtf2", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                    bsPopover(id = "qtf2", title = "Pathways",
                              content = paste0("Only pathways related to metabolic syndrome were considered from each databases."),
                              placement = "right", 
                              trigger = "hover",
                              options = list(container = "body")),
                    plotOutput("plot_st2", height = 200), width = 6, height = 300
                    )
                  ),
                fluidRow(
                  #box(plotOutput("plot_st3",height = 400),width = 12,height = 500),
                  box(plotOutput("plot_st4", height = 400),
                      bsButton("qtf1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                      bsPopover(id = "qtf1", title = "F1 scores",
                                content = paste0('a data frame as RDS file with two columns. The first column should be named values (the F1 scores)',
                                                 '. The second column should be named networks (a similar name for all rows i.e example/outPut/F1_scores.rds).'),
                                placement = "right", 
                                trigger = "hover",
                                options = list(container = "body")),
                      fluidRow(
                        box(fileInput("F1_scores_input", label = "Add new F1-scores", accept = c("rds","A list",".rds")), width = 6, solidHeader = T),
                        box(textInput("new_F1_Colr_input", label = "Color for the new layer", value = "orange"), width = 6, solidHeader = T)
                        ),
                      width = 12,height = 580)
                  )
                ),        
        

### Tab 2: Toxicogenomics pipeline ---------------------------------------------

        tabItem(tabName = 'p1',
                headerPanel(""),
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
                        numericInput('number_cpu_input', 'number of cpus', 4, min=1, step = 1),
                        bsButton("qtjob", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                        bsPopover(id = "qtjob", title = "Job name",
                                  content = paste("Please enter a job name. This name will be used during visualization of the relults."),
                                  placement = "right", 
                                  trigger = "hover",
                                  options = list(container = "body")),
                        textInput('job_name_input', label = 'Job name', value = 'my_job'),
                        verbatimTextOutput('status_lbl'), br(),
                        height = 400, width = 6, solidHeader = T),
                    
                    box(
                      bsButton("qtnet", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                      bsPopover(id = "qtnet", title = "Network",
                                content = paste("Gene-gene co-expression network as three column (gene1, gene2, edge parameter) data frame in rds format. Maximum file size: 100mb.",
                                                 "Note: Use only Entrez IDs for gene. The edge parameter could be edge weight generated from WGCNA or wTO package.", sep = " " ),
                                placement = "right", 
                                trigger = "hover",
                                options = list(container = "body")),
                      fileInput('network_input', label = 'Gene-Gene network', accept = c("rds", "A data frame with three columns gene1 gene2 weight/p_value", ".rds")),
                      bsButton("qtedc", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                      bsPopover(id = "qtedc", title = "EDCs",
                                content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                 'to that EDC as Entrez gene IDs'),
                                placement = "right", 
                                trigger = "hover",
                                options = list(container = "body")),
                      fileInput('edc_input', label = 'List of EDCs and its MIEs',accept = c("rds","A list", ".rds")),
                      bsButton("qtdec", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                      bsPopover(id = "qtdec", title = "Negative controls",
                                content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                 'to that negative control as Entrez gene IDs'),
                                placement = "right", 
                                trigger = "hover",
                                options = list(container = "body")),
                      fileInput('decoy_input',label = 'List of negative controls and MIEs', accept = c("rds","A list",".rds")),
                      height = 400, width = 6, solidHeader = T)
                    
                  ),
                  
                  box(
                    bsButton("qtopt", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                    bsPopover(id = "qtopt", title = "Parameter selection",
                              content = paste("Selection of optimal parameters for RWR-GSEA (Optional step). The step uses pareto based solution to", 
                                              "find the optimal (read minimal) proportion of edges with highest weight from the actual network to be", 
                                              "used during RWR and the number of genes with highest visit probability to be selected post-RWR for", 
                                              "GSEA. These values are also optimized to maximize the clustering of EDCs and negative controls based on their respective MIEs.",
                                              "Here, enter all possible values to be considered during optimization.",
                                              "The combination with the highest silhouette is preferably used.", sep=" "),
                              placement = "right", 
                              trigger = "hover",
                              options = list(container = "body")),
                    fluidRow(
                      box(textInput('edge_set_input', label = "Proportion of network edges to consider for Random Walk (%):", value = '2,3,5,10'), width = 6, solidHeader = T),
                      box(textInput('gene_set_input', label = "Number of top genes to select after Random Walk :", value = '200,500,1000'), width = 6, solidHeader = T)
                    ),
                    actionButton('pareto_btn',label = 'Optimize'), 
                    tableOutput('pareto_tabl'),
                    width = 12, collapsible = T, collapsed = T, title = "Selction of optimal parameters for RWR-GSEA"),
                  
                  box(
                    fluidRow(
                      box( bsButton("qfinal_gene", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                           bsPopover(id = "qfinal_gene", title = "Final gene set for enrichment",
                                     content = paste("The number of genes with highest visit probability to be selected after RWR for gene set enrichment analysis.",
                                                     "The field is automatically updated if optimization step is used.", sep = " "),
                                     placement = "right", 
                                     trigger = "hover",
                                     options = list(container = "body")),
                           textInput('final_gene_input',label = 'Final number of genes after random walk:', value = '500'),width = 6, solidHeader = T),
                      box(bsButton("qfinal_edge", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                          bsPopover(id = "qfinal_edge", title = "Edge percent from the network",
                                    content = paste("The proportion of edges with highest weight from the actual network  to be used for RWR.",
                                                     "The field is automatically updated if optimization step is used.", sep = " "),
                                    placement = "right", 
                                    trigger = "hover",
                                    options = list(container = "body")),
                          textInput('final_edge_input',label = 'Final Edge  proportion from the network:',value = '.05'), width = 6, solidHeader = T)
                    ),
                    fluidRow(
                      box(numericInput('k_input', 'K-fold cross validation', 5, min=2, max = 20, step = 1), width = 6, solidHeader = T),
                      box(numericInput('repeat_input', 'Cross validation repeats',1, min=1, max = 10, step = 1),width = 6,solidHeader = T)
                      ),
                    actionButton('network_btn',label = 'Run RWR-FGSEA-GLM'),
                    plotOutput('f1_plt'),
                    width = 12, collapsed = T, collapsible = T, title = "Training and validation of elastic net based classifier"),
                  
                  box(
                    selectInput('export_input', 'Select data set to export as rds file:',
                                choices = c("Pathway activation scores",
                                            "Model parameters",
                                            "F1 scores",
                                            "Elastic net coefficients",
                                            "predicted_items")),
                    downloadButton('export_btn', 'Export selected item'), width = 12, title = 'Export', collapsible = T, collapsed = T),
                  
                  box(
                    fluidRow(
                      box(fileInput('test_compounds_input',label = 'A list of MIEs for unknown compounds', accept = c("rds","A list",".rds")), width = 6, solidHeader = T),
                      box(fileInput('all_model_parameters_input', label = 'Model and paramters for prediction', accept = c("rds","A list",".rds")), width = 6, solidHeader = T)),
                    actionButton('predict_btn', label = 'Predict Unknown compounds'),
                    tableOutput('prediction_tab1'), width = 12, collapsible = T, title = 'Prediction of new compounds', collapsed = T)
                 
                   #end of side bar panel2
                  , width = 12)),
        
### Tab 3: Putative pathways ---------------------------------------------------        
        
         tabItem(tabName = 'p2',
                 headerPanel("Profiling the molecular activity of EDCs / Molecular activity profiling of EDCs"), 
                 sidebarPanel(
                   bsButton("q2", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q2", title = "GLM coefficients",
                             content = paste("Select the minimum GLM coefficient for visualization"),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'GLM_coef',
                               label='GLM Coefficient cut-off',
                               min = 0, max = 1, value = 0, step = 0.01, round = F),
                   bsButton("q3", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q3", title = "Pathway activation score",
                             content = paste("Select the minimum pathway activation score for visualization"),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'NES',
                               label = "Pathway activation score cut-off",
                               min = 0, max = 1.7, value = 0, step = 0.01, round = F),
                   hr(),
                   selectInput(inputId = "pathway_category_input", 
                               label = "Pathway databases (Hold Ctrl to select multiple databases)", 
                               choices = c('KEGG','REACTOME','WIKI','GO'), 
                               multiple=TRUE, selectize=FALSE),
                   hr(),
                   selectInput(inputId = 'data_layer_input', 
                               label = 'Data layers (Hold Ctrl to select multiple data layers)', 
                               state.name, multiple = TRUE, selectize = FALSE),
                   actionButton(inputId = 'calc2', label = 'Show'),
                   hr(),
                   br(),
                   bsButton("q1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "q1", title = "New data layer",
                             content = paste0('A data frame as rds file with four columns is needed, the first column should be ', 
                                              'the name of data layer ','the second column should be the average of pathway activations scores ',
                                              'across all EDCs for that pathway , the third column is the GLM coefficients for the pathway ',
                                              'and the last column is the name of the pathway. You can select any ',
                                              'names for the columns but the order should be as above'),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('new_data_glm_input', label = 'Add new data layer', accept = c("rds", "A list", ".rds"))
                   ), #end of side bar panel2
                 mainPanel(
                   plotOutput('plt_ptway',
                              height = 900,
                              dblclick = 'ptw_dbl_click',
                              brush = brushOpts(id = 'ptw_brush',
                                                resetOnNew = T)),
                   downloadButton('export_btn_pathways', 'Export Plot data as csv'),
                   br(),
                   hr()
                   ) # end of main panel2
                 ),  # end tabpanel 2 prediction form MIEs
        
### Tab 4: Class probabilities for EDCs ----------------------------------------     

         tabItem(tabName = 'p3',         
                 headerPanel("EDC class probability"), 
                 sidebarPanel( 
                   bsButton("qtedc1", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtedc1", title = "Class probability for CTD chemicals",
                             content = paste("Enter the compound(s) for which class probabilities to visualize.",
                                               "The data layers selected below will be used for the prediction of average and harmonic sum of EDC scores.",  
                                               "Hold Ctrl key to select multiple data layers.", 
                                               "Note: Maximum 5 compounds can be viewed/compared at once.", sep = " "),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   
                   
                   selectizeInput(inputId = 'cmpname', 
                                  label = 'CAS, MESH ID, Compound name',
                                  choices = c('1962-83-0'),
                                  options = list(maxOptions = 10, maxItems = 5)),
                   selectInput(inputId = 'edc_score_layer_input', 
                               label = 'Data layers', 
                               c('PPI_STRINGdb'),
                               multiple=TRUE, selectize = FALSE),
                 
                   
                     #checkboxInput('chkbox_most_informatiave',
                   #'Select the most correlated networks with ToxCast',F),   # most informative 

                   actionButton(inputId = 'calc', label = 'Show on plot'),
                   br(),
                   hr(),
                   downloadButton('export_btn_edcscores','Export plot data as csv'),
                   # helpText("Note: You can compare the class probabilities ",
                   #          "as well as the EDC scores for maximum 5 compounds ")
                   hr(),

                   bsButton("qtfile", label = "", icon = icon("question"), style = "info", size = "extra-small"),
                   bsPopover(id = "qtfile", title = "Class probability for other chemicals",
                             content = paste("For chemicals not in CTD, its MIEs can be used as input to calculate the EDC scores.",
                                             "Enter the name of the compound and the list of MIEs as Entrez gene ID in the respective fields below.",
                                             "Select the data layers based on which to calculate the probability scores and hit the",
                                             strong("Calculate from MIEs"), "button to calculate.", br(),
                                             "Notes: 1) For this module to run, the file", 
                                             em("all_precompiled_pipeline.RDSS"), "should exist in the folder large_file.",
                                             "If missing, the", strong("Calculate from MIEs"), "button will not appear.", br(),
                                             "2) Only one compound at a time.", br(),
                                             "3) The function might return an error if the MIE is not among the data space genes.", sep = " "),
                             placement = "right", 
                             trigger = "hover",
                             options = list(container = "body")),
                   textInput(inputId = "txt_input_newcompound_name", 
                             label = "Name of the compound",
                             value = "new_compound"),

                   textInput(inputId = "txt_input_mies", 
                             label = "MIEs of the compound", 
                             value = "4617,4654,4656,5077,25937"),
                   useShinyjs(),
                   actionButton(inputId = 'mie2classprob_btn', label = 'Calculate from MIEs'),
                   hr()
                   ), # End of side bar panel 1
                 mainPanel(
                  box(tableOutput('table_edc_scores'), collapsible = T, collapsed = T, width = 12, title = "Average and harmonic sum of EDC scores"),
                  box(
                    plotOutput('plot_class_prob_scores',
                              height = 500,
                              dblclick = 'plot_class_prob_scores_dbl_click',
                              brush = brushOpts(id = 'class_prob_scores_brush', resetOnNew = T)
                              ), 
                    width = 12, collapsible = T, title = "Plot of compound class probability across data layers")
                  #plotOutput('plot_edc_score', height = 300),
                  ) # End of main panel 1
                 ),

### Tab 5: Evaluation of ToxPi -------------------------------------------------   


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





#Updated by: Arindam Ghosh

