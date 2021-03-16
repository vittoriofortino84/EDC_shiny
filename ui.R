#ui.R

### Libraries to import --------------------------------------------------------

library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
library(DT)

### Sidebar menu ---------------------------------------------------------------

shinyUI(
  dashboardPage(
    dashboardHeader(title = "EDTox:GUI"),
    dashboardSidebar(
      sidebarMenu(id = "tabs",
                  menuItem("Home", tabName = "tab_home", icon = icon("dashboard")),
                  menuItem("Summary", tabName = "tab_dashboard", icon = icon("th")),
                  menuItem("Toxicogenomics Pipeline", tabName = "tab_pipeline", icon = icon('th')),
                  menuItem("Pathway activation scores", tabName = "tab_pathway", icon = icon('th')),
                  menuItem("EDC-class probability", tabName = "tab_edcScore", icon = icon('th')),
                  menuItem("Comparison with ToxPi scores", tabName = "tab_toxpiScore", icon = icon('th'))
                  #menuItem("Prediction from MIEs", tabName = "p5", icon = icon('th'))
                )
      ),
      dashboardBody(
       tags$script(HTML("$('body').addClass('fixed');")), #To fix header and sidebar
       tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))), #To keep fixed background (esp in tab2)
        
       tabItems(
        

### Tab 0: Home ----------------------------------------------------------------
# Home page with pipeline diagram and text explaining the pipeline

         tabItem(tabName = "tab_home",
                 fluidRow(box(img(style = "max-width: 100%; width: 70%; height: auto", src = "images/EDTox_Pipeline.jpg"),
                          p(style = "font-size: 125%; text-align: justify; margin: 0.5% 1% 0.5% 1%", "The EDTox:GUI is an R Shiny appliation that provides an interactive graphical user interface to utilize the EDTox dataspace", 
                              "for prediction of endocrine disruption potential of a chemical compound and prediction of possible adverse outcome pathways (AOP).",
                              "The R Shiny application serves two purposes. First, it can be used to generate models for classification of chemicals",
                              "as EDCs or non-EDCs based on networks derieved from toxicogenomics data.",
                              "And the second, it can be used for infering endocrine disruption potential based on the EDTox dataspace"), width = 20), 
                          align = "center")
                 ),

### Tab 1: Summary -------------------------------------------------------------

         tabItem(tabName = "tab_dashboard",
                 fluidRow(
                   box(plotOutput("plot_st1", height = 250), width = 6, height = 315),
                   box(
                     title = list(strong("Pathways used in the pipeline"), bsButton("qt_path", label = "", icon = icon("question"), style = "info", size = "extra-small")),#---
                     bsPopover(id = "qt_path", title = "Pathways",
                               content = paste0("Only pathways related to metabolic syndrome were considered from each databases."),
                               placement = "right",
                               trigger = "hover",
                               options = list(container = "body")),
                     plotOutput("plot_st2", height = 250), width = 6, height = 315
                     )
                   ),
                 fluidRow(
                   #box(plotOutput("plot_st3",height = 400),width = 12,height = 500),
                   box(plotOutput("plot_st4", height = 400),
                       bsPopover(id = "qt_f1Score", title = "F1 scores",
                                 content = paste0('A data frame as RDS file with two columns. The first column should be named values (the F1 scores)',
                                                  '. The second column should be named networks (a similar name for all rows i.e example/outPut/F1_scores.rds).'),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fluidRow(
                         box(
                           fileInput("F1_scores_input", 
                                     label = list("Add new F1-scores", bsButton("qt_f1Score", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                                     accept = c("rds","A list",".rds")), 
                           width = 6, solidHeader = T),
                         box(textInput("new_F1_Colr_input", label = "Color for the new layer", value = "orange"), width = 6, solidHeader = T)
                         ),
                       width = 12, height = 580)
                   )
                 ),


### Tab 2: Toxicogenomics pipeline ---------------------------------------------

         tabItem(tabName = "tab_pipeline",
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
                         bsPopover(id = "qt_jobName", title = "Job name",
                                   content = paste("Please enter a job name. This name will be used during visualization of the relults."),
                                   placement = "right",
                                   trigger = "hover",
                                   options = list(container = "body")),
                         textInput('job_name_input',
                                   label = list("Job name", bsButton("qt_jobName", label = "", icon = icon("question"), style = "info", size = "extra-small")),  #---
                                   value = 'my_job'),
                         verbatimTextOutput('status_lbl'), br(),
                         height = 400, width = 6, solidHeader = T),

                     box(
                       bsPopover(id = "qt_net", title = "Network",
                                 content = paste("Gene-gene co-expression network as three column (gene1, gene2, edge parameter) data frame in rds format. Maximum file size: 100mb.",
                                                  br(), "Note: Use only Entrez IDs for gene. The edge parameter could be edge weight generated from WGCNA or wTO package.", sep = " " ),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('network_input',
                                 label = list("Gene-Gene network", bsButton("qt_net", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                                 accept = c("rds", "A data frame with three columns gene1 gene2 weight/p_value", ".rds")),
                       bsPopover(id = "qt_edc", title = "EDCs",
                                 content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                  'to that EDC as Entrez gene IDs'),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('edc_input',
                                 label = list("List of EDCs and its MIEs", bsButton("qt_edc", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                                 accept = c("rds","A list", ".rds")),
                       bsPopover(id = "qt_decoy", title = "Negative controls",
                                 content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                  'to that negative control as Entrez gene IDs'),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('decoy_input',
                                 label = list("List of negative controls and MIEs", bsButton("qt_decoy", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                                 accept = c("rds","A list",".rds")),
                       height = 400, width = 6, solidHeader = T)

                   ),

                   box(
                     bsPopover(id = "qt_optimize", title = "Parameter selection",
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
                       box(textInput('gene_set_input', label = "Number of top genes to select after Random Walk:", value = '200,500,1000'), width = 6, solidHeader = T)
                     ),
                     actionButton('pareto_btn',label = 'Optimize'),
                     tableOutput('pareto_tabl'),
                     width = 12, collapsible = T, collapsed = T, 
                     title = list("Selction of optimal parameters for RWR-FGSEA", bsButton("qt_optimize", label = "", icon = icon("question"), style = "info", size = "extra-small"))),#---

                   box(
                     fluidRow(
                       box( 
                         bsPopover(id = "qt_rwrGene", title = "Final gene set for enrichment",
                                   content = paste("The number of genes with highest visit probability to be selected after RWR for gene set enrichment analysis.",
                                                    "The field is automatically updated if optimization step is used.", sep = " "),
                                   placement = "right",
                                   trigger = "hover",
                                   options = list(container = "body")),
                         textInput('final_gene_input',
                                   label = list("Number of top genes to select after Random Walk:", bsButton("qt_rwrGene", label = "", icon = icon("question"), style = "info", size = "extra-small")), #--- 
                                   value = '500'),
                         width = 6, solidHeader = T),
                       box(
                         bsPopover(id = "qt_rwrEdge", title = "Edge proportion from the network",
                                     content = paste("The proportion of edges with highest weight from the actual network  to be used for RWR.",
                                                      "The field is automatically updated if optimization step is used.", sep = " "),
                                     placement = "right",
                                     trigger = "hover",
                                     options = list(container = "body")),
                           textInput('final_edge_input',
                                     label = list("Proportion of network edges to consider for Random Walk (%):", bsButton("qt_rwrEdge", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                                     value = '5'), 
                         width = 6, solidHeader = T)
                     ),
                     fluidRow(
                       box(numericInput('k_input', 'K-fold cross validation', 5, min=2, max = 20, step = 1), width = 6, solidHeader = T),
                       box(numericInput('repeat_input', 'Cross validation repeats',1, min=1, max = 10, step = 1),width = 6,solidHeader = T)
                       ),
                     actionButton('network_btn',label = 'Run RWR-FGSEA-GLM'),
                     plotOutput('f1_plt'),
                     width = 12, collapsed = T, collapsible = T, title = "Training and validation of elastic net based classifier"),

                   box(
                     selectInput(inputId = 'export_input', 
                                 label = "Select data set to export as rds file:",
                                 choices = c("Pathway activation scores",
                                             "Model parameters",
                                             "F1 scores",
                                             "Elastic net coefficients",
                                             "Predicted scores")),
                     
                     downloadButton('export_btn', 'Export selected item'), 
                     width = 12, title = 'Export', collapsible = T, collapsed = T),

                   box(
                     fluidRow(
                       box(fileInput('test_compounds_input',label = 'A list of MIEs for unknown compounds', accept = c("rds","A list",".rds")), width = 6, solidHeader = T),
                       box(fileInput('all_model_parameters_input', label = 'Model and paramters for prediction', accept = c("rds","A list",".rds")), width = 6, solidHeader = T)),
                     actionButton('predict_btn', label = 'Predict unknown compounds'),
                     tableOutput('prediction_tab1'), width = 12, collapsible = T, title = 'Prediction of new compounds', collapsed = T)

                    #end of side bar panel2
                   , width = 12)),

### Tab 3: Molecular activity profiling of EDCs --------------------------------

         tabItem(tabName = "tab_pathway",
                 headerPanel("Molecular activity profiling of EDCs"),
                 sidebarPanel(
                   bsPopover(id = "qt_glmCoef", title = "GLM coefficients",
                             content = paste("Select the minimum GLM coefficient for visualization"),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'GLM_coef',
                               label = list("GLM Coefficient cut-off ", bsButton("qt_glmCoef", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                               min = 0, max = 1, value = 0, step = 0.01, round = F),
                   bsPopover(id = "qt_pathScore", title = "Pathway activation score",
                             content = paste("Select the minimum pathway activation score for visualization"),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'NES',
                               label = list("Pathway activation score cut-off", bsButton("qt_pathScore", label = "", icon = icon("question"), style = "info", size = "extra-small")), #---
                               min = 0, max = 1.7, value = 0, step = 0.01, round = F),
                   hr(),
                   selectInput(inputId = "pathway_category_input",
                               label = "Pathway databases (Hold Ctrl to select multiple databases)",
                               choices = c('KEGG','REACTOME','WIKI','GO'),
                               multiple=TRUE, selectize=FALSE),
                   hr(),
                   selectInput(inputId = 'data_layer_input',
                               label = "Data layers (Hold Ctrl select multiple data layers)",
                               state.name, multiple = TRUE, selectize = FALSE),
                   actionButton(inputId = 'calc2', label = 'Show'),
                   hr(),
                   br(),
                   bsPopover(id = "qt_newLayer", title = "New data layer(s)",
                             content = paste("New data layer can be imported as a four column data frame in rds format.",
                                               "The columns should list the name of the data layer, the average pathway activation scores across all EDCs for the pathway,",
                                               "the GLM coefficient for the pathway and the name of the pathway, respectively in the order mentioned.", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('new_data_glm_input', 
                             label = list("Add new data layer", bsButton("qt_newLayer", label = "", icon = icon("question"), style = "info", size = "extra-small")),  #---
                             accept = c("rds", "A list", ".rds"))
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

### Tab 4: EDC class probability -----------------------------------------------

         tabItem(tabName = "tab_edcScore",
                 headerPanel("EDC-class probability"),
                 sidebarPanel(
                   p(strong("Class probability for CTD chemicals "), bsButton("qt_scoreCTD", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                   bsPopover(id = "qt_scoreCTD", title = "Class probability for CTD chemicals",
                             content = paste("Enter the compound(s) for which class probabilities to visualize.",
                                               "The data layers selected below will be used for the calculation of average and harmonic sum of EDC-class probabilities.",
                                               "Hold", strong("Ctrl"), "key to select multiple data layers.", br(),
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
                   p(strong("Class probability for other chemicals "), bsButton("qt_scoreOther", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                   bsPopover(id = "qt_scoreOther", title = "Class probability for other chemicals",
                             content = paste("For chemicals not in CTD, its MIEs can be used as input to calculate the EDC-class probabilities.",
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
                  box(tableOutput('table_edc_scores'), collapsible = T, collapsed = T, width = 12, title = "Average and harmonic sum of EDC-class probabilities"),
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

### Tab 5: Comparison with ToxPi Scores ----------------------------------------

         tabItem(tabName = "tab_toxpiScore",
                 headerPanel("Comparison with ToxPi Scores"),
                 sidebarPanel(
                   bsPopover(id = "qt_dataLayer", title = "Data layers",
                             content = paste("Select data layers to be used for calculation of average EDC scores.",
                                                "Hold", strong("Ctrl"), "to select multiple data layers.", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   #---selectInput(inputId = 'toxpi_layer_input', label = 'Data layers', c('PPI_STRINGdb'), multiple = T, selectize = F),
                   selectInput(inputId = 'toxpi_layer_input', 
                               label = list("Data layers", bsButton("qt_dataLayer", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                               c('PPI_STRINGdb'), multiple = T, selectize = F),
                   sliderInput(inputId = 'slider_toxpi_plt_toxpi_cutoff',
                               label='Minimum ToxPi score:',
                               min = 0.000, max = 0.4, value = 0.1, step = 0.01, round = F),
                   sliderInput(inputId = 'slider_toxpi_plt_edc_score_cutoff',
                               label ='Minimum EDC score:',
                               min = 0.000, max = 1, value = 0.8, step = 0.01, round = F),
                   sliderInput(inputId = 'slider_click_toxpi',
                               label = 'Plot click prescision:',
                               min = 0.000, max = 0.05, value = 0.005, step = 0.0001, round = F),
                   actionButton(inputId = 'toxpi_btn', label = 'Calculate'),
                   br(),
                   hr(),
                   actionButton(inputId = 'toxpiBtn_refresh', label = 'Refresh'),
                   br(),
                   hr(),
                   downloadButton('export_btn_toxpi', 'Export plot data as csv')
                   #verbatimTextOutput("txt_toxpi_selected_layers",placeholder = F)
                   ),
                 mainPanel(
                   box(plotOutput('plot_toxpi',click = "toxpi_plot_click"
                                  #dblclick = 'toxpi_plot_dbl_click',
                                  #brush = brushOpts(id = 'toxpi_brush',
                                  #resetOnNew = T)
                                  ) ,collapsible = T, width = 12),
                   verbatimTextOutput("txt_toxpi_click"),
                   #verbatimTextOutput("txt_toxpi_compounds"),
                   #tags$head(tags$style("#txt_toxpi_compounds{color:blue; font-size:14px; font-style:bold; overflow-y:scroll; max-height: 100px; background: ghostwhite;}")),
                   dataTableOutput('table_toxpi')
                   ) # End of main panel 1
                 ), # end of tab item toxcast

### Tab 6: mie 2 class prob -------------------------------------------------
         tabItem(tabName = 'p5',
                 headerPanel('Predcition of class probabilities from MIEs'),
                 sidebarPanel(),
                 mainPanel()
                 ) # End of prediction from mIEs

### Tab #: Any new tabs here

      ) # End of all tabsets panel
    )
  )
)





#Updated by: Arindam Ghosh
