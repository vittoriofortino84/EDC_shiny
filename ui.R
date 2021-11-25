#ui.R

### Libraries to import --------------------------------------------------------

library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
library(DT)

### Sidebar menu of the app ---------------------------------------------------------------
 
shinyUI(
  dashboardPage(
    dashboardHeader(title = "EDTox"),
    dashboardSidebar(
      sidebarMenu(id = "tabs",
                  menuItem("Home", tabName = "tab_home", icon = icon("dashboard")),
                  menuItem("Summary", tabName = "tab_dashboard", icon = icon("th")),
                  menuItem("Toxicogenomics pipeline", tabName = "tab_pipeline", icon = icon('th')),
                  menuItem("Pathway activation scores", tabName = "tab_pathway", icon = icon('th')),
                  menuItem("EDC-class probability scores", tabName = "tab_edcScore", icon = icon('th')),
                  menuItem("Comparison with ToxPi scores", tabName = "tab_toxpiScore", icon = icon('th'))
                )
      ),
      dashboardBody(
       tags$script(HTML("$('body').addClass('fixed');")), #To fix header and sidebar
       tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }'))), #To keep fixed background (esp in tab2)
             fluidPage(
      									
																  
    
 
    ),
		
		
		
       tabItems(
        

### Tab 1: Home ----------------------------------------------------------------
# Home page with pipeline diagram and text explaining the pipeline
          tabItem(tabName = "tab_home",
		                    fluidRow(column(width = 12, includeHTML("www/Home.html"))),
				              ),

### Tab 2: Summary -------------------------------------------------------------

         tabItem(tabName = "tab_dashboard",
                 fluidRow(
                   box(title = strong("Distribution of EDC-class probability scores of CTD chemicals"),
                     plotOutput("plot_st1", height = 250), width = 6, height = 315),
                   box(
                     title = list(strong("Pathways used in the pipeline"), bsButton("qt_path", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                     bsPopover(id = "qt_path", title = "Pathways",
                               content = paste0("Only gene sets related to adverse outcome pathways were retrieved from GO."),
                               placement = "right",
                               trigger = "hover",
                               options = list(container = "body")),
                     plotOutput("plot_st2", height = 250), width = 6, height = 315
                     )
                   ),
                 fluidRow(
                   box(title = strong("Accuracy of trained classifiers"),
                       plotOutput("plot_st4", height = 400),
                       bsPopover(id = "qt_f1Score", title = "F1 scores",
                                 content = paste("A data frame as RDS file with two columns. The first column should be named values (the F1 scores).",
                                                  "The second column should be named networks (a similar name for all rows)",
                                                 "The new F1 scores would be appended at the end of the plot.", sep = " "),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fluidRow(
                         box(
                           fileInput("F1_scores_input", 
                                     label = list("Add classification performance of new classifier", bsButton("qt_f1Score", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                                     accept = c("rds","A list",".rds")), 
                           width = 6, solidHeader = T),
                         #box(textInput("new_F1_Colr_input", label = "Color for the new layer", value = "orange"), width = 6, solidHeader = T)
                         ),
                       width = 12, height = 580)
                   )
                 ),


### Tab 3: Toxicogenomics pipeline ---------------------------------------------

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
                         numericInput('number_cpu_input', 'Number of CPUs', 1, min=1, step = 1),
                         bsPopover(id = "qt_jobName", title = "Classifier name",
                                   content = paste("Please enter a classifier name. This name will be used during visualization of the relults."),
                                   placement = "right",
                                   trigger = "hover",
                                   options = list(container = "body")),
                         textInput('job_name_input',
                                   label = list("Classifier name", bsButton("qt_jobName", label = "", icon = icon("question"), style = "info", size = "extra-small")),  
                                   value = 'my_job'),
                         verbatimTextOutput('status_lbl'), br(),
                         height = 400, width = 6, solidHeader = T),

                     box(
                       bsPopover(id = "qt_net", title = "Network",
                                 content = paste("Gene-gene co-expression network as three column (gene1, gene2, edge parameter) data frame in rds format. Maximum file size: 50mb.",
                                                  br(), "Note: Use only Entrez IDs for gene. The edge parameter could be edge weight generated from WGCNA or wTO package. To download example: refer to Home page", sep = " " ),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('network_input',
                                 label = list("Gene-Gene network", bsButton("qt_net", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                                 accept = c("rds", "A data frame with three columns gene1 gene2 weight/p_value", ".rds")),
                       bsPopover(id = "qt_edc", title = "EDCs",
                                 content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                  'to that EDC as Entrez gene IDs'),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('edc_input',
                                 label = list("List of EDCs and its MIEs", bsButton("qt_edc", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                                 accept = c("rds","A list", ".rds")),
                       bsPopover(id = "qt_decoy", title = "Negative controls",
                                 content = paste0('Named R list as rds file where each item contain the MIEs releated  ',
                                                  'to that negative control as Entrez gene IDs'),
                                 placement = "right",
                                 trigger = "hover",
                                 options = list(container = "body")),
                       fileInput('decoy_input',
                                 label = list("List of negative controls and MIEs", bsButton("qt_decoy", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
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
                     title = list("Selection of optimal parameters for RWR-FGSEA", bsButton("qt_optimize", label = "", icon = icon("question"), style = "info", size = "extra-small"))),

                   box(
                     fluidRow(
                       box(
                         bsPopover(id = "qt_rwrEdge", title = "Edge proportion from the network",
                                     content = paste("The proportion of edges with highest weight from the actual network  to be used for RWR.",
                                                      "The field is automatically updated if optimization step is used.", sep = " "),
                                     placement = "right",
                                     trigger = "hover",
                                     options = list(container = "body")),
                           textInput('final_edge_input',
                                     label = list("Proportion of network edges to consider for Random Walk (%):", bsButton("qt_rwrEdge", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                                     value = '5'), 
                         width = 6, solidHeader = T),                       
                       box( 
                           bsPopover(id = "qt_rwrGene", title = "Final gene set for enrichment",
                                     content = paste("The number of genes with highest visit probability to be selected after RWR for gene set enrichment analysis.",
                                                     "The field is automatically updated if optimization step is used.", sep = " "),
                                     placement = "right",
                                     trigger = "hover",
                                     options = list(container = "body")),
                           textInput('final_gene_input',
                                     label = list("Number of top genes to select after Random Walk:", bsButton("qt_rwrGene", label = "", icon = icon("question"), style = "info", size = "extra-small")),  
                                     value = '500'),
                           width = 6, solidHeader = T)
                       ),
                     fluidRow(
                       box(numericInput('k_input', 'K-fold cross validation', 5, min=2, max = 20, step = 1), width = 6, solidHeader = T),
                       box(numericInput('repeat_input', 'Cross validation repeats',1, min=1, max = 10, step = 1),width = 6,solidHeader = T)
                       ),
                     actionButton('network_btn',label = 'Run RWR-FGSEA-GLM'),
                     plotOutput('f1_plt'),
                     width = 12, collapsed = T, collapsible = T, title = "Training and validation of GLM based classifier"),

                   
                   bsPopover(id = "qt_export", title = "Export options",
                             content = paste(strong("model_parameters:"), "For parameters of the trained model", br(),
                                             strong("cls_perfs:"), "For F1 scores of the trained model", br(),
                                             strong("edc_moas:"), "For GLM elastic-net coefficients and pathway activation scores", br(),
                                             strong("predicted_edcscores:"), "For EDC-class probabilities predicted using the newly trained classifier", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   
                   box(
                     selectInput(inputId = 'export_input', 
                                 label = list("Select data set to export as rds file:", bsButton("qt_export", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                                 choices = c(
					     #"Pathway activation scores",
                                             "model_parameters",
                                             "cls_perfs",
                                             "edc_moas",
                                             "predicted_edcscores")),
                     
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

### Tab 4: Molecular activity profiling of EDCs --------------------------------

         tabItem(tabName = "tab_pathway",
                 headerPanel("Molecular activity profiling of EDCs"),
                 sidebarPanel(
                   bsPopover(id = "qt_glmCoef", title = "GLM coefficients",
                             content = paste("Select the minimum GLM coefficient"),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'GLM_coef',
                               label = list("GLM coefficient threshold ", bsButton("qt_glmCoef", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                               min = 0, max = 1, value = 0, step = 0.01, round = F),
                   bsPopover(id = "qt_pathScore", title = "Pathway activation score",
                             content = paste("Select the minimum pathway activation score"),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   sliderInput(inputId = 'NES',
                               label = list("Pathway activation score threshold ", bsButton("qt_pathScore", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                               min = 0, max = 2, value = 0, step = 0.01, round = F),
                   hr(),
                   selectInput(inputId = "pathway_category_input",
                               label = "Pathway databases (Hold Ctrl to select multiple databases)",
                               choices = c('KEGG','REACTOME','WIKI','GO'),selected=c('KEGG','REACTOME','WIKI','GO'),
                               multiple=TRUE, selectize=FALSE),
                   hr(),
                   selectInput(inputId = 'data_layer_input',
                               label = "Classifiers (Hold Ctrl select multiple classifiers)",
                               state.name, multiple = TRUE, selectize = FALSE),
                   actionButton(inputId = 'calc2', label = 'Update Plot'),
                   hr(),
                   br(),
                   bsPopover(id = "qt_newLayer", title = "Pathway activation scores from new classifier",
                             content = paste("Pathway activation scores from a new classifier can be imported as a four column data frame in rds format.",
                                               "The columns should list the name of the classifier, the average pathway activation scores across all EDCs for the pathway,",
                                               "the GLM coefficient for the pathway and the name of the pathway, respectively in the order mentioned.", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   fileInput('new_data_glm_input', 
                             label = list("Add new classifier", bsButton("qt_newLayer", label = "", icon = icon("question"), style = "info", size = "extra-small")),  
                             accept = c("rds", "A list", ".rds"))
                   ), #end of side bar panel2
                 mainPanel(
                   plotOutput('plt_ptway',
                              height = 900,
                              dblclick = 'ptw_dbl_click',
                              brush = brushOpts(id = 'ptw_brush',
                                                resetOnNew = T)),
                   downloadButton('export_btn_pathways', 'Export plot data'),
                   br(),
                   hr()
                   ) # end of main panel2
                 ),  # end tabpanel 2 prediction form MIEs

### Tab 5: EDC class probability scores-----------------------------------------------

         tabItem(tabName = "tab_edcScore",
                 headerPanel("EDC-class probability scores"),
                 sidebarPanel(
                   p(strong("Compile EDC-class probability scores for CTD chemicals "), bsButton("qt_scoreCTD", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                   bsPopover(id = "qt_scoreCTD", title = "Compile EDC-class probability scores for chemicals annotated in Comparative Toxicogenomics Database (CTD)",
                             content = paste("Enter the compound(s) for which class probabilities to visualize.",
                                               "The classifiers selected below will be used to generate the EDC-class probability scores.",
                                               "Hold", strong("Ctrl"), "key to select multiple data layers.", br(),
                                               "Note: Maximum 5 compounds can be viewed/compared at once.", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   box(

                   selectizeInput(inputId = 'addtocmpname',
                                  label = 'Library of compound names',
                                  choices = c(''),
                                  options = list(maxOptions = 10, maxItems = 5)), 
                   actionButton(inputId = 'comp_dic_btn',
                                label = 'Fetch compound list'),collapsible = T,collapsed = T,width = 14,title="Compound Search",solidHeader=T),
hr(),

			       textInput('cmpname',
			                  label = 'CAS, MeSH ID, Compound name',value = 'Bisphenol B,C030298,1962-83-0'),



                   selectInput(inputId = 'edc_score_layer_input',
                               label = 'Classifiers',
                               c('PPI_STRINGdb'),
                               multiple=TRUE, selectize = FALSE),


                  actionButton(inputId = 'calc',
                                label = 'Show on plot'),br(),hr(),
                  bsPopover(id = "qallcompile", title = "Scores for all CTD chemicals",
                            content = paste0(" Calculate EDC-class probability scores for all the chemicals in CTD based on all the classifiers.",
                                             "The classifiers selected above will be used to compute the average and harmonic scores for each compound.",
                                             "Use the", strong("Export"), "button below to export the results as csv file.",
                                             "Note: Compilation of the scores might take time", sep = " "),
                            placement = "right", 
                            trigger = "hover",
                            options = list(container = "body")),																													 
                            br(),hr(),
                  downloadButton('export_all_btn_edcscores', 'Export scores for all CTD chemicals'),
			            checkboxInput('harmonicExport', 'include Harmonic EDC score', value = F, width = NULL),
                       hr(),
              
                   p(strong("Compile EDC-class probability scores for new chemicals "), bsButton("qt_scoreOther", label = "", icon = icon("question"), style = "info", size = "extra-small")),
                   bsPopover(id = "qt_scoreOther", title = "Compile EDC-class probability scores for new chemicals",
                             content = paste("Enter the name of the compound and the comma separated list of MIEs as Entrez gene ID or Symbols in the respective fields below.i.e", strong("4617,4654,ERK3,4656,5077,25937"),
                                             "Select the classifiers based on which to calculate the probability scores and hit the",
                                             strong("Calculate from MIEs"), "button to calculate.", br(),
                                             "", br(),
                                             "", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   textInput(inputId = "txt_input_newcompound_name",
                             label = "Name of the new compound",
                             value = "new_compound"),

                   textInput(inputId = "txt_input_mies",
                             label = "MIEs of the new compound",
                             value = ""),
                   useShinyjs(),
                   actionButton(inputId = 'mie2classprob_btn', label = 'Calculate from MIEs'),
                   hr()
                   ), # End of side bar panel 1
                 mainPanel(
                  box(tableOutput('table_edc_scores'), collapsible = T, collapsed = T, width = 12, title = "Average and harmonic sum of EDC-class probability scores"),
                  box(
                    plotOutput('plot_class_prob_scores',
                              height = 500,
                              dblclick = 'plot_class_prob_scores_dbl_click',
                              brush = brushOpts(id = 'class_prob_scores_brush', resetOnNew = T)
                              ),
                    width = 12, collapsible = T, title = "Plot of compound EDC-class probability scores for selected classifiers"),
					      
                       downloadButton('export_btn_edcscores','Export plot data')
                  #plotOutput('plot_edc_score', height = 300),
                  ) # End of main panel 1
                 ),

### Tab 6: Comparison with ToxPi Scores ----------------------------------------

         tabItem(tabName = "tab_toxpiScore",
                 headerPanel("Comparison with ED-based ToxPi Scores"),
                 sidebarPanel(
                   bsPopover(id = "qt_dataLayer", title = "Classifiers",
                             content = paste("Select the classifiers to be used for calculation of average EDC-class probability scores.",
                                                "Hold", strong("Ctrl"), "to select multiple data layers.", sep = " "),
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")),
                   selectInput(inputId = 'toxpi_layer_input', 
                               label = list("Classifiers", bsButton("qt_dataLayer", label = "", icon = icon("question"), style = "info", size = "extra-small")), 
                               c('PPI_STRINGdb'), multiple = T, selectize = F),
                   sliderInput(inputId = 'slider_toxpi_plt_toxpi_cutoff',
                               label='Minimum ToxPi score:',
                               min = 0.000, max = 0.5, value = 0.1, step = 0.01, round = F),
                   sliderInput(inputId = 'slider_toxpi_plt_edc_score_cutoff',
                               label ='Minimum EDC-class probability score:',
                               min = 0.000, max = 1, value = 0.8, step = 0.01, round = F),
                   sliderInput(inputId = 'slider_click_toxpi',
                               label = 'Plot click prescision:',
                               min = 0.000, max = 0.05, value = 0.005, step = 0.0001, round = F),
                   actionButton(inputId = 'toxpi_btn', label = 'Calculate'),
                   br(),
                   hr(),

                   downloadButton('export_btn_toxpi', 'Export plot data')
                   #verbatimTextOutput("txt_toxpi_selected_layers",placeholder = F)
                   ),
                 mainPanel(
                   box(plotOutput('plot_toxpi', click = "toxpi_plot_click"
                                  #dblclick = 'toxpi_plot_dbl_click',
                                  #brush = brushOpts(id = 'toxpi_brush',
                                  #resetOnNew = T)
                                  ), collapsible = T, width = 12, title = "Average EDC-class probabilities VS ToxPi scores"),
                   verbatimTextOutput("txt_toxpi_click"),
                   #verbatimTextOutput("txt_toxpi_compounds"),
                   #tags$head(tags$style("#txt_toxpi_compounds{color:blue; font-size:14px; font-style:bold; overflow-y:scroll; max-height: 100px; background: ghostwhite;}")),
                   dataTableOutput('table_toxpi')
                   ) # End of main panel 1
                 ) # end of tab item toxcast


### Tab #: Any new tabs here

      ) # End of all tabsets panel
    )
  )
)
