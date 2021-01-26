# EDC_shiny
## Run EDC_shiny
## Install Dependencies
#### Install CRAN dependencies

```r
cran_pkgs <- c("doParallel", "foreach", "igraph", "plyr", "shiny","shinyjs", "shinyBS", "shinydashboard", "DT", "ggplot2","RColorBrewer", "dnet","tidyr","ggpubr","rPref","caret","glmnet","shape")

cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]

if(length(cran_pkgs.inst)>0){ 
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."));  
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE);print("Installed!!!") 
  }
}
```

#### Install Bioconductor dependencies for R versions < 3.5
```r
source("http://bioconductor.org/biocLite.R")
bioc_pkgs <- c( "supraHex","hexbin",  "fgsea","Rgraphviz")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))];

if(length(bioc_pkgs.inst)>0){
  source("http://bioconductor.org/biocLite.R")
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"));  
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."));  
    biocLite(pkg, suppressUpdates=TRUE); print("Installed!!!")
  }
}
```
#### Install Bioconductor dependencies for R versions > 3.5
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
bioc_pkgs <- c( "supraHex","hexbin",  "fgsea","Rgraphviz")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))];
if(length(bioc_pkgs.inst)>0){
 print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"));  
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."));  
    BiocManager::install(pkg)
    print("Installed!!!")
  }
}
```
## Launch from GitHub

```r
runGitHub("EDC_shiny", "vittoriofortino84", subdir="EDC_shiny")
## Using the archived file
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.tar.gz", subdir="EDC-shiny")
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.zip", subdir="EDC-shiny")
```
## Launch locally
```
git clone https://github.com/vittoriofortino84/EDC_shiny EDC_shiny_clone
```

## Launch from R

```r
## Run by using runApp()
setwd("~/EDC-shiny")
shiny::runApp(appDir = getwd(), port = getOption("shiny.port"),
              launch.browser = getOption("shiny.launch.browser", interactive()),
              host = getOption("shiny.host", "127.0.0.1"), workerId = "",
              quiet = FALSE, display.mode = c("auto", "normal", "showcase"), 
              test.mode = getOption("shiny.testmode", FALSE))
```

## EDC_shiny Execution Steps
### Summary Tab
![alt text](https://github.com/vittoriofortino84/EDC_shiny/blob/master/example/figure/Summary_tab.png)


- presents an overview of the 24 toxicogenomics data layers, the pathways used in this application as well as the distribution of the predicted toxicity scores across all 24 data layers for 12 k compounds in CTD.
- The F1 accuracy of different data layers are compared using one boxplot. 
- The input file bellow the F1 scores plot permits the user to upload and compare F1 accuracy score from a new generated machine learning based classifier. 
### Toxicogenomics Pipeline Tab
![alt text](https://github.com/vittoriofortino84/EDC_shiny/blob/master/example/figure/Toxicogenomics_pipeline.png)

- In this tab the concept of gene-gene network, random walk with restart and gene set enrichment analysis will be used for classification of EDCs and negative controls (decoys). Meanwhile, it is possible to predict the EDC class probability for new compounds using the generated classifier model.
The user can select the CPU usage for the shiny application in the box named number of CPUs from the main panel. A proper name for the Job name box is highly recommended to make it possible for further analyses of the results including comparison of F1-scores, GLM coefficients and pathway activation scores between previously generated data layers and the one generated with this application.
Gene-gene networks should be one RDS file with three columns; the first two columns represent the edges as the nodes related to gene co-expression network and the third column is the weights or topological overlaps between the edges. Different packages including wTO or wGCNA can be utilized to prepare the gene-gene networks for EDC-shiny application.
For the next two boxes including list of EDCs and MIEs and list of negative controls and MIEs, the user should provide the list of compounds and their MIEs as RDS files for the compounds used as training set (EDCs and negative controls).

![alt text](https://github.com/vittoriofortino84/EDC_shiny/blob/master/example/figure/toxicogenomics_pipeline2.png)
- To increase the homogeneity of the networks and enhance robustness of enrichment analysis, the user can select to use pareto solution for optimization. In this step the pipeline will be optimized based on three parameters: (1) selection of the top K (.02, .03, .05, .1) proportion of the edges with respect to their topological overlap or correlation; (2) selection of the top N (200, 500, 1000) genes (based on the probabilities calculated by the RWR) with the highest proximities with respect to ED-MIEs; (3) silhouette scores to measure the EDC-cluster cohesion based on the set of N genes selected with the RWR. The final combination of K and N values giving the highest silhouette score within each network will be suggested by the EDC shiny application for the further random walk with restart and gene set enrichment steps of the pipeline.

![alt text](https://github.com/vittoriofortino84/EDC_shiny/blob/master/example/figure/toxicogenomics_pipeline3.png)
- In RWR-GLM panel, the two parameters related to number of the selected genes after random walk and the proportion of the edges from the network will be automatically updated after termination of pareto solution. However, the user can manipulate these parameters as needed. By Running RWR-FGSEA-GLM, the MIEs for EDCs and negative controls will participate in random walk with restart simulation experiments using the optimized networks and the resulting sorted vector of genes will be used in fast gene set enrichment analysis. EDC-shiny acts as the interface to dnet, FGSEA and Caret packages to do all the experiments. The pathway activation scores from the gene set enrichment analysis will be used by elastic-net generalized linear approach to setup models for classification of EDCs and negative controls. The accuracy of the obtained machine learning based classifier will be evaluated by F1-scores using k-fold cross validation. The package doParallel is used for speeding up the cross validation experiments and the user can select the number of k folds and the cross-validation repeats (cross validation repeat > 1 is not recommended on common desktop or laptop computers).


