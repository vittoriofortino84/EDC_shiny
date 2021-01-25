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

```
#### Launch locally
## Clone the git repository
#### git clone https://github.com/vittoriofortino84/EDC_shiny EDC_shiny_clone
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
