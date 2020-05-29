# EDC_shiny
EDC_shiny Execution Steps
Run EDC_shiny
Install Dependencies
#Install CRAN dependencies
cran_pkgs <- c(  "doParallel", "foreach", "igraph", "plyr", "shiny",
                 "shinyjs", "shinyBS", "shinydashboard",  "DT", 
                  "ggplot2","RColorBrewer", "dnet","tidyr")

cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
    print("Installed!!!")
  }
}
#Install Bioconductor dependencies

# if(!"GOSemSim" %in% rownames(installed.packages())){
#   print("Installing GOSemSim from GitHub!")
#   devtools::install_github("GuangchuangYu/GOSemSim")
# }
source("http://bioconductor.org/biocLite.R")
bioc_pkgs <- c( "supraHex","hexbin",  "fgsea")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
  source("http://bioconductor.org/biocLite.R")
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."))
    biocLite(pkg, suppressUpdates=TRUE)
    print("Installed!!!")
  }
}

# Launch from GitHub
# Load 'shiny' library
library(shiny)
Using runGitHub
runGitHub("EDC_shiny", "vittoriofortino84", subdir="INfORM-app")
# Using the archived file
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.tar.gz", subdir="EDC-shiny")
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.zip", subdir="EDC-shiny")
# Launch locally
# Clone the git repository
git clone https://github.com/vittoriofortino84/EDC_shiny EDC_shiny_clone
# Run by using runApp()
setwd("~/EDC-shiny")