##List of packages to be installed in specific path

#Giving exampleof Catarina's custom path

#In order to set the custom R library location I follow the instructions on this video:  
#https://www.youtube.com/watch?v=sbp5Q8niTho
CWRlib<-"/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/Rlib"


packages_cw<- c("remotes",
     "TMB", 
    "tmbstan"
    "rstan",
    "RTMB",
    "ggplot2", 
    "tidyr", 
    "stringi",
    "dplyr",
    "cowplot",
    "here",
    "gsl",
    "magrittr",
    "mvtnorm",
    "Rcpp",
    "RcppArmadillo",
    "scales",
    "sn",
    "synchrony",
    "zoo",
    "tibble",
    "stats",
    "Rdpack",
    "here" 


    )


#install CRAN packages
install.packages(packages_cw)


#install custom (in development) packages 
remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
remotes::install_github("Pacific-salmon-assess/samSim",  force=TRUE)


#Proceed to run the code in slurm_scripts/parallel_run
# adjust slurm_scripts submit.sh accordingly