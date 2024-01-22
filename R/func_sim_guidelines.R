#=================================================
#run closed loop simulations
#=================================================

library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)


samsim_tv <- function(path=".",simfile,u, n){
  
  #gudelines scenarios example
  if(simfile==1){
    simPars<- read.csv(paste0(path,"/data/guidelines/SimPars1.5.csv"))
    cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars1.5.csv"))
  
  }else if(simfile==2){
     simPars <- read.csv(paste0(path,"/data/guidelines/Simpars2.0.csv"))
     cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
  }
  

  genericRecoverySim(simPar=simParsu[u,], 
                      cuPar=cuPar, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir=paste0(path,"/gdlout"), 
                      nTrials=n, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



  return(NULL)
}



