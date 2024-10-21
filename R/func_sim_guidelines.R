#=================================================
#run closed loop simulations
#=================================================

library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)


samsim_tv <- function(path=".",outname='ss_out',simfile,u, n){
  
  #gudelines scenarios example
  if(simfile==1){
    simPars<- read.csv(paste0(path,"/data/guidelines/SimPars2.0.csv"))
    cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
  
  }else if(simfile==2){
     simPars <- read.csv(paste0(path,"/data/guidelines/SimPars2.1.csv"))
     cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
  }else if(simfile==3){
  simPars <- read.csv(paste0(path,"/data/guidelines/SimPars2.2.csv"))
  cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
  }
  else if(simfile==4){
    simPars <- read.csv(paste0(path,"/data/guidelines/SimPars2.3.csv"))
    cuPar <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
  }
  

  genericRecoverySim(simPar=simPars[u,], 
                      cuPar=cuPar, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir=paste0(path,"/",outname), 
                      nTrials=n, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



  return(NULL)
}



