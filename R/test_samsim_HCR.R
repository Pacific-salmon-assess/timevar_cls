#=================================================
#run closed loop simulations
#=================================================

library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)

path="."
  
  #gudelines scenarios example
  
    simPars <- read.csv(paste0(path,"/data/guidelines/SimPars2.4.csv"))
 
      cuPars <- read.csv(paste0(path,"/data/guidelines/CUPars2.0.csv"))
   

u=5

  genericRecoverySim(simPar=simPars[u,], 
                      cuPar=cuPars, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="test", 
                      nTrials=1, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)





 HCRresult<-readRDS(paste0("./test/SamSimOutputs/simData/",
                                       simPars$nameOM[u],"/", 
                                       simPars$scenario[u],"/",
                                       paste(simPars$nameOM[u],"_", simPars$nameMP[u], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout


SRresult<-readRDS(paste0("./test/SamSimOutputs/simData/", 
                                  simPars$nameOM[u],"/",
                                  simPars$scenario[u],"/",
                                  paste(simPars$nameOM[u],"_", simPars$nameMP[u], "_", "CUsrDat.RData",sep="")))$srDatout
  


print(SRresult,n=70)

head(SRresult)

  unique(HCRresult$year)


  HCRresult[HCRresult$year>54,]

    sum(HCRresult$totalCatch, na.rm=T)