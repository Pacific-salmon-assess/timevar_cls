#=================================================
#run closed loop simulations
#=================================================

#library(samEst)
#library(samSim)
#library(ggplot2)
#library(dplyr)


samsim_tv <- function(outpath="out",simPars,cuPars,u, n){
  
  #gudelines scenarios example
    
  

  genericRecoverySim(simPar=simPars[u,], 
                      cuPar=cuPars, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir=outpath, 
                      nTrials=n, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)

  return(NULL)
}



