#library(samEst)
#library(samSim)

# Manually define the inputs for row 3 (which contains u = 864)
#outpath <- "all_scenarios"
#simPars_path <- "../data/cls/SimPars.csv"
#cuPars_path <- "../data/cls/CUPars.csv"
#u_val <- 864
#n_trials <- 1000

# Load the specific CSVs
#cuPar <- read.csv(cuPars_path)
#simPars_all <- read.csv(simPars_path)

# Manually execute the function that failed
# This matches the function call seen in your original traceback
#genericRecoverySim(
#    simPar = simPars_all[u_val, ],
#    cuPar = cuPar,
#    catchDat = NULL,
#    srDat = NULL,
#    variableCU = FALSE,
#    ricPars = NULL,
#    larkPars = NULL,
#    cuCustomCorrMat = NULL,
#    outDir = outpath,
#    nTrials = n_trials,
#    makeSubDirs = TRUE,
#    random = FALSE,
#    uniqueProd = TRUE,
#    uniqueSurv = FALSE
#)
library(samEst)
library(samSim)
library(here)
source("R/func_sim.R")
id<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1
cuPar<-read.csv("../data/cls/CUPars.csv")
simPars<-read.csv("../data/cls/SimPars.csv")

#nrow(simPars)
#pars<-data.frame(outpath="all_scenarios",simPars="../data/cls/SimPars.csv",cuPars="../data/cls/CUPars.csv",u=id,n=1000)


samsim_tv(outpath="test",simPars="../data/cls/SimPars.csv",cuPars="../data/cls/CUPars.csv",u=id,n=1000)
