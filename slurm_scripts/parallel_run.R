
library(samEst)
library(samSim)
library(here)
source("../R/func_sim.R")
id<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1
cuPar<-read.csv("../data/cls/CUPars.csv")
simPars<-read.csv("../data/cls/SimPars.csv")


samsim_tv(outpath="test",simPars="../data/cls/SimPars.csv",cuPars="../data/cls/CUPars.csv",u=id,n=1000)

