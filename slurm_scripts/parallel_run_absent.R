
library(samEst)
library(samSim)
library(here)
source("../R/func_sim.R")
absent<-c(359, 466,  791, 1018, 1030, 1042, 1054
)
id<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1
cuPar<-read.csv("../data/cls/CUPars.csv")
simPars<-read.csv("../data/cls/SimPars.csv")


samsim_tv(outpath="absent",simPars="../data/cls/SimPars.csv",cuPars="../data/cls/CUPars.csv",u=absent[id],n=1000)

