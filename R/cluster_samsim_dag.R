
#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)


library(samEst)
library(samSim)
library(rslurm)
library(here)

source("R/func_sim_guidelines.R")
#source("R/genericRecoverySimulator.R")


#simPars1.5 <- read.csv("data/guidelines/SimPars1.5.csv")
#cuPar1.5 <- read.csv("data/guidelines/CUPars1.5.csv")

simPars1 <- read.csv("data/guidelines/SimPars2.0.csv")
simPars2 <- read.csv("data/guidelines/SimPars2.1.csv")
simPars3 <- read.csv("data/guidelines/SimPars2.2.csv")
cuPar <- read.csv("data/guidelines/CUPars2.0.csv")
#here()

#base case 
#tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
#  a=5,
87#  u=1)
  
 



tst=samsim_tv(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims",simfile=2,outname='test',u=2,n=1)

a=2
tst.d <- readRDS(paste0("test/SamSimOutputs/simData/",
                                      simPars1$nameOM[a],"/", 
                                      simPars1$scenario[a],"/",
                                      paste(simPars1$nameOM[a],"_", simPars1$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

tst.d2 <- readRDS(paste0("test/SamSimOutputs/simData/", 
                                 simPars1$nameOM[a],"/",
                                 simPars1$scenario[a],"/",
                                 paste(simPars1$nameOM[a],"_", simPars1$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(1,nrow(simPars1))),
                 outname='umsy_track',
                 u=c(seq_len(nrow(simPars1))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls1',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")



pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(2,nrow(simPars2))),
                 outname='bm_track',
                 u=c(seq_len(nrow(simPars2))),
                 n=500)


sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls2',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(2,nrow(simPars3))),
                 outname='umsy_bm_track',
                 u=c(seq_len(nrow(simPars3))),
                 n=500)hjh

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls3',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")yui