
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

simPars1 <- read.csv("data/guidelines/SimPars2.0.csv") #90% umsy tracking - no abundance benchmark
simPars2 <- read.csv("data/guidelines/SimPars2.1.csv") #65% ER tracking - abundance benchmark reduces ER 90%
simPars3 <- read.csv("data/guidelines/SimPars2.2.csv")#90% umsy tracking - abundance benchmark reduces ER 90%
simPars4 <- read.csv("data/guidelines/SimPars2.3.csv") #90% umsy tracking - fixed abundance benchmark reduces ER 90%

cuPar <- read.csv("data/guidelines/CUPars2.0.csv")
#here()samsim_tv

#base case 
  
 





pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(1,nrow(simPars1))),
                 outname='umsy_track',
                 u=c(seq_len(nrow(simPars1))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls1',
                       nodes = 18, cpus_per_node = 1, submit = FALSE,
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
                 simfile=c(rep(3,nrow(simPars3))),
                 outname='umsy_bm_track',
                 u=c(seq_len(nrow(simPars3))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls3',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(4,nrow(simPars4))),
                 outname='umsy_fixedbm',
                 u=c(seq_len(nrow(simPars4))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls4',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")