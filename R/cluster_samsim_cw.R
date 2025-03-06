
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

#simPars1 <- read.csv("data/guidelines/SimPars2.0.csv") #90% umsy tracking - no abundance benchmark
#simPars2 <- read.csv("data/guidelines/SimPars2.1.csv") #65% ER tracking - abundance benchmark reduces ER 90%
simPars3 <- read.csv("data/guidelines/SimPars2.2.csv")#90% umsy tracking - abundance benchmark reduces ER 90%
#simPars4 <- read.csv("data/guidelines/SimPars2.3.csv") #90% umsy tracking - fixed abundance benchmark reduces ER 90%


#samsim_tv(path='.',outname='test3',simfile=1,u=4,n=30)




pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/caw001/results/timevar_cls",
                 simfile=c(rep(3,nrow(simPars3))),
                 outname='umsy_bm_track',
                 u=c(seq_len(nrow(simPars3))),
                 n=500)

samsim_tv(simfile=3,cuPar=1,u=25,n=20,outname="test")


sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls3',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/caw001/results/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/caw001/Rlib/4.3")



pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(4,nrow(simPars4))),
                 cuPar=c(rep(1,nrow(simPars2))),
                 outname='umsy_fixedbm',
                 u=c(seq_len(nrow(simPars4))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls4',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")

cuPar2 <- read.csv("data/guidelines/CUPars1.2.csv")


pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                  simfile=c(rep(2,nrow(simPars2))),
                  cuPar=c(rep(2,nrow(simPars2))),
                  outname='bm_track_cu1.2',
                  u=c(seq_len(nrow(simPars2))),
                  n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'ss_cls2_cupar1.2',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                 simfile=c(rep(4,nrow(simPars4))),
                 cuPar=c(rep(2,nrow(simPars2))),
                 outname='umsy_fixedbm_cu1.2',
                 u=c(seq_len(nrow(simPars2))),
                 n=500)

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'ss_cls4_cupar1.2',
                       nodes = 28, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/dag004/results/cl-sims/",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/dag004/Rlib/4.1")

