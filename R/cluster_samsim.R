
#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)


library(samEst)
library(samSim)
library(rslurm)

source("R/func_sim_guidelines.R")


simPars1.5 <- read.csv("data/guidelines/SimPars1.5.csv")
cuPar1.5 <- read.csv("data/guidelines/CUPars1.5.csv")

simPars2.0 <- read.csv("data/guidelines/Simpars2.0.csv")
cuPar2.0 <- read.csv("data/guidelines/CUPars2.0.csv")
#here()

#base case 
#tst<-tmb_func(path="/fs/vnas_Hdfo/comda/caw001/Documents/tvsimest/cluster-tvsimest",
#  a=5,
#  u=1)
  
 

pars<-data.frame(path="..",
  simfile=c(rep(1,nrow(simPars1.5)),rep(2,nrow(simPars2.0))),
  u=c(seq_len(nrow(simPars1.5)),seq_len(nrow(simPars2.0))),
  n=1000)


sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls',
                    nodes = 28, cpus_per_node = 1, submit = FALSE,
                    pkgs=c("samEst","samSim"),
                    rscript_path = "/home/caw001/Documents/tvsimest/timevar_cls",
                    libPaths="/gpfs/fs7/dfo/hpcmc/comda/caw001/Rlib/4.1")





