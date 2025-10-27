
#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)


library(samEst)
library(samSim)
library(rslurm)
library(here)

source("R/func_sim_guidelines.R")

cuPar <- read.csv("data/guidelines/CUPars2.0.csv")



simPars_all <- read.csv("data/guidelines/SimParsHCR_all.csv")


#samsim_tv(path='.',outname='test3',simfile=1,u=4,n=30)




pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                 simfile=c(rep(3,nrow(simPars_all))),
                 cuPar=1,
                 outname='all_HCR',
                 u=c(seq_len(nrow(simPars_all))),
                 n=1000)

samsim_tv(path=".", simfile=3,cuPar=1,u=25,n=20,outname="test")
samsim_tv(path=".", simfile=5,cuPar=1,u=25,n=20,outname="test")
samsim_tv(path=".", simfile=6,cuPar=1,u=25,n=20,outname="test")


sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls3',
                       nodes = 198, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")



#old HCR ones with regimes of 10 years and staring trends in year 10 and shortere shift time series
#3-Step ER - green: Umsy, amber:50% UMSY, red:10% ER
simPars3 <- read.csv("data/guidelines/SimPars2.2.csv")


#samsim_tv(path='.',outname='test3',simfile=1,u=4,n=30)




pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                 simfile=c(rep(3,nrow(simPars3))),
                 cuPar=1,
                 outname='umsy_bm_track',
                 u=c(seq_len(nrow(simPars3))),
                 n=1000)

samsim_tv(path=".", simfile=3,cuPar=1,u=25,n=20,outname="test")
samsim_tv(path=".", simfile=5,cuPar=1,u=25,n=20,outname="test")
samsim_tv(path=".", simfile=6,cuPar=1,u=25,n=20,outname="test")


sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls3',
                       nodes = 66, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")



#2-Step ER non-precautionary - green and amber: Umsy red:10% ER

simPars5 <- read.csv("data/guidelines/SimPars2.4.csv")# umsy tracking - abundance benchmark reduces ER 90%

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                 simfile=c(rep(5,nrow(simPars5))),
                 cuPar=1,
                 outname='umsy_bm_track_noamber',
                 u=c(seq_len(nrow(simPars5))),
                 n=1000)




sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls5',
                       nodes = 66, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")

#2-Step ER - precautionary - green : Umsy, amber and red:10% ER
simPars6 <- read.csv("data/guidelines/SimPars2.5.csv")#

pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                 simfile=c(rep(6,nrow(simPars6))),
                 cuPar=1,
                 outname='umsy_bm_track_low_er_amber_red',
                 u=c(seq_len(nrow(simPars6))),
                 n=500)

simPars6[9,]
samsim_tv(path=".", simfile=6,cuPar=1,u=9,n=500,outname='umsy_bm_track_low_er_amber_red')
 
sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_cls6',
                       nodes = 66, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")




res <- get_slurm_out(sjobcls, outtype = 'table', wait = TRUE)