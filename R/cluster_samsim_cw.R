
#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)

#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)

library(samEst)
library(samSim)
library(rslurm)
library(here)

source("R/func_sim.R")



cuPar <- read.csv("data/cls/CUPars.csv")
simPars_all <- read.csv("data/cls/SimPars.csv")




 omworks<-list()
for(u in seq_len(nrow(simPars_all))){
   omworks[[u]]<- tryCatch(

           samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=u,n=3),
             #warning = function(w) {
   
              # return(paste("A warning occurred:", conditionMessage(w)))
             #},
           
             error = function(e) { return(paste0("error in line ",u)) }
          )
   }

 samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=316,n=3)
samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=421,n=3)

samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=529,n=3)

samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=637,n=3)

samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=745,n=3)


simPars_all[316,]
#"error in line 313" "error in line 314" "error in line 315"
# [4] "error in line 316" "error in line 317" "error in line 318"
# [7] "error in line 319" "error in line 320" "error in line 321"
#[10] "error in line 322" "error in line 323" "error in line 324"
#[13] "error in line 421" "error in line 422" "error in line 423"
#[16] "error in line 424" "error in line 425" "error in line 426"
#[19] "error in line 427" "error in line 428" "error in line 429"
#[22] "error in line 430" "error in line 431" "error in line 432"
#[25] "error in line 529" "error in line 530" "error in line 531"
#[28] "error in line 532" "error in line 533" "error in line 534"
#[31] "error in line 535" "error in line 536" "error in line 537"
#[34] "error in line 538" "error in line 539" "error in line 540"
#[37] "error in line 637" "error in line 638" "error in line 639"
#[40] "error in line 640" "error in line 641" "error in line 642"
#[43] "error in line 643" "error in line 644" "error in line 645"
#[46] "error in line 646" "error in line 647" "error in line 648"
#[49] "error in line 745" "error in line 746" "error in line 747"
#[52] "error in line 748" "error in line 749" "error in line 750"
#[55] "error in line 751" "error in line 752" "error in line 753"
#[58] "error in line 754" "error in line 755" "error in line 756"

#debug which ones die not work
unlist(omworks)

simPars_all[864,]


samsim_tv(outpath="test",simPars="data/cls/SimPars.csv",cuPars="data/cls/CUPars.csv",u=35,n=3)


pars<-data.frame(outpath="all_scenarios",
                 simPars="../data/cls/SimPars.csv",
                 cuPars="../data/cls/CUPars.csv",
                 u=c(seq_len(nrow(simPars_all))),
                 n=1000)

    

sjobcls <- slurm_apply(samsim_tv, pars, jobname = 'samsim_clspaul',
                       nodes = 300, cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")




#stuff that failed

#rows that failed
failrun<-c(  70, 274, 349, 350, 355, 356, 357, 358, 379, 380, 511, 512, 513, 514)

#  run on second tim



simPars_all[failrun,"scenario"]

pars_fail<-data.frame(outpath="all_scenarios",
                 simPars="../data/cls/SimPars.csv",
                 cuPars="../data/cls/CUPars.csv",
                 u=failrun,
                 n=1000)


sjobcls_fail <- slurm_apply(samsim_tv, pars, jobname = 'samsim_fail',
                       nodes = length(failrun), cpus_per_node = 1, submit = FALSE,
                       pkgs=c("samEst","samSim","here"),
                       rscript_path = "/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                       libPaths="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/R_4.3_ubuntu2404/x86_64-pc-linux-gnu-library/4.3")



#=old sims stuff
cuPar <- read.csv("data/guidelines/CUPars2.0.csv")
simPars_all <- read.csv("data/guidelines/SimParsHCR_all.csv")


samsim_tv(path='.',outname='test3',simfile=7,u=39,n=3,cuPar=1)




pars<-data.frame(path="/gpfs/fs7/dfo/hpcmc/pfm/spfm100/caw001/timevar_cls",
                 simfile=c(rep(3,nrow(simPars_all))),
                 cuPar=1,
                 outname='all_HCR',
                 u=c(seq_len(nrow(simPars_all))),
                 n=1000)

head(pars)

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