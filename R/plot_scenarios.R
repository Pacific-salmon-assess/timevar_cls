#=================================================
#run closed loop simulations
#=================================================

#
#remotes::install_github("Pacific-salmon-assess/samEst",  force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", force=TRUE)

library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringi)

source("R/util_funcs.R")

simPars <- read.csv("data/cls/SimPars.csv")
cuPar <- read.csv("data/cls/CUPars.csv")

head(simPars)
simPars_test<-simPars[simPars$nameMP=="1yr_autocorr_HCR3_retro",]


hcrDatalist<-list()
srData<- list()


for(a in seq_len(nrow(simPars_test))){
#a=42
  print(paste("scenario", a))
  genericRecoverySim(simPar=simPars_test[a,], 
                      cuPar=cuPar, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="scndraw", 
                      nTrials=1, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)
  hcrDatalist[[a]] <-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/scndraw/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
                           paste(simPars_test$nameOM[a],"_", simPars_test$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout



hcrDatalist[[a]]$scenario<-simPars_test$scenario[a]
hcrDatalist[[a]]$nameOM<-simPars_test$nameOM[a]
hcrDatalist[[a]]$nameMP<-simPars_test$nameMP[a]

srData[[a]]<-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/scndraw/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
                         paste(simPars_test$nameOM[a],"_", simPars_test$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

srData[[a]]$scenario<-simPars_test$scenario[a]
srData[[a]]$nameOM<-simPars_test$nameOM[a]
srData[[a]]$nameMP<-simPars_test$nameMP[a]

}


hcrdat <- data.table::rbindlist(hcrDatalist)#do.call(rbind,hcrDatalist)
srdat<-data.table::rbindlist(srData)#srdat<- do.call(rbind,srData)


hcrdat<-hcrdat[hcrdat$year>50,]
srdat<-srdat[srdat$year>50,]


hcrdat<-add_status_format(hcrdat)

#regex to split MP into its parts
pat <- "^([^_]+)_((?:both(?:_tv_u_smsy)?|autocorr|rwa))_(HCR[1-4])_(retro|forecast)$"

splitMP <- stri_match_first_regex(srdat$nameMP, pat)


srdat$freq_assess     = splitMP[,2]
srdat$rp_type         = splitMP[,3]
srdat$hcr             = splitMP[,4]
srdat$management_type = splitMP[,5]
 

splitMPhcrdat <- stri_match_first_regex(hcrdat$nameMP, pat)


hcrdat$freq_assess     = splitMPhcrdat[,2]
hcrdat$rp_type         = splitMPhcrdat[,3]
hcrdat$hcr             = splitMPhcrdat[,4]
hcrdat$management_type = splitMPhcrdat[,5]


#combine both


shared_cols<-names(hcrdat)[names(hcrdat) %in% names(srdat)]
shared_cols<-shared_cols[-3]

alldat <- inner_join(hcrdat, srdat, by = shared_cols)


#add scenario classification
#make df with year, sceanrio, and value of alph aor beta 

#plot param trajectories
head(alldat)




ggplot(alldat,
      aes(x=year,y= alpha))+
    facet_grid(~nameOM )+
    geom_line(alpha=.2)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 4))+
  scale_colour_viridis_c(end=.8) +
  scale_fill_viridis_c(end=.8) 



#make curves




ggplot(alldat,
      aes(x=spawners,y= recruits,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatdbg$capacity)*2))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 

