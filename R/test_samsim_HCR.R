#=================================================
#run closed loop simulations
#=================================================

library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)
library(stringi)
source("R/util_funcs.R")
path="."
  
  #gudelines scenarios example
  
simPars <- read.csv("data/cls/SimPars.csv")
 cuPars <- read.csv("data/cls/CUPars.csv") #cu pars
   
names(simPars)   

unique(simPars$scenario)
u=5

compsimPars<-simPars[simPars$scenario %in% c(
  "stationarylAR1_10yr_rwa_HCR1_forecast",
                         "stationarylAR1_10yr_rwa_HCR3_forecast",
                        "stationarylAR1_10yr_rwa_HCR2_forecast",
                        "stationarylAR1_10yr_rwa_HCR4_forecast",
                        "stationarylAR1_10yr_rwa_HCR1_retro",
                      "stationarylAR1_10yr_rwa_HCR3_retro",
                        "stationarylAR1_10yr_rwa_HCR2_retro",
                        "stationarylAR1_10yr_rwa_HCR4_retro"),]




for(u in 1:8){
  genericRecoverySim(simPar=compsimPars[u,], 
                      cuPar=cuPars, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="banana", 
                      nTrials=30, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)

}
  
simPar=compsimPars[u,]

ifelse(is.null(simPar$cvERSMU),NA,simPar$cvERSMU)

cuPar$cvER

cuPar=cuPars 
catchDat=NULL 
srDat=NULL
variableCU=FALSE 
ricPars=NULL 
larkPars=NULL 
cuCustomCorrMat= NULL
outDir="banana" 
nTrials=1 
makeSubDirs=TRUE 
random=FALSE 
uniqueProd=TRUE
uniqueSurv=FALSE


 #read these in, check if they are the same. 


hcrDatalist<-list()
srData<-list()

for(a in seq_len(nrow(compsimPars))){

  hcrDatalist[[a]] <- tryCatch(readRDS(paste0("./banana/SamSimOutputs/simData/",
                                       compsimPars$nameOM[a],"/", 
                                       compsimPars$scenario[a],"/",
                                       paste(compsimPars$nameOM[a],"_", compsimPars$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
                                       ,
  
                                
              error = function(e) {
                  message( e$message)
                  -99
              }
  )
  


  hcrDatalist[[a]]$scenario<-compsimPars$scenario[a]
  hcrDatalist[[a]]$nameOM<-compsimPars$nameOM[a]
  hcrDatalist[[a]]$nameMP<-compsimPars$nameMP[a]

  
  srData[[a]] <- readRDS(paste0("./banana/SamSimOutputs/simData/", 
                                  compsimPars$nameOM[a],"/",
                                  compsimPars$scenario[a],"/",
                                  paste(compsimPars$nameOM[a],"_", compsimPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
   
  srData[[a]]$scenario<-compsimPars$scenario[a]
  srData[[a]]$nameOM<-compsimPars$nameOM[a]
  srData[[a]]$nameMP<-compsimPars$nameMP[a]
  
}

hcrdatdbg <- data.table::rbindlist(hcrDatalist)#do.call(rbind,hcrDatalist)
srdatdbg<-data.table::rbindlist(srData)#srdat<- do.call(rbind,srData)

hcrdatdbg<-hcrdatdbg[hcrdatdbg$year>50,]
srdatdbg<-srdatdbg[srdatdbg$year>50,]

hcrdatdbg<-add_status_format(hcrdatdbg)

#regex to split MP into its parts
pat <- "^([^_]+)_((?:both(?:_tv_u_smsy)?|autocorr|rwa))_(HCR[1-4])_(retro|forecast)$"

splitMP <- stri_match_first_regex(srdatdbg$nameMP, pat)


srdatdbg$freq_assess     = splitMP[,2]
srdatdbg$rp_type         = splitMP[,3]
srdatdbg$hcr             = splitMP[,4]
srdatdbg$management_type = splitMP[,5]
 
splitMPhcrdat <- stri_match_first_regex(hcrdatdbg$nameMP, pat)


hcrdatdbg$freq_assess     = splitMPhcrdat[,2]
hcrdatdbg$rp_type         = splitMPhcrdat[,3]
hcrdatdbg$hcr             = splitMPhcrdat[,4]
hcrdatdbg$management_type = splitMPhcrdat[,5]
 

head(srdatdbg)
srdatdbg_shift <- srdatdbg[srdatdbg$year>51,]
srdatdbg_shift$spawners_lastyr <- srdatdbg$spawners[srdatdbg$year<120]
srdatdbg_shift$obsSpawners_lastyr <- srdatdbg$obsSpawners[srdatdbg$year<120]
srdatdbg_shift$recruits_lastyr <- srdatdbg$recruits[srdatdbg$year<120]  
  



ggplot(srdatdbg_shift,
      aes(x=spawners_lastyr,y= ER,colour=hcr, fill = hcr,
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


ggplot(srdatdbg_shift,
      aes(x=obsSpawners_lastyr,y= ER,colour=hcr, fill = hcr,
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


ggplot(srdatdbg_shift,
      aes(x=recruits_lastyr,y= ER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatdbg$capacity)*5))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 




ggplot(srdatdbg, aes(x=year, spawners,
    colour=freq_assess, fill = freq_assess,
    group = iteration) )+
  geom_line(linewidth = 1)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners")


  summspwdat <- srdat_plot_freq_assess|>
  group_by(year, freq_assess, management_type, hcr) |>
  summarise(
    q10 = quantile(spawners, 0.10),
    q50 = quantile(spawners, 0.50),
    q90 = quantile(spawners, 0.90),
    .groups = "drop"
  )
    
  spawn_plotlist_freq_assess[[rp]]<-ggplot(summspwdat, aes(x=year, q50,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 
  labs(x = "Year", y = "Spawners", 
    title = paste("Spawner Abundance for scenario",scn[sc],"and ref pts from",rps[rp], "model"))
