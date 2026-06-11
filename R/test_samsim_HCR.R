#=================================================
#run closed loop simulations
#=================================================

library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)

path="."
  
  #gudelines scenarios example
  
simPars <- read.csv("data/cls/SimPars.csv")
 cuPars <- read.csv("data/cls/CUPars.csv") #cu pars
   
names(simPars)   

unique(simPars$scenario)
u=5

compsimPars<-simPars[simPars$scenario %in% c("decLinearProd2to0.5_10yr_rwa_HCR3_forecast",
                        "decLinearProd2to0.5_10yr_rwa_HCR2_forecast",
                        "decLinearProd2to0.5_10yr_rwa_HCR4_forecast"),]


u=1

for(u in 1:3){
  genericRecoverySim(simPar=compsimPars[u,], 
                      cuPar=cuPars, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="banana", 
                      nTrials=10, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)

}
  
simPar=compsimPars[u,]
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
 

ggplot(srdat, aes(x=year, spawners,
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
