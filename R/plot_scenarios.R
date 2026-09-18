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
  #genericRecoverySim(simPar=simPars_test[a,], 
  #                    cuPar=cuPar, 
  #                    catchDat=NULL, 
  #                    srDat=NULL,
  #                    variableCU=FALSE, 
  #                    ricPars=NULL, 
  #                    larkPars=NULL, 
  #                    cuCustomCorrMat= NULL,
  #                    outDir="scndraw", 
  #                    nTrials=1, 
  #                    makeSubDirs=TRUE, 
  #                    random=FALSE, 
  #                    uniqueProd=TRUE,
  #                    uniqueSurv=FALSE)
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
alldat$typeOM<-"tv productivity"
alldat$typeOM[alldat$nameOM%in%c("decLinearcap0.25x","regCap0.25x","incLinearcap2x")]<-"tv capacity"
alldat$typeOM[alldat$nameOM%in%c("decLinearProd0.25xShiftCap1.5x","decLinearProd0.5xdecLinearCap0.5x")]<-"tv both"
alldat$typeOM[alldat$nameOM%in%c("stationary0.4AR1","stationary0.8AR1")]<-"stationary"

#add scenario classification
#make df with year, sceanrio, and value of alph aor beta 

#plot param trajectories
head(alldat)

unique(alldat$nameOM)

paramdf <- melt(
  alldat,
  id.vars = c("year", "scenario", "freq_assess", "rp_type", "hcr", "nameOM", "nameMP","typeOM"),
  measure.vars = c("capacity", "alpha"),
  variable.name = "param.name",
  value.name = "param"
)

head(paramdf)

paramdf$param.name[param.name=="alpha",]<-"productivity"
unique(paramdf$nameOM)

paramdf$typeOM <- factor(paramdf$typeOM, levels = c("stationary", "tv productivity","tv capacity", "tv both"))
paramdf$nameOM <- factor(paramdf$nameOM, levels = c("stationary0.4AR1",
                                                     "stationary0.8AR1",
                                                     "decLinearProd0.25x",
                                                      "regProd0.5x",
                                                      "regProd0.25x",
                                                      "incLinearProd1.5x",
                                                      "decLinearcap0.25x",
                                                      "incLinearcap2x",
                                                      "regCap0.25x",
                                                        "decLinearProd0.25xShiftCap1.5x",
                                                         "decLinearProd0.5xdecLinearCap0.5x")

unique(paramdf$nameOM)
custom_labeller <- function(x) {
  ifelse(is.na(custom_labels[x]), x, custom_labels[x])
}

custom_labels <- c(
   "stationary0.4AR1"="stationary0.4AR1",
   "stationary0.8AR1"= "stationary0.8AR1",
  "decLinearcap0.25x"  ="decLinearcap0.25x",                 
  "decLinearProd0.25x" ="decLinearProd0.25x" ,             
  "incLinearProd1.5x" ="incLinearProd1.5x",
  "incLinearcap2x" = "incLinearcap2x",
  "regProd0.5x" ="regProd0.5x",
  "regProd0.25x"="regProd0.25x",
  "regCap0.25x"="regCap0.25x",
  "decLinearProd0.25xShiftCap1.5x"="decLinearProd0.25x\nShiftCap1.5x",
   "decLinearProd0.5xdecLinearCap0.5x"= "decLinearProd0.5x\ndecLinearCap0.5x",
  "decLinearProd0.25xShiftCap1.5x" = "Declining Linear Prod (0.25x)\nShift Capacity (1.5x)",
  "stationary0.4AR1" = "Stationary (0.4 AR1)"
  # add one entry per unique value in param.name / nameOM / whatever you're faceting on
)


ggplot(paramdf,
aes(x=year,y= param,colour=typeOM ))+
facet_grid(param.name~nameOM ,scales = "free_y", labeller = as_labeller(custom_labeller ))+
geom_line(linewidth=2, alpha=0.8)+
  theme_minimal(base_size=16)+
  labs(y = "Parameter") +
  theme(
    panel.grid.minor = element_blank(),  # removes minor gridlines
    panel.grid.major = element_blank(),   # removes major gridlines too, if desired
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top", 
    axis.ticks = element_line(color = "black"),
     axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 



#make curves

simData<-list()
actualSR<-list()
for(a in seq_len(nrow(simPars_test))){

  simData[[a]] <- readRDS(paste0("scndraw/SamSimOutputs/simData/", 
                          simPars_test$nameOM[a],"/",
                          simPars_test$scenario[a],"/",
                          paste(simPars_test$nameOM[a],"_", 
                          simPars_test$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout


  dat<-simData[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat$scenario <- simPars_test$scenario[a]
  simData[[a]]<-dat

  S <- seq(0,600000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPars_test$scenario[a],
      nameOM=simPars_test$nameOM[a])

}




srdf<-do.call(rbind,actualSR)
head(srdf)
SRexample<-  ggplot(srdf) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scencode)

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



