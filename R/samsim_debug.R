#=================================================
#run closed loop simulations
#=================================================

library(here)
setwd(paste0("C:\\Users\\worc\\Documents\\timevar\\samSim\\R","/.."))
devtools::document()
#devtools::load_all()
here()

#install samsim 
#
#remotes::install_github("Pacific-salmon-assess/samEst",  force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)

library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringi)

source("R/util_funcs.R")

simPars <- read.csv("data/cls/SimPars_debug.csv")
cuPar <- read.csv("data/cls/CUPars_debug.csv")

hcrDatalist<-list()
srData<- list()
simPars_test<-simPars[simPars$nameOM=="stationarylAR1"&grepl("autocorr", simPars$scenario) & grepl("5yr", simPars$scenario),]


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
                      outDir="hcrtest", 
                      nTrials=10, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)
  hcrDatalist[[a]] <-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/hcrtest/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
                           paste(simPars_test$nameOM[a],"_", simPars_test$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout



hcrDatalist[[a]]$scenario<-simPars_test$scenario[a]
hcrDatalist[[a]]$nameOM<-simPars_test$nameOM[a]
hcrDatalist[[a]]$nameMP<-simPars_test$nameMP[a]

srData[[a]]<-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/hcrtest/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
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

unique(srdat$year)

srdat_shift <- srdat[srdat$year>51,]
srdat_shift$spawners_lastyr <- srdat$spawners[srdat$year<120]
srdat_shift$obsSpawners_lastyr <- srdat$obsSpawners[srdat$year<120]
srdat_shift$recruits_lastyr <- srdat$recruits[srdat$year<120]  
srdat_shift <- srdat_shift[srdat_shift$management_type=="retro",]  

  
srdat_forecast <- srdat[srdat$management_type=="forecast",]

srdat_forecast$spawners_lastyr <- srdat_forecast$spawners
srdat_forecast$obsSpawners_lastyr <- srdat_forecast$obsSpawners
srdat_forecast$recruits_lastyr <- srdat_forecast$recruits  
  
srdat_comb<-rbind(srdat_forecast,srdat_shift)  

plot(srdat_forecast$spawners_lastyr,srdat_forecast$targetER)



ggplot(srdat_comb,
      aes(x=spawners_lastyr,y= targetER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY*.8), color = "black",alpha=.2)+
    geom_vline(aes(xintercept = sMSY*.8), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdat_forecast$capacity)*2))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 


srdat_forecast$forecastlow<-srdat_forecast$sGen/(1-0.1)
srdat_forecast$forecasthigh<-srdat_forecast$sMSY*.8/(1-srdat_forecast$uMSY)

head(srdat_forecast)

ggplot(srdat_forecast,
      aes(x=forecastRunsize,y= targetER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = forecastlow), color = "darkgoldenrod2")+
    geom_vline(aes(xintercept = forecasthigh), color = "darkgoldenrod2")+
    geom_vline(aes(xintercept = sMSY*.8), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdat_forecast$capacity)*2))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 




shared_cols<-names(hcrdat)[names(hcrdat) %in% names(srdat)]
shared_cols<-shared_cols[-3]

alldat <- inner_join(hcrdat, srdat, by = shared_cols)


alldat$forecastlow<-alldat$sGen/(1-0.1)
alldat$forecasthigh<-alldat $sMSY*.8/(1-alldat $uMSY)

names(alldat)
head(alldat)

calcTAC_fixedER(rec=100000000,canER=0.6,amER=0,ppnMixVec=1,cvERcan=0.0000000001,cvERam=0,maxER=.99)

ggplot(alldatramp ,
      aes(x=spawners,y= targetER,colour=year, fill = year,
    group = hcr))+
      facet_grid(~freq_assess)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY*.8), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  scale_colour_viridis_c(end=.8) +
  scale_fill_viridis_c(end=.8) 



ggplot(alldatramp ,
      aes(x=year,y= upperObsBM, group=iteration))+
    geom_line(alpha=0.2)
    geom_point(alpha=.2)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, max(alldatramp$capacity)))
  #(end=.8) +
  #scale_fill_viridis_c(end=.8) 



alldattest<-alldatramp[alldatramp$targetER<0.11,]
names(alldattest)


df<-alldattest[1]
reder<-.1

testER.iter <-reder+ max(df$UmsyBM-reder,0)*(max(df$forecastRunsize*(1-df$UmsyBM)-df$lowerObsBM,0)/(df$upperObsBM-df$lowerObsBM))
                df$forecastRunsize*(1-df$UmsyBM)<df$upperObsBM
           


                for(o in 1:20){
                  testER.iter <-reder+ max(df$UmsyBM-reder,0)*(max(df$forecastRunsize*(1-testER.iter)-df$lowerObsBM,0)/(df$upperObsBM-df$lowerObsBM))  
                  df$forecastRunsize*(1-testER.iter)<df$upperObsBM
                  print(testER.iter)
                  if(df$forecastRunsize*(1-testER.iter)>df$upperObsBM){
                    testER.iter<-testER.iter*1.1
                    #testER.iter <-reder+ max(df$UmsyBM-reder,0)*(max(df$forecastRunsize*(1-testER.iter)-df$lowerObsBM,0)/(df$upperObsBM-df$lowerObsBM))  
                    #print(testER.iter)
                    #break
                  }
                  #if(testER.iter<reder){testER.iter<-reder}
                  
                
                }


forecastupperBM<-df$upperObsBM/(1-df$UmsyBM)
forecastupperBM*(1-df$UmsyBM)

                  df$forecastRunsize*(1-testER.iter)<df$upperObsBM

alldatramp<-alldat[alldat$hcr=="HCR4"&alldat$management_type=="forecast"&
                   alldat$forecastRunsize*(1-alldat$targetER)>alldat$upperObsBM&
                   alldat$forecastRunsize*(1-alldat$UmsyBM)<alldat$upperObsBM,]
nrow(alldatramp)


alldatramp$estER<-alldatramp$UmsyBM*((alldatramp$forecastRunsize*(1-alldatramp$UmsyBM)-alldatramp$lowerObsBM)/(alldatramp$upperObsBM-alldatramp$lowerObsBM))

weirdb<-alldatramp[alldatramp$estER<alldatramp$targetER,]

df<-weirdb[1,]
names(df)
df$forecastRunsize
df$upperObsBM

df$forecastRunsize*(1-df$UmsyBM)<df$upperObsBM
df$forecastRunsize*(1-df$UmsyBM)<df$lowerObsBM

df$uMSyEst
df$UmsyBM
df$targetER

df$UmsyBM*(df$forecastRunsize*(1-df$UmsyBM)/df$upperObsBM)

((df$forecastRunsize*(1-df$UmsyBM)-df$lowerObsBM)/(df$upperObsBM-df$lowerObsBM))


ggplot(alldatramp ,
      aes(x=targetER,y= estER,colour=year, fill = year,
    group = hcr))+
geom_abline(slope = 1, intercept = 0, color = "black") +
  coord_equal()+
      facet_grid(~freq_assess)+
    geom_point(alpha=.2)+
  theme_minimal(base_size=16)+
  scale_colour_viridis_c(end=.8) +
  scale_fill_viridis_c(end=.8) 


ggplot(alldatramp ,
      aes(x=spawners,y= targetER,colour=year, fill = year,
    group = hcr))+
      facet_grid(~freq_assess)+
    #geom_vline(aes(xintercept = upperObsBM,alpha=year), color = "black" )+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY*.8), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatcomb$capacity)*2))+
  scale_colour_viridis_c(end=.8) +
  scale_fill_viridis_c(end=.8) 




testalldat<-alldat[alldat$targetER<0.5,]



testalldat[1:20,]

nrow(testalldat)
shared_cols<-c("iteration", "year", "scenario","nameOM",    "nameMP"    )


testhcrdat <- semi_join(hcrdat, testsrdat, by = shared_cols)
nrow(testhcrdat)


jhsr<-left_join(testhcrdat[,-c("CU")],testsrdat[,-c("CU")])

hcrdat$UmsyBM[1:30]
jhsr$year

ggplot(jhsr,
      aes(x=spawners,y= ER,colour=year, fill = year,
    group ))+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatcomb$capacity)*2))+
  scale_colour_viridis_c(end=.8) +
  scale_fill_viridis_c(end=.8) 



df<-jhsr[1,]
names(df)
df$forecastRunsize
df$upperObsBM

df$forecastRunsize*(1-df$uMSyEst)<df$upperObsBM
df$forecastRunsize*(1-df$uMSyEst)<df$lowerObsBM
df$uMSyEst
df$UmsyBM



estER<-df$UmsyBM*(df$forecastRunsize*(1-df$UmsyBM)/df$upperObsBM)

#foreRecRY[y, k]*(1-trendCanER.iter[y,k])<=lowerObsBM[y-1,k]
df$forecastRunsize*(1-estER)<df$lowerObsBM

df$targetER/df$uMSyEst

bmUMSY[y-1,k,n]*(foreRecRY[y, k]*(1-bmUMSY[y-1,k,n])/(upperObsBM[y-1,k]))
               
df$targetER


ggplot(srdat_forecast,
      aes(x=spawners,y= ER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatcomb$capacity)*2))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 



fora<-exp(qnorm(runif(10000, 0.0001, 0.9999),log(1), .2))

summary(fora)
hist(fora)




ggplot(srdat_forecast,
      aes(x=forecastRunsize,y= ER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdat$capacity)*4))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 




ggplot(srdat_forecast,
      aes(x=forecastRunsize,y=Runsize,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  #coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdat$capacity)*4))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 




ggplot(srdatcomb,
      aes(x=spawners,y= ER,colour=hcr, fill = hcr,
    group = hcr))+
      facet_grid(management_type~freq_assess+hcr)+
    geom_point(alpha=.2)+
    geom_vline(aes(xintercept = sMSY), color = "blue")+
    geom_vline(aes(xintercept = sGen), color = "blue")+
    geom_hline(aes(yintercept = uMSY), color = "blue")+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1),xlim=c(0,max(srdatcomb$capacity)*2))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 






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

