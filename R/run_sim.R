#=================================================
#run closed loop simulations
#=================================================

#setwd(paste0("C:\\Users\\worc\\Documents\\timevar\\samSim\\R","/.."))
#devtools::document()
#devtools::load_all()

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)


library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)



#harrison example
simPars <- read.csv("data/HAR/SimPars.csv")
cuPar <- read.csv("data/HAR/CUPars.csv")
#here()
hcrDatalist<-list()
srData<- list()

for(a in seq_len(nrow(simPars))){

  #a=3
  genericRecoverySim(simPar=simPars[a,], 
                      cuPar=cuPar, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="../timevar_cls/outs", 
                      nTrials=20, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



  hcrDatalist[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                           paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist[[a]]$scenario<-simPars$scenario[a]
   hcrDatalist[[a]]$nameOM<-simPars$nameOM[a]
   hcrDatalist[[a]]$nameMP<-simPars$nameMP[a]


   srData[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout


  srData[[a]]$scenario<-simPars$scenario[a]
  srData[[a]]$nameOM<-simPars$nameOM[a]
  srData[[a]]$nameMP<-simPars$nameMP[a]
}

  





hcrdat <- do.call(rbind,hcrDatalist)
srdat<- do.call(rbind,srData)


hcrdat$aboveLowerassess<-"correct below upper BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==1]<-"correct above upper BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==0& hcrdat$aboveLowerObsBM==1]<-"wrong optimistic"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==0]<-"wrong pessimistic"



ggplot(hcrdat)+
geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
ggtitle("proportion above lower BM")+
theme_bw(14)



hcrdat$aboveUpperassess<-"correct below upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==1]<-"correct above upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==0& hcrdat$aboveUpperObsBM==1]<-"wrong optimistic"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==0]<-"wrong pessimistic"


ggplot(hcrdat)+
geom_bar(aes(x=year, fill=factor(aboveUpperassess)),position = "fill")+
scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
ggtitle("poportion above upper BM")+
theme_bw(14)




head(hcrdat)
head(srdat)

hcrdat<-hcrdat[hcrdat$year>50,]
summary(hcrData)

aboveUpperBM<- aggregate(hcrdat$aboveUpperBM, 
                          list(year=hcrdat$year,
                            scenario=hcrdat$scenario),
                           function(x){sum(x/length(x))})

aboveUpperObsBM<- aggregate(hcrdat$aboveUpperObsBM, 
                          list(year=hcrdat$year,
                            scenario=hcrdat$scenario),
                           function(x){sum(x/length(x))})

upperBMdf<-data.frame(year=aboveUpperBM$year,
                     scenario=aboveUpperBM$scenario,
                      aboveUpperBM=c(aboveUpperBM$x,aboveUpperObsBM$x),
                      type=rep(c("true","est"), each=nrow(aboveUpperBM)))

head(upperBMdf)

ggplot(upperBMdf)+
geom_line(aes(x=year,y=aboveUpperBM, color=type))+
scale_color_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
theme_bw(14)



ggplot(upperBMdf)+
geom_bar(aes(x=year,y=aboveUpperBM, color=type, fill=type, color=type), 
    stat="identity", position="dodge")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
theme_bw(14)


aboveLowerBM<- aggregate(hcrdat$aboveLowerBM, 
                          list(year=hcrdat$year,
                            scenario=hcrdat$scenario),
                           function(x){sum(x/length(x))})

aboveLowerObsBM<- aggregate(hcrdat$aboveLowerObsBM, 
                          list(year=hcrdat$year,
                            scenario=hcrdat$scenario),
                           function(x){sum(x/length(x))})

lowerBMdf<-data.frame(year=aboveLowerBM$year,
                     scenario=aboveLowerBM$scenario,
                     aboveLowerBM=c(aboveLowerBM$x,aboveLowerObsBM$x),
                      type=rep(c("true","est"), each=nrow(aboveLowerBM)))



ggplot(lowerBMdf)+
geom_line(aes(x=year,y=aboveLowerBM, color=type))+
scale_color_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
theme_bw(14)




ggplot(lowerBMdf)+
geom_bar(aes(x=year,y=aboveLowerBM, fill=type, color=type), 
    stat="identity", position="dodge")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
coord_cartesian(y=c(.85, 1))+
theme_bw(14)



srData <- readRDS(paste0("C:/Users/worc/Documents/timevar/SBC_chinook_cls/outs/SamSimOutputs/simData/", simPars$nameOM[a],"/",simPars$scenario[a],"/",
                         paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

head(srData)
srData<-srData[srData$year>50,]


print(srData,n=200)

ggplot(srData)+
geom_line(aes(x=year,y=ER, color=iteration,group=iteration))+
scale_color_viridis_c(begin=.1, end=.8) +
theme_bw(14)
