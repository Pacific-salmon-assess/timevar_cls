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


here::here()
#gudelines scenarios example
simPars1.5 <- read.csv("data/guidelines/SimPars1.5.csv")
cuPar1.5 <- read.csv("data/guidelines/CUPars1.5.csv")

simPars2.0 <- read.csv("data/guidelines/SimPars2.0.csv")
cuPar2.0 <- read.csv("data/guidelines/CUPars2.0.csv")
#here()
hcrDatalist1.5<-list()
srData1.5<- list()

#1.5 initial alpha scenarios
for(a in seq_len(nrow(simPars1.5))){

  #a=3
  genericRecoverySim(simPar=simPars1.5[a,], 
                      cuPar=cuPar1.5, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="../timevar_cls/gdlout", 
                      nTrials=5, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



   hcrDatalist1.5[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/gdlout/SamSimOutputs/simData/",
   simPars1.5$nameOM[a],"/", 
   simPars1.5$scenario[a],"/",
   paste(simPars1.5$nameOM[a],"_", simPars1.5$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist1.5[[a]]$scenario<-simPars1.5$scenario[a]
   hcrDatalist1.5[[a]]$nameOM<-simPars1.5$nameOM[a]
   hcrDatalist1.5[[a]]$nameMP<-simPars1.5$nameMP[a]


   srData1.5[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/gdlout/SamSimOutputs/simData/", 
                         simPars1.5$nameOM[a],"/",
                         simPars1.5$scenario[a],"/",
                         paste(simPars1.5$nameOM[a],"_", simPars1.5$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout


  srData1.5[[a]]$scenario<-simPars1.5$scenario[a]
  srData1.5[[a]]$nameOM<-simPars1.5$nameOM[a]
  srData1.5[[a]]$nameMP<-simPars1.5$nameMP[a]
}

  
hcrDatalist2.0<-list()
srData2.0<- list()
#1.5 initial alpha scenarios
for(a in seq_len(nrow(simPars2.0))){

  #a=3
  genericRecoverySim(simPar=simPars2.0[a,], 
                      cuPar=cuPar2.0, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="../timevar_cls/gdlout", 
                      nTrials=5, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



   hcrDatalist2.0[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/gdlout/SamSimOutputs/simData/",
   simPars2.0$nameOM[a],"/", 
   simPars2.0$scenario[a],"/",
   paste(simPars2.0$nameOM[a],"_", simPars2.0$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist2.0[[a]]$scenario<-simPars2.0$scenario[a]
   hcrDatalist2.0[[a]]$nameOM<-simPars2.0$nameOM[a]
   hcrDatalist2.0[[a]]$nameMP<-simPars2.0$nameMP[a]


   srData2.0[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/gdlout/SamSimOutputs/simData/", 
                         simPars2.0$nameOM[a],"/",
                         simPars2.0$scenario[a],"/",
                         paste(simPars2.0$nameOM[a],"_", simPars2.0$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout


  srData2.0[[a]]$scenario<-simPars2.0$scenario[a]
  srData2.0[[a]]$nameOM<-simPars2.0$nameOM[a]
  srData2.0[[a]]$nameMP<-simPars2.0$nameMP[a]
}



head(srdat1.5)

print(srdat1.5, n=100)
hcrdat1.5 <- do.call(rbind,hcrDatalist1.5)
srdat1.5<- do.call(rbind,srData1.5)

hcrdat1.5<-hcrdat1.5[hcrdat1.5$year>54,]
srdat1.5<-srdat1.5[srdat1.5$year>54,]

hcrdat2.0 <- do.call(rbind,hcrDatalist2.0)
srdat2.0<- do.call(rbind,srData2.0)

hcrdat2.0<-hcrdat2.0[hcrdat2.0$year>54,]
srdat2.0<-srdat2.0[srdat2.0$year>54,]


hcrdat<-rbind(hcrdat1.5,hcrdat2.0)
srdat<-rbind(srdat1.5,srdat2.0)

print(hcrdat, n=100)
print(srdat, n=100)

hcrdat$aboveLowerassess<-"correct below lower BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==1]<-"correct above lower BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==0& hcrdat$aboveLowerObsBM==1]<-"wrong optimistic"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==0]<-"wrong pessimistic"

head(hcrdat1.5)

ggplot(hcrdat)+
geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
ggtitle("proportion above lower BM")+
theme_bw(14)



hcrdat$aboveUpperassess<-"correct below upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==1]<-"correct above upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==0& hcrdat$aboveUpperObsBM==1]<-"wrong optimistic"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==0]<-"wrong pessimistic"


ggplot(hcrdat)+
geom_bar(aes(x=year, fill=factor(aboveUpperassess)),position = "fill")+
scale_fill_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
ggtitle("poportion above upper BM")+
theme_bw(14)





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
