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
library(cowplot)




here::here()
#gudelines scenarios example
simPars1.5 <- read.csv("data/guidelines/SimPars1.5.csv")
cuPar1.5 <- read.csv("data/guidelines/CUPars1.5.csv")

simPars2.0 <- read.csv("data/guidelines/Simpars2.0.csv")
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
                      outDir="../timevar_cls/cterout", 
                      nTrials=25, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



   hcrDatalist1.5[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/cterout/SamSimOutputs/simData/",
   simPars1.5$nameOM[a],"/", 
   simPars1.5$scenario[a],"/",
   paste(simPars1.5$nameOM[a],"_", simPars1.5$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist1.5[[a]]$scenario<-simPars1.5$scenario[a]
   hcrDatalist1.5[[a]]$nameOM<-simPars1.5$nameOM[a]
   hcrDatalist1.5[[a]]$nameMP<-simPars1.5$nameMP[a]


   srData1.5[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/cterout/SamSimOutputs/simData/", 
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
                      outDir="../timevar_cls/cterout", 
                      nTrials=25, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)



   hcrDatalist2.0[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/cterout/SamSimOutputs/simData/",
   simPars2.0$nameOM[a],"/", 
   simPars2.0$scenario[a],"/",
   paste(simPars2.0$nameOM[a],"_", simPars2.0$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist2.0[[a]]$scenario<-simPars2.0$scenario[a]
   hcrDatalist2.0[[a]]$nameOM<-simPars2.0$nameOM[a]
   hcrDatalist2.0[[a]]$nameMP<-simPars2.0$nameMP[a]


   srData2.0[[a]] <- readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/cterout/SamSimOutputs/simData/", 
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


hcrdat$aboveLowerassess<-"correct below lower BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==1]<-"correct above lower BM"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==0& hcrdat$aboveLowerObsBM==1]<-"wrong optimistic"
hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==0]<-"wrong pessimistic"

head(hcrdat1.5)

abovelower<-ggplot(hcrdat)+
geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
ggtitle("proportion above lower BM")+
theme_bw(14)
ggsave("figs/LB_assess.png",plot=abovelower)


hcrdat$aboveUpperassess<-"correct below upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==1]<-"correct above upper BM"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==0& hcrdat$aboveUpperObsBM==1]<-"wrong optimistic"
hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==0]<-"wrong pessimistic"
summary(as.factor(hcrdat$aboveUpperassess))

unique(hcrdat$scenario)


aboveupper<-ggplot(hcrdat)+
geom_bar(aes(x=year, fill=as.factor(aboveUpperassess)),position = "fill")+
scale_fill_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
ggtitle("poportion above upper BM")+
theme_bw(14)
aboveupper
ggsave("figs/UB_assess.png",plot=aboveupper)

#==================================================
#red amber green.


names(hcrdat)


hcrdat$status<-NA
hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==1&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==1]<-"green"
hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==0]<-"amber"
hcrdat$status[hcrdat$aboveLowerBM==0&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==0&
              hcrdat$aboveUpperObsBM==0]<-"red"


hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==1&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==0]<-"pessimistic green-> amber"


hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==1]<-"optimistic amber -> green"


hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==0&
              hcrdat$aboveUpperObsBM==0]<-"pessimistic amber -> red"

hcrdat$status[hcrdat$aboveLowerBM==0&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==0]<-"optimistic red -> amber"

hcrdat$status[hcrdat$aboveLowerBM==1&
              hcrdat$aboveUpperBM==1&
              hcrdat$aboveLowerObsBM==0&
              hcrdat$aboveUpperObsBM==0]<-"pessimistic green -> red"

hcrdat$status[hcrdat$aboveLowerBM==0&
              hcrdat$aboveUpperBM==0&
              hcrdat$aboveLowerObsBM==1&
              hcrdat$aboveUpperObsBM==1]<-"optimistic red -> green"


unique(hcrdat$status)

hcrdat$status_agg<-hcrdat$status
hcrdat$status_agg[hcrdat$status %in% c( "pessimistic green-> amber",                       
  "pessimistic amber -> red",  
 "pessimistic green -> red" )]<- "pessimistic"
hcrdat$status_agg[hcrdat$status %in% c( "optimistic amber -> green", 
  "optimistic red -> amber",  
  "optimistic red -> green")]<- "optimistic"


hcrdat$status_agg<-factor(hcrdat$status_agg,levels=c("red",
                                                     "amber",
                                                     "green",
                                                     "pessimistic",
                                                     "optimistic"))

statusColsall <- c("#9A3F3F","gray10","gray90","#DFD98D","gray20","gray80","#8EB687","gray30","gray70")
hcrdat$status<-factor(hcrdat$status,levels=c("red",
                                                      "optimistic red -> amber",
                                                      "pessimistic amber -> red",                                                      
                                                     "amber",
                                                     "optimistic amber -> green",
                                                     "pessimistic green-> amber",
                                                     "green",
                                                     "optimistic red -> green",
                                                     "pessimistic green -> red"
                                                     ))




st_agg<-ggplot(hcrdat)+
geom_bar(aes(x=year, fill=status_agg),position = "fill")+
scale_fill_manual(values = 
statusCols)+
facet_grid(nameOM~nameMP)+
ggtitle("status")+
theme_bw(12)+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
st_agg



st_agg_all<-ggplot(hcrdat[hcrdat$year>69,])+
geom_bar(aes(x=nameMP, fill=status_agg),position = "fill")+
scale_fill_manual(values = 
statusCols)+
facet_grid(nameOM~.)+
ggtitle("status")+
theme_bw(12)+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


plot_grid(
  st_agg, st_agg_all,
  align = "h", axis = "bt",
  rel_widths = c(.75,.25)
)



ggplot(hcrdat)+
geom_bar(aes(x=year, fill=status),position = "fill")+
scale_fill_manual(values = 
statusColsall)+
facet_grid(nameOM~nameMP)+
ggtitle("status")+
theme_bw(14)




#==================================================
#true vs obs benchmarks

BMcomp<-data.frame(year=hcrdat$year,
  iteration=hcrdat$iteration,
  typeest=as.factor(rep(c("truth","obs"),each=length(hcrdat$year))),
  lowerBM=c(hcrdat$lowerBM,hcrdat$lowerObsBM),
  upperBM=c(hcrdat$upperBM,hcrdat$upperObsBM),
  nameOM=hcrdat$nameOM,
  nameMP=hcrdat$nameMP
  )


ggplot(BMcomp)+
geom_line(aes(x=year,y=upperBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
theme_bw(14)


ggplot(BMcomp)+
geom_line(aes(x=year,y=lowerBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d(begin=.1, end=.8) +
facet_grid(nameOM~nameMP)+
coord_cartesian(ylim=c(0,70000))+
theme_bw(14)

