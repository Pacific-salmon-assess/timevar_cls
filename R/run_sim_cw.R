#=================================================
#run closed loop simulations
#=================================================

#setwd(paste0("C:\\Users\\worc\\Documents\\timevar\\samSim\\R","/.."))
#devtools::document()
#devtools::load_all()

#install samsim 
#
#remotes::install_github("Pacific-salmon-assess/samEst",  force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)


library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)
library(data.table)
library(stringi)

source("R/util_funcs.R")

#harrison example
simPars <- read.csv("data/cls/SimPars.csv")
cuPar <- read.csv("data/cls/CUPars.csv")


hcrDatalist<-list()
srData<- list()


simPars_test<-simPars[simPars$nameOM=="decLinearProd2to0.5",]

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
                      outDir="outs", 
                      nTrials=5, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)
  hcrDatalist[[a]] <-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/outs/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
                           paste(simPars_test$nameOM[a],"_", simPars_test$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout



hcrDatalist[[a]]$scenario<-simPars_test$scenario[a]
hcrDatalist[[a]]$nameOM<-simPars_test$nameOM[a]
hcrDatalist[[a]]$nameMP<-simPars_test$nameMP[a]

srData[[a]]<-readRDS(paste0("C:/Users/worc/Documents/timevar/timevar_cls/outs/SamSimOutputs/simData/", simPars_test$nameOM[a],"/",simPars_test$scenario[a],"/",
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

frass<-unique(srdat$freq_assess )

 spawn_plotlist<-list()
  catch_plotlist<-list()
  u_plotlist<-list()

  smsy_plotlist<-list()
  umsy_plotlist<-list()
  sgen_plotlist<-list()


for(rp in seq_along(frass)){
    
    srdat_plot_freq_assess<-srdat[  srdat$freq_assess==frass[rp],]

    hcrdat_plot_freq_assess<-hcrdat[ hcrdat$freq_assess==frass[rp],]

   
     
  spawn_plotlist[[rp]]<-ggplot(srdat_plot_freq_assess, aes(x=year, y=spawners,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
 # stat_summary(
 #   fun.data = function(x) {
 #     qs <- quantile(x, c(0.10, 0.5, 0.90))
 #     data.frame(ymin = qs[1],ymax = qs[3])
 #   },
 #   geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Spawner Abundance for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp], "model"))

  #catch_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=totalCatch,
  #  colour=freq_assess, fill = freq_assess,
  #  group = freq_assess)) +
  #geom_line(aes(group=iteration), alpha=0.3)+
  #stat_summary(
  #  fun.data = function(x) {
  #    qs <- quantile(x, c(0.10, 0.5, 0.90))
  #    data.frame(ymin = qs[1],ymax = qs[3])
  #  },
  #  geom = "ribbon",alpha = 0.2,colour = NA)+
  #stat_summary(
  #  fun = median, geom = "line", linewidth = 1) +
  #facet_grid(management_type~hcr)+
  #theme_minimal(base_size=16)+
  #coord_cartesian(ylim = c(0, 600000))+
  #scale_colour_viridis_d(end=.8) +
  #scale_fill_viridis_d(end=.8) +
  #labs(x = "Year", y = "Spawners", 
  #  title = paste("Total catch for scenario",unique(simPars_test$nameOM),"and ref pts from",unique(simPars_test$nameOM), "model"))

 
catch_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=totalCatch,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
 
  stat_summary(
    fun.data = function(x) {
      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  facet_grid(rp_type+management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 600000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

 u_plotlist[[rp]]<-ggplot(srdat_plot_freq_assess, aes(x=year, y=ER,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {
      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_freq_assess, aes(year,uMSY), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

  smsy_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=upperObsBM/.8,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_freq_assess, aes(year,sMSY), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y =  expression(paste(S[MSY])), 
    title = paste( "Smsy estimates for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

  umsy_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=UmsyBM,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_freq_assess, aes(year,uMSY), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = expression(paste(U[MSY])), 
    title = paste("Umsy estimates for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))


  sgen_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=lowerObsBM,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_freq_assess, aes(year,sGen), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 50000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = expression(paste(S[gen])), 
    title = paste("Sgen estimates for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

  }



ggplot(srdat_plot_freq_assess, aes(x=year, y=spawners,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {
      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 
  labs(x = "Year", y = "Spawners", 
    title = paste("Spawner Abundance for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

  catch_plotlist[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=totalCatch,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {
      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) 
  labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",unique(simPars_test$nameOM),"and assessment every",frass[rp]))

}



all_plots <- c(spawn_plotlist, catch_plotlist 
   status_plotlist,smsy_plotlist, umsy_plotlist,sgen_plotlist)
pdf(paste0("figs_brainstorm/",scn[sc],"_freqassess_plots.pdf"), width = 16, height = 12)
invisible(lapply(all_plots, print))
dev.off()



ggplot(hcrdat, aes(x=year, y=totalCatch,
    colour=as.factor(iteration), group= as.factor(iteration),
    group = freq_assess)) +
  geom_line(alpha = 0.2)+
  #stat_summary(
  #  fun = median, geom = "line", linewidth = 1) +
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 200000))+
  #scale_colour_viridis_c(end=.8) +
  #scale_fill_viridis_c(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",simPars$scenario[a]))






  
hcrDatalist[[a]]$aboveLowerassess<-"correct below upper BM"
hcrDatalist$aboveLowerassess[hcrDatalist$aboveLowerBM==1& hcrDatalist$aboveLowerObsBM==1]<-"correct above upper BM"
hcrDatalist$aboveLowerassess[hcrDatalist$aboveLowerBM==0& hcrDatalist$aboveLowerObsBM==1]<-"wrong optimistic"
hcrDatalist$aboveLowerassess[hcrDatalist$aboveLowerBM==1& hcrDatalist$aboveLowerObsBM==0]<-"wrong pessimistic"

names(hcrDatalist)

head(hcrDatalist$uMSyEst)
#plot UMSY
hcrDatalist[55:60,]

ggplot(hcrDatalist)+
geom_boxplot(aes(x=as.factor(year), y=uMSyEst))+
#scale_color_viridis_d(begin=.1, end=.8) +
#scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
#ggtitle("proportion above lower BM")+
theme_bw(14)



ggplot(hcrDatalist)+
geom_boxplot(aes(x=as.factor(year), y=sMSYEst))+
#scale_color_viridis_d(begin=.1, end=.8) +
#scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
#ggtitle("proportion above lower BM")+
theme_bw(14)


#plot alpha and 

#plot benchmarks

ggplot(hcrDatalist)+
geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
facet_wrap(~scenario)+
ggtitle("proportion above lower BM")+
theme_bw(14)








hcrdat <- do.call(rbind,hcrDatalist)
srdat<- do.call(rbind,srData)



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
