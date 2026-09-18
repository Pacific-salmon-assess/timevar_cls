#install samsim
#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)
#data.table


library(samEst)

library(samSim)
library(ggplot2)
library(dplyr)
library(cowplot)
library(data.table)
library(stringi)

source("R/util_funcs.R")

statusCols <- c("#9A3F3F","#DFD98D","#8EB687","gray75","gray25")

#guidelines scenarios - load data####
#simPars_um <- read.csv("data/guidelines/SimPars2.0.csv") #simpars for ER tracking umsy, no EG
cuPar <- read.csv("data/cls/CUPars.csv") #cu pars
simPars <- read.csv("data/cls/SimPars.csv") #simpars for ER tracking umsy, assessed EG, stepped HCR

hcrDatalist<-list()
srData<-list()



for(a in seq_len(nrow(simPars))){
  

  hcrDatalist[[a]] <- tryCatch(readRDS(paste0("./allscn/SamSimOutputs/simData/",
                                       simPars$nameOM[a],"/", 
                                       simPars$scenario[a],"/",
                                       paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
                                       ,
  
                                
              error = function(e) {
                  message( e$message)
                  -99
              }
  )
  


  hcrDatalist[[a]]$scenario<-simPars$scenario[a]
  hcrDatalist[[a]]$nameOM<-simPars$nameOM[a]
  hcrDatalist[[a]]$nameMP<-simPars$nameMP[a]

  
  srData[[a]] <- readRDS(paste0("./allscn/SamSimOutputs/simData/", 
                                  simPars$nameOM[a],"/",
                                  simPars$scenario[a],"/",
                                  paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
   
  srData[[a]]$scenario<-simPars$scenario[a]
  srData[[a]]$nameOM<-simPars$nameOM[a]
  srData[[a]]$nameMP<-simPars$nameMP[a]
  
}




###

#check which ones are misssing

presencelist<-list()

for(a in seq_len(nrow(simPars))){
  

  presencelist[[a]] <- tryCatch(readRDS(paste0("./test/SamSimOutputs/simData/",
                                       simPars$nameOM[a],"/", 
                                       simPars$scenario[a],"/",
                                       paste(simPars$nameOM[a],"_", simPars$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
                                       ,
  
                                
              error = function(e) {
                  message( e$message)
                  -99
              }
  )
  

}

(1:nrow(simPars))[sapply(presencelist, function(x) identical(x, -99))]


# this is taking too long, may need to filter by scenario

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
 
#need to add one column for each
#HCR
#asses_freq
#rp_type
#management_type forecast/retrospective


#compute AAV

aavdf <- hcrdat %>%
  group_by(scenario, iteration, nameOM, nameMP,  hcr, rp_type, freq_assess, management_type) %>%
  summarise(
    aav = sum(abs(diff(totalCatch))) / sum(totalCatch),
    .groups = "drop"
  )


hcrdat$freq_assess<-factor(hcrdat$freq_assess,levels=c("1yr","5yr","10yr"))
srdat$freq_assess<-factor(srdat$freq_assess,levels=c("1yr","5yr","10yr"))
aavdf$freq_assess<-factor(aavdf$freq_assess,levels=c("1yr","5yr","10yr"))

hcrdat$wsp.status<-factor(hcrdat$wsp.status,levels=c("red","amber","green"))

#add in scenarios

hcrdat$paramvary<-"alpha"
hcrdat$paramvary[hcrdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
hcrdat$paramvary[hcrdat$nameOM%in%c("decLinearcap0.25","incLinearcap2","regCap0.25")]<-"beta"

srdat$paramvary<-"alpha"
srdat$paramvary[srdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
srdat$paramvary[srdat$nameOM%in%c("decLinearcap0.25","incLinearcap2","regCap0.25")]<-"beta"

aavdf$paramvary<-"alpha"
aavdf$paramvary[aavdf$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
aavdf$paramvary[aavdf$nameOM%in%c("decLinearcap0.25","incLinearcap2","regCap0.25")]<-"beta"


hcrdat$paramdirection<-"decrease"
hcrdat$paramdirection[hcrdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
hcrdat$paramdirection[hcrdat$nameOM%in%c("incLinearcap2","incLinearProd2to3" ,"incLinearcap2")]<-"increase"

srdat$paramdirection<-"decrease"
srdat$paramdirection[srdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
srdat$paramdirection[srdat$nameOM%in%c("incLinearcap2","incLinearProd2to3" ,"incLinearcap2")]<-"increase"

aavdf$paramdirection<-"decrease"
aavdf$paramdirection[aavdf$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"none"
aavdf$paramdirection[aavdf$nameOM%in%c("incLinearcap2","incLinearProd2to3" ,"incLinearcap2")]<-"increase"


hcrdat$changetype<-"linear"
hcrdat$changetype[hcrdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"stable"
hcrdat$changetype[hcrdat$nameOM%in%c("regProd2to1", "regProd2to0.5","regCap0.25")]<-"regime"

srdat$changetype<-"linear"
srdat$changetype[srdat$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"stable"
srdat$changetype[srdat$nameOM%in%c("regProd2to1", "regProd2to0.5","regCap0.25")]<-"regime"

aavdf$changetype<-"linear"
aavdf$changetype[aavdf$nameOM%in%c("stationarylAR1","stationaryhAR1")]<-"stable"
aavdf$changetype[aavdf$nameOM%in%c("regProd2to1", "regProd2to0.5","regCap0.25")]<-"regime"


#save efforts half way
saveRDS(aavdf, "data/cls/aavdf_large_scn.rds")
saveRDS(hcrdat, "data/cls/hcrdat_large_scn.rds")
saveRDS(srdat, "data/cls/srdat_large_scn.rds")
#last run on 27/Jul/2026


#=================================================
#If above code has been run
#aavdf<-readRDS("data/cls/aavdf_large_scn.rds")
#hcrdat<-readRDS("data/cls/hcrdat_large_scn.rds")
#srdat<-readRDS("data/cls/srdat_large_scn.rds")



scn<-unique(hcrdat$nameOM)
#scn<-c("regProd2to0.5",
#"regCap0.25",
#"decLinearProd2to0.5",
#"decLinearcap0.25",
#"stationarylAR1")

#comparisons to do: compare each one of the managemet procedurs
#hcrdat$freq_assess     
#hcrdat$rp_type         
#hcrdat$hcr             
#hcrdat$management_type 



#For 4 scenarios

#in terms of status, spawner abundance, catch aav  

#Let AI do the rest of the cross comparisons

#compare assessment frequency
#which assess_freq is better?

rps<-unique(srdat$rp_type)


for(sc in seq_along(scn)){
  #sc<-4
  spawn_plotlist_freq_assess<-list()
  catch_plotlist_freq_assess<-list()
  aav_plotlist_freq_assess<-list()
  status_plotlist_freq_assess<-list()
  smsy_plotlist_freq_assess<-list()
  umsy_plotlist_freq_assess<-list()
  sgen_plotlist_freq_assess<-list()


  for(rp in seq_along(rps)){
    #rp=2
    srdat_plot_freq_assess<-srdat[srdat$nameOM%in%scn[sc]&
                  srdat$rp_type==rps[rp],]

    hcrdat_plot_freq_assess<-hcrdat[hcrdat$nameOM%in%scn[sc]&
                  hcrdat$rp_type==rps[rp],]

    aav_plot_freq_assess<-aavdf[aavdf$nameOM==scn[sc]&
                  aavdf$rp_type==rps[rp],]

    
  summspwdat <- srdat_plot_freq_assess|>
  group_by(year, freq_assess, management_type, hcr) |>
  summarise(
    q10 = quantile(spawners, 0.10),
    q50 = quantile(spawners, 0.50),
    q90 = quantile(spawners, 0.90),
    sMSY= unique(sMSY),
    .groups = "drop"
  )
    #summspwdat
  ylimlow<-min(summspwdat$q10)
  ylimhigh<-max(summspwdat$q90)

  spawn_plotlist_freq_assess[[rp]]<-ggplot(summspwdat, aes(x=year, q50,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1)+
  geom_line(aes(x=year, sMSY),colour="black")+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(ylimlow, ylimhigh))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Spawner Abundance for scenario",scn[sc],"and ref pts from",rps[rp], "model"))


  #catch data
  
  summcatdat <- hcrdat_plot_freq_assess|>
  group_by(year, freq_assess, management_type, hcr) |>
  summarise(
    q10 = quantile(totalCatch, 0.10),
    q50 = quantile(totalCatch, 0.50),
    q90 = quantile(totalCatch, 0.90),
    .groups = "drop"
  )
  ylimlow<-min(summcatdat$q10)
  ylimhigh<-max(summcatdat$q90)
  catch_plotlist_freq_assess[[rp]]<-ggplot(summcatdat, aes(x=year, y=q50,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1)+
    facet_grid(management_type~hcr)+
    theme_minimal(base_size=16)+
    coord_cartesian(ylim = c(ylimlow, ylimhigh))+
    scale_colour_viridis_d(end=.8) +
    scale_fill_viridis_d(end=.8) +
    labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",scn[sc],"and ref pts from",rps[rp], "model"))

  
   
  aav_plotlist_freq_assess[[rp]]<-ggplot(aav_plot_freq_assess)+
  geom_boxplot(aes(x = interaction(freq_assess, management_type),y=aav,colour=freq_assess,fill=management_type), alpha=.6, outliers=FALSE)+
  facet_grid(~hcr)+
  theme_minimal(base_size=16)+
  scale_colour_viridis_d(end=.8)+
  scale_fill_viridis_d(end=.8)+
   labs(x = "Year", y = "AAV",
    title = paste("AAV for scenario",scn[sc],"and ref pts from",rps[rp], "model"))+
       theme(axis.text.x = element_text(angle = 45, hjust = 1))



  status_plotlist_freq_assess[[rp]]<-ggplot(hcrdat_plot_freq_assess)+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = statusCols)+
  facet_grid(management_type+freq_assess~hcr)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  ggtitle( paste("status for",scn[sc],"and ref pts from",rps[rp], "model"))+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom")



  
  summsmsydat <- hcrdat_plot_freq_assess|>
  group_by(year, freq_assess, management_type, hcr) |>
  summarise(
    sMSYq10 = quantile(upperObsBM/.8, 0.10),
    sMSYq50 = quantile(upperObsBM/.8, 0.50),
    sMSYq90 = quantile(upperObsBM/.8, 0.90),
    .groups = "drop"
  )
    

  smsy_plotlist_freq_assess[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=upperObsBM/.8,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
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
  coord_cartesian(ylim = c(min(hcrdat_plot_freq_assess$upperObsBM/.8), max(hcrdat_plot_freq_assess$upperObsBM/.8)))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y =  expression(paste(S[MSY])), 
    title = paste( "Smsy estimates for scenario",scn[sc],"and ref pts from",rps[rp], "model"))


  umsy_plotlist_freq_assess[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=UmsyBM,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
  stat_summary(
    fun.data = function(x) { qs <- quantile(x, c(0.10, 0.5, 0.90))
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
    title = paste("Umsy estimates for scenario",scn[sc],"and ref pts from",rps[rp], "model"))


  sgen_plotlist_freq_assess[[rp]]<-ggplot(hcrdat_plot_freq_assess, aes(x=year, y=lowerObsBM,
    colour=freq_assess, fill = freq_assess,
    group = freq_assess)) +
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
  coord_cartesian(ylim = c(min(hcrdat_plot_freq_assess$lowerObsBM), max(hcrdat_plot_freq_assess$lowerObsBM)))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = expression(paste(S[gen])), 
    title = paste("Sgen estimates for scenario",scn[sc],"and ref pts from",rps[rp], "model"))

  }

all_plots <- c(spawn_plotlist_freq_assess, catch_plotlist_freq_assess, aav_plotlist_freq_assess, 
   status_plotlist_freq_assess,smsy_plotlist_freq_assess, umsy_plotlist_freq_assess,sgen_plotlist_freq_assess)
pdf(paste0("figs_brainstorm/assess_freq_comparison/",scn[sc],"_and_",rps[rp],"_freqassess_plots.pdf"), width = 16, height = 12)
invisible(lapply(all_plots, print))
dev.off()
  

 
}
   
#what is the best set of reference points? 
fqs<-unique(srdat$freq_assess)

for(sc in seq_along(scn)){
  #sc<-4
  spawn_plotlist_refpoint<-list()
  catch_plotlist_refpoint<-list()
  aav_plotlist_refpoint<-list()
  status_plotlist_refpoint<-list()
  smsy_plotlist_refpoint<-list()
  umsy_plotlist_refpoint<-list()
  sgen_plotlist_refpoint<-list()

  #spawn_plotlist_hcr<-list()
  #catch_plotlist_hcr<-list()
  #aav_plotlist_hcr<-list()
  #status_plotlist_hcr<-list()
  #smsy_plotlist_hcr<-list()
  #umsy_plotlist_hcr<-list()
  #sgen_plotlist_hcr<-list()

  for(fa in seq_along(fqs)){
    #fa=2
    srdat_plot_refpoint<-srdat[srdat$nameOM%in%scn[sc]&
                  srdat$freq_assess==fqs[fa],]

    hcrdat_plot_refpoint<-hcrdat[hcrdat$nameOM%in%scn[sc]&
                  hcrdat$freq_assess==fqs[fa],]

    aav_plot_refpoint<-aavdf[aavdf$nameOM==scn[sc]&
                  aavdf$freq_assess==fqs[fa],]

 
  
    
  summspwdat <- srdat_plot_refpoint|>
  group_by(year, rp_type, management_type, hcr) |>
  summarise(
    q10 = quantile(spawners, 0.10),
    q50 = quantile(spawners, 0.50),
    q90 = quantile(spawners, 0.90),
    sMSY= unique(sMSY),
    .groups = "drop"
  )
  ylimlow<-min(summspwdat$q10)
  ylimhigh<-max(summspwdat$q90)
   
  #summspwdat
  spawn_plotlist_refpoint[[fa]]<-ggplot(summspwdat, aes(x=year, q50,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1)+
  geom_line(aes(x=year, sMSY),colour="black")+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(ylimlow, ylimhigh))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = "Spawners", 
    title = paste("Spawner Abundance for scenario",scn[sc],"and assessment every",fqs[fa]))


  


  #catch data
  
  summcatdat <- hcrdat_plot_refpoint|>
  group_by(year, rp_type, management_type, hcr) |>
  summarise(
    q10 = quantile(totalCatch, 0.10),
    q50 = quantile(totalCatch, 0.50),
    q90 = quantile(totalCatch, 0.90),
    .groups = "drop"
  )

  ylimlow<-min(summcatdat$q10)
  ylimhigh<-max(summcatdat$q90)
  catch_plotlist_refpoint[[fa]]<-ggplot(summcatdat, aes(x=year, y=q50,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1)+
    facet_grid(management_type~hcr)+
    theme_minimal(base_size=16)+
    coord_cartesian(ylim = c(ylimlow, ylimhigh))+
    scale_colour_viridis_d(end=.8) +
    scale_fill_viridis_d(end=.8) +
    labs(x = "Year", y = "Spawners", 
    title = paste("Total catch for scenario",scn[sc],"and assessment every",fqs[fa]))

   
  aav_plotlist_refpoint[[fa]]<-ggplot(aav_plot_refpoint)+
     geom_boxplot(aes(x = interaction(rp_type, management_type),y=aav,colour=rp_type, fill= management_type),alpha=.6, outliers=FALSE)+
     facet_grid(~hcr)+
     theme_minimal(base_size=16)+
     scale_colour_viridis_d(end=.8)+
     scale_fill_viridis_d(end=.8)+
      labs(x = "Year", y = "AAV",
       title = paste("AAV for scenario",scn[sc],"and assessment every",fqs[fa]))+
       theme(axis.text.x = element_text(angle = 45, hjust = 1))

  


  status_plotlist_refpoint[[fa]]<-ggplot(hcrdat_plot_refpoint)+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = statusCols)+
  facet_grid(management_type+rp_type~hcr)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  ggtitle( paste("status for",scn[sc],"and assessment every",fqs[fa]))+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "bottom")



  
  summsmsydat <- hcrdat_plot_refpoint|>
  group_by(year, rp_type, management_type, hcr) |>
  summarise(
    sMSYq10 = quantile(upperObsBM/.8, 0.10),
    sMSYq50 = quantile(upperObsBM/.8, 0.50),
    sMSYq90 = quantile(upperObsBM/.8, 0.90),
    .groups = "drop"
  )
    

  smsy_plotlist_refpoint[[fa]]<-ggplot(hcrdat_plot_refpoint, aes(x=year, y=upperObsBM/.8,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_refpoint, aes(year,sMSY), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(min(hcrdat_plot_refpoint$upperObsBM/.8), max(hcrdat_plot_refpoint$upperObsBM/.8)))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y =  expression(paste(S[MSY])), 
    title = paste( "Smsy estimates for scenario","and assessment every",fqs[fa]))


  umsy_plotlist_refpoint[[fa]]<-ggplot(hcrdat_plot_refpoint, aes(x=year, y=UmsyBM,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) { qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_refpoint, aes(year,uMSY), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 1))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = expression(paste(U[MSY])), 
    title = paste("Umsy estimates for scenario",scn[sc],"and assessment every",fqs[fa]))


  sgen_plotlist_refpoint[[fa]]<-ggplot(hcrdat_plot_refpoint, aes(x=year, y=lowerObsBM,
    colour=rp_type, fill = rp_type,
    group = rp_type)) +
  stat_summary(
    fun.data = function(x) {      qs <- quantile(x, c(0.10, 0.5, 0.90))
      data.frame(ymin = qs[1],ymax = qs[3])
    },
    geom = "ribbon",alpha = 0.2,colour = NA)+
  stat_summary(
    fun = median, geom = "line", linewidth = 1) +
  geom_line(data=srdat_plot_refpoint, aes(year,sGen), color="black", linewidth=1.2)+
  facet_grid(management_type~hcr)+
  theme_minimal(base_size=16)+
  coord_cartesian(ylim = c(0, 50000))+
  scale_colour_viridis_d(end=.8) +
  scale_fill_viridis_d(end=.8) +
  labs(x = "Year", y = expression(paste(S[gen])), 
    title = paste("Sgen estimates for scenario",scn[sc],"and assessment every",fqs[fa]))

  }

all_plots <- c(spawn_plotlist_refpoint, catch_plotlist_refpoint, aav_plotlist_refpoint, 
   status_plotlist_refpoint,smsy_plotlist_refpoint, umsy_plotlist_refpoint,sgen_plotlist_refpoint)
pdf(paste0("figs_brainstorm/ref_point_comparison/",scn[sc],"_",fqs[fa],"_refpoint_plots.pdf"), width = 16, height = 12)
invisible(lapply(all_plots, print))
dev.off()
  


 
}



#################################################

