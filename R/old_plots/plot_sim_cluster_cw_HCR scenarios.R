#install samsim
#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk-hatch", force=TRUE)


library(samEst)
library(samSim)
library(ggplot2)
library(dplyr)
library(cowplot)

source("R/util_funcs.R")

#guidelines scenarios - load data####
#simPars_um <- read.csv("data/guidelines/SimPars2.0.csv") #simpars for ER tracking umsy, no EG
cuPar <- read.csv("data/guidelines/CUPars2.0.csv") #cu pars
simPars_HCR1 <- read.csv("data/guidelines/SimPars2.2.csv") #simpars for ER tracking umsy, assessed EG, stepped HCR
simPars_HCR2 <- read.csv("data/guidelines/SimPars2.4.csv") #simpars for ER tracking umsy, assessed EG, lowE R if red
simPars_HCR3 <- read.csv("data/guidelines/SimPars2.5.csv") #simpars for ER tracking umsy, assessed EG lowE R if amber or red

hcrDatalist_HCR1<-list()
hcrDatalist_HCR2<-list()
hcrDatalist_HCR3<-list()

srData_HCR1<- list()
srData_HCR2<- list()
srData_HCR3<- list()

for(a in seq_len(nrow(simPars_HCR1))){
  
   hcrDatalist_HCR1[[a]] <- readRDS(paste0("./umsy_bm_track/SamSimOutputs/simData/",
                                       simPars_HCR1$nameOM[a],"/", 
                                       simPars_HCR1$scenario[a],"/",
                                       paste(simPars_HCR1$nameOM[a],"_", simPars_HCR1$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
   hcrDatalist_HCR2[[a]] <- readRDS(paste0("./umsy_bm_track_noamber/SamSimOutputs/simData/",
                                       simPars_HCR2$nameOM[a],"/", 
                                       simPars_HCR2$scenario[a],"/",
                                       paste(simPars_HCR2$nameOM[a],"_", simPars_HCR2$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
   hcrDatalist_HCR3[[a]] <- readRDS(paste0("./umsy_bm_track_low_er_amber_red/SamSimOutputs/simData/",
                                       simPars_HCR3$nameOM[a],"/", 
                                       simPars_HCR3$scenario[a],"/",
                                       paste(simPars_HCR3$nameOM[a],"_", simPars_HCR3$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  


  hcrDatalist_HCR1[[a]]$scenario<-simPars_HCR1$scenario[a]
  hcrDatalist_HCR1[[a]]$nameOM<-simPars_HCR1$nameOM[a]
  hcrDatalist_HCR1[[a]]$nameMP<-simPars_HCR1$nameMP[a]

  hcrDatalist_HCR2[[a]]$scenario<-simPars_HCR2$scenario[a]
  hcrDatalist_HCR2[[a]]$nameOM<-simPars_HCR2$nameOM[a]
  hcrDatalist_HCR2[[a]]$nameMP<-simPars_HCR2$nameMP[a]

  hcrDatalist_HCR3[[a]]$scenario<-simPars_HCR3$scenario[a]
  hcrDatalist_HCR3[[a]]$nameOM<-simPars_HCR3$nameOM[a]
  hcrDatalist_HCR3[[a]]$nameMP<-simPars_HCR3$nameMP[a]
  
  
  
  srData_HCR1[[a]] <- readRDS(paste0("./umsy_bm_track/SamSimOutputs/simData/", 
                                  simPars_HCR1$nameOM[a],"/",
                                  simPars_HCR1$scenario[a],"/",
                                  paste(simPars_HCR1$nameOM[a],"_", simPars_HCR1$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  


  srData_HCR2[[a]] <- readRDS(paste0("./umsy_bm_track_noamber/SamSimOutputs/simData/", 
                                  simPars_HCR2$nameOM[a],"/",
                                  simPars_HCR2$scenario[a],"/",
                                  paste(simPars_HCR2$nameOM[a],"_", simPars_HCR2$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  

  srData_HCR3[[a]] <- readRDS(paste0("./umsy_bm_track_low_er_amber_red/SamSimOutputs/simData/", 
                                 simPars_HCR3$nameOM[a],"/",
                                 simPars_HCR3$scenario[a],"/",
                                 paste(simPars_HCR3$nameOM[a],"_", simPars_HCR3$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  
  srData_HCR1[[a]]$scenario<-simPars_HCR1$scenario[a]
  srData_HCR1[[a]]$nameOM<-simPars_HCR1$nameOM[a]
  srData_HCR1[[a]]$nameMP<-simPars_HCR1$nameMP[a]

  srData_HCR2[[a]]$scenario<-simPars_HCR2$scenario[a]
  srData_HCR2[[a]]$nameOM<-simPars_HCR2$nameOM[a]
  srData_HCR2[[a]]$nameMP<-simPars_HCR2$nameMP[a]

  srData_HCR3[[a]]$scenario<-simPars_HCR3$scenario[a]
  srData_HCR3[[a]]$nameOM<-simPars_HCR3$nameOM[a]
  srData_HCR3[[a]]$nameMP<-simPars_HCR3$nameMP[a]
}

hcrdat_HCR1 <- do.call(rbind,hcrDatalist_HCR1)
srdat_HCR1<- do.call(rbind,srData_HCR1)

hcrdat_HCR2 <- do.call(rbind,hcrDatalist_HCR2)
srdat_HCR2<- do.call(rbind,srData_HCR2)

hcrdat_HCR3 <- do.call(rbind,hcrDatalist_HCR3)
srdat_HCR3 <- do.call(rbind,srData_HCR3)

hcrdat_HCR1<-hcrdat_HCR1[hcrdat_HCR1$year>50,]
srdat_HCR1<-srdat_HCR1[srdat_HCR1$year>50,]
srdat_HCR1$HCR<-1

hcrdat_HCR2<-hcrdat_HCR2[hcrdat_HCR2$year>50,]
srdat_HCR2<-srdat_HCR2[srdat_HCR2$year>50,]
srdat_HCR2$HCR<-2

hcrdat_HCR3<-hcrdat_HCR3[hcrdat_HCR3$year>50,]
srdat_HCR3<-srdat_HCR3[srdat_HCR3$year>50,]
srdat_HCR3$HCR<-3

srdat_all<-rbind(srdat_HCR1,srdat_HCR2,srdat_HCR3)


scnwboth<-unique(hcrdat_HCR1$scenario)[grep(pattern="both", unique(hcrdat_HCR1$scenario))]
nameMPboth<-gsub("rwa", "both", hcrdat_HCR1$nameMP[hcrdat_HCR1$scenario%in%scnwboth])
hcrdat_HCR1$nameMP[hcrdat_HCR1$scenario%in%scnwboth]<-nameMPboth

scnwboth<-unique(hcrdat_HCR3$scenario)[grep(pattern="both", unique(hcrdat_HCR3$scenario))]
nameMPboth<-gsub("rwa", "both", hcrdat_HCR3$nameMP[hcrdat_HCR3$scenario%in%scnwboth])
hcrdat_HCR3$nameMP[hcrdat_HCR3$scenario%in%scnwboth]<-nameMPboth


srwboth<-unique(srdat_HCR1$scenario)[grep(pattern="both", unique(srdat_HCR1$scenario))]
nameMPboth_sr<-gsub("rwa", "both", srdat_HCR1$nameMP[srdat_HCR1$scenario%in%srwboth])
srdat_HCR1$nameMP[srdat_HCR1$scenario%in%srwboth]<-nameMPboth_sr


srwboth<-unique(srdat_HCR3$scenario)[grep(pattern="both", unique(srdat_HCR3$scenario))]
nameMPboth_sr<-gsub("rwa", "both", srdat_HCR3$nameMP[srdat_HCR3$scenario%in%srwboth])
srdat_HCR3$nameMP[srdat_HCR3$scenario%in%srwboth]<-nameMPboth_sr

hcrdat1=format_dat(hcrdat_HCR1)
hcrdat2=format_dat(hcrdat_HCR2)
hcrdat3=format_dat(hcrdat_HCR3)

hcrdat1$HCR<-1
hcrdat2$HCR<-2
hcrdat3$HCR<-3

hcrdat_all<-rbind(hcrdat1,hcrdat2,hcrdat3)


hcrdat1$plotOM<-dplyr::recode_factor(factor(hcrdat1$plotOM),
                                     "AR1 low" = "Static - 0.4AR1",
                                     "AR1 high" = "Static - 0.8AR1",
                                     "cap * 0.5"="-50%Cap", 
                                     "cap * 0.25"="-75%Cap",
                                     "cap * 1.5"="+50%Cap",
                                     "prod 2 -> 3"="+50%Prod",
                                     "prod 2 -> 1.5"="-25%Prod",     
                                     "prod 2 -> 1"="-50%Prod",
                                     "prod 2 -> 0.5"="-75%Prod",
                                     "regProd2to1"='-50%PRegime',
                                     "regProd2to0.5"='-75%PRegime')

hcrdat2$plotOM<-dplyr::recode_factor(factor(hcrdat2$plotOM),
                                     "AR1 low" = "Static - 0.4AR1",
                                     "AR1 high" = "Static - 0.8AR1",
                                     "cap * 0.5"="-50%Cap", 
                                     "cap * 0.25"="-75%Cap",
                                     "cap * 1.5"="+50%Cap",
                                     "prod 2 -> 3"="+50%Prod",
                                     "prod 2 -> 1.5"="-25%Prod",     
                                     "prod 2 -> 1"="-50%Prod",
                                     "prod 2 -> 0.5"="-75%Prod",
                                     "regProd2to1"='-50%PRegime',
                                     "regProd2to0.5"='-75%PRegime')


hcrdat3$plotOM<-dplyr::recode_factor(factor(hcrdat3$plotOM),
                                     "AR1 low" = "Static - 0.4AR1",
                                     "AR1 high" = "Static - 0.8AR1",
                                     "cap * 0.5"="-50%Cap", 
                                     "cap * 0.25"="-75%Cap",
                                     "cap * 1.5"="+50%Cap",
                                     "prod 2 -> 3"="+50%Prod",
                                     "prod 2 -> 1.5"="-25%Prod",     
                                     "prod 2 -> 1"="-50%Prod",
                                     "prod 2 -> 0.5"="-75%Prod",
                                     "regProd2to1"='-50%PRegime',
                                     "regProd2to0.5"='-75%PRegime')



hcrdat_all$plotOM<-dplyr::recode_factor(factor(hcrdat_all$plotOM),
                                     "AR1 low" = "Static - 0.4AR1",
                                     "AR1 high" = "Static - 0.8AR1",
                                     "cap * 0.5"="-50%Cap", 
                                     "cap * 0.25"="-75%Cap",
                                     "cap * 1.5"="+50%Cap",
                                     "prod 2 -> 3"="+50%Prod",
                                     "prod 2 -> 1.5"="-25%Prod",     
                                     "prod 2 -> 1"="-50%Prod",
                                     "prod 2 -> 0.5"="-75%Prod",
                                     "regProd2to1"='-50%PRegime',
                                     "regProd2to0.5"='-75%PRegime')



hcrdat_all$plotMP=dplyr::case_match(factor(hcrdat_all$plotMP),
                         "10yr_autocorr_HCR1"~ "Stationary (10y)",
                         "5yr_autocorr_HCR1" ~ "Stationary (5y)", 
                         "10yr_rwa_HCR1" ~ "Time-varying (10y)",    
                         "5yr_rwa_HCR1" ~ "Time-varying (5y)",     
                         "10yr_both_HCR1" ~ "Mixed (10y)" ,  
                         "5yr_both_HCR1" ~ "Mixed (5y)", 
                         "10yr_autocorr_HCR2" ~ "Stationary (10y)",
                         "5yr_autocorr_HCR2" ~ "Stationary (5y)", 
                         "10yr_rwa_HCR2" ~ "Time-varying (10y)",   
                         "5yr_rwa_HCR2" ~ "Time-varying (5y)",     
                         "10yr_both_HCR2" ~ "Mixed (10y)",  
                         "5yr_both_HCR2" ~ "Mixed (5y)",    
                         "10yr_autocorr_HCR3" ~ "Stationary (10y)",
                         "5yr_autocorr_HCR3" ~ "Stationary (5y)",
                         "10yr_rwa_HCR3" ~ "Time-varying (10y)",    
                         "5yr_rwa_HCR3" ~ "Time-varying (5y)",    
                         "10yr_both_HCR3" ~ "Mixed (10y)", 
                         "5yr_both_HCR3" ~ "Mixed (5y)" )



#visual assets
statusCols <- c("#9A3F3F","#DFD98D","#8EB687","gray75","gray25")
statusColsall <- c("#9A3F3F","gray90","gray10","#DFD98D","gray80","gray20","#8EB687","gray70","gray30")

#Figures: true status####
head(hcrdat3)
## assessed umsy & abundace benchmarks####
st_agg2<-ggplot(hcrdat2[hcrdat2$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")


st_agg3<-ggplot(hcrdat3[hcrdat3$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")


#stclass<-plot_grid(
#  st_agg, st_agg_all,
#  align = "h", axis = "bt",
#  rel_widths = c(.75,.25)
#)

#st_agg


st_agg_all<-ggplot(hcrdat_all[hcrdat_all$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('')+
  xlab('')+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



ggsave("figs_AFS/true_status_ERtar0.9UMSY_0.1BM.png",plot=st_agg_all,width=14,height=8.5)
#ggsave("outputs/figs/true_status_ERtar0.9UMSY_0.1BM.png",plot=stclass,width=14,height=8.5)





#Figures: true vs obs benchmarks - abundance scenarios####

#Fishery metric calcs - total catch, CV catch,####

##total catch and variance in annual catch - from 10y after burnin to end####

catch1_c=hcrdat_all[hcrdat_all$year>59&hcrdat_all$year<111,] %>% 
group_by(plotOM,plotMP,iteration) %>% 
summarize(m.catch=exp(mean(log(totalCatch)))) %>% 
srdat_all$plotOM<-dplyr::recode_factor(factor(srdat_all$nameOM),
                                     "stationarylAR1"  = "Static - 0.4AR1",
                                     "stationaryhAR1" = "Static - 0.8AR1",
                                     "decLinearcap0.50"="-50%Cap", 
                                     "decLinearcap0.25"="-75%Cap",
                                     "incLinearcap1.5"="+50%Cap",
                                     "incLinearProd2to3"="+50%Prod",
                                     "decLinearProd2to1.5"="-25%Prod",     
                                     "decLinearProd2to1"="-50%Prod",
                                     "decLinearProd2to0.5"="-75%Prod",
                                     "regProd2to1"='-50%PRegime',
                                     "regProd2to0.5"='-75%PRegime')


srdat_all$plotMP<-dplyr::recode_factor(factor(srdat_all$nameMP),
                                     "10yr_autocorr_adaptER" = "Stationary (10y)",
                                     "10yr_rwa_adaptER" = "Time-varying (10y)",
                                     "10yr_both_adaptER" = "Mixed (10y)", 
                                     "5yr_autocorr_adaptER" = "Stationary (5y)",
                                     "5yr_rwa_adaptER" = "Time-varying (5y)",
                                     "5yr_both_adaptER" = "Mixed (5y)")

group_by(plotOM) %>% summarize(max.catch=max(m.catch))

head(hcrdat_all)

#max catch
catch1_c=hcrdat_all[hcrdat_all$year>59&hcrdat_all$year<111,] %>% 
         group_by(plotOM,plotMP,HCR,iteration) %>% 
         summarize(m.catch=exp(mean(log(totalCatch)))) %>% 
         group_by(plotOM) %>% summarize(max.catch=max(m.catch))



catch_c=hcrdat_all[hcrdat_all$year>59&hcrdat_all$year<111,] %>% 
                group_by(plotOM,plotMP,HCR,iteration) %>% 
                summarize(total.catch=sum(log(totalCatch)),
                          m.catch=exp(mean(log(totalCatch))),
                          cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))),
                          aav.catch=median(abs(totalCatch[-1]-totalCatch[-length(totalCatch)])/
                            (totalCatch[-1]+totalCatch[-length(totalCatch)])))

#scaled catch
catch_c$scale.ann.catch=catch_c$m.catch/catch1_c$max.catch[match(catch_c$plotOM,catch1_c$plotOM)]



#proportion above benchmarks####
bmdat_c= hcrdat_all[hcrdat_all$year>59&hcrdat_all$year<111,] %>% 
group_by(plotOM,plotMP,HCR,iteration)  %>% 
summarize(p.red=sum(aboveLowerBM)/n(),
          p.bm=sum(aboveUpperBM)/n())
#aqui

bmdat_tot_c=hcrdat_all[hcrdat_all$year>59&hcrdat_all$year<111,] %>% 
group_by(plotOM,plotMP,HCR)  %>% 
summarize(p.red=sum(aboveLowerBM)/n(),
  p.bm=1-sum(aboveUpperBM)/n())

#ER series####

meanER_c=srdat_all[srdat_all$year<111,] %>% group_by(plotOM,plotMP,year,HCR) %>% 
reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_c$scnmp.t=paste(meanER_c$plotOM,meanER_c$plotMP,meanER_c$year,sep="_")
meanER_c=distinct(meanER_c,scnmp.t,.keep_all=T)


#spawners through time####
spawners_c1=srdat_all[srdat_all$year>59&srdat_all$year<111,] %>% group_by(plotOM,plotMP,year,HCR) %>% 
reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_c1$mpy=paste(spawners_c1$plotOM,spawners_c1$plotMP,spawners_c1$year)
spawners_c1=distinct(spawners_c1,mpy,.keep_all=TRUE)

#total spawners####
spawners_c2=srdat_all[srdat_all$year>59&srdat_all$year<111,] %>% 
  group_by(plotOM,plotMP,iteration,HCR) %>% 
  reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),
          smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_c2$mpy=paste(spawners_c2$plotOM,spawners_c2$plotMP,spawners_c2$iteration)
spawners_c2=distinct(spawners_c2,mpy,.keep_all=TRUE)


unique(spawners_c2$plotMP)

spawners_c2$plotMP<-dplyr::case_match(factor(spawners_c2$plotMP),
                         "10yr_autocorr_HCR1"~ "Stationary (10y)",
                         "5yr_autocorr_HCR1" ~ "Stationary (5y)", 
                         "10yr_rwa_HCR1" ~ "Time-varying (10y)",    
                         "5yr_rwa_HCR1" ~ "Time-varying (5y)",     
                         "10yr_both_HCR1" ~ "Mixed (10y)" ,  
                         "5yr_both_HCR1" ~ "Mixed (5y)", 
                         "10yr_autocorr_HCR2" ~ "Stationary (10y)",
                         "5yr_autocorr_HCR2" ~ "Stationary (5y)", 
                         "10yr_rwa_HCR2" ~ "Time-varying (10y)",   
                         "5yr_rwa_HCR2" ~ "Time-varying (5y)",     
                         "10yr_both_HCR2" ~ "Mixed (10y)",  
                         "5yr_both_HCR2" ~ "Mixed (5y)",    
                         "10yr_autocorr_HCR3" ~ "Stationary (10y)",
                         "5yr_autocorr_HCR3" ~ "Stationary (5y)",
                         "10yr_rwa_HCR3" ~ "Time-varying (10y)",    
                         "5yr_rwa_HCR3" ~ "Time-varying (5y)",    
                         "10yr_both_HCR3" ~ "Mixed (10y)", 
                         "5yr_both_HCR3" ~ "Mixed (5y)" )

#Figures: catch violin plots####

catch_plot_scc<-ggplot(catch_c,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  #scale_fill_manual(values=c('dodgerblue4','goldenrod4','dodgerblue1','goldenrod1', 'orangered4','orangered1')) +
  geom_violin()+
  #  geom_boxplot(width=0.1, color="white", alpha=0.5) +
  stat_summary(fun= median,
               geom = "crossbar", 
               width = 0.98,
               position = position_dodge(width = .1),
  )+
  facet_grid(plotOM~.)+
  ylab('Scaled Mean Annual Catch')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none")

catch_plot_scc
#ggsave("outputs/figs/catch_violins_ERtar0.9UMSY_0.1BM.png",plot=catch_plot_scc,width=8,height=12)


#Scenario subsets####

unique(hcrdat_all$plotMP)

#total catch and variance in annual catch
hcrdat_all.s=subset(hcrdat_all,plotOM %in%  c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))

hcrdat_all.s$plotOM=droplevels(hcrdat_all.s$plotOM)
hcrdat_all.s=subset(hcrdat_all.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)","Mixed (5y)"))

hcrdat_all.s$plotMP=dplyr::recode_factor(factor(hcrdat_all.s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying",
                                      "Mixed (5y)"="Mixed")
hcrdat_all.s$HCR

## assessed umsy & abundace benchmarks####
st_agg.s<-ggplot(hcrdat_all.s[hcrdat_all.s$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~HCR+plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
st_agg.s

st_agg_all.s<-ggplot(hcrdat_all.s[hcrdat_all.s$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~HCR)+
  ylab('')+
  xlab('')+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

stclass.s<-plot_grid(
  st_agg.s, st_agg_all.s,
  align = "h", axis = "bt",
  rel_widths = c(.75,.25)
)

stclass.s
ggsave("figs_HCR/status_stepwise.png",plot=stclass.s,width=14,height=8)



#subset of catch data

catch_c_s=subset(catch_c,plotOM %in% levels(hcrdat_all.s$plotOM))
catch_c_s=subset(catch_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)',"Mixed (5y)"))
catch_c_s$plotOM=droplevels(catch_c_s$plotOM)
catch_c_s$plotMP=dplyr::recode_factor(factor(catch_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying",
                                         "Mixed (5y)"="Mixed")





#subset of benchmark data

bmdat_c_s=subset(bmdat_c,plotOM %in% levels(hcrdat_all.s$plotOM))
bmdat_c_s=subset(bmdat_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)',"Mixed (5y)"))
bmdat_c_s$plotOM=droplevels(bmdat_c_s$plotOM)
bmdat_c_s$plotMP=dplyr::recode_factor(factor(bmdat_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying",
                                        "Mixed (5y)"="Mixed")


#subset of spawner data



spawners_c_s=subset(spawners_c2,plotOM %in% levels(hcrdat_all.s$plotOM))
spawners_c_s=subset(spawners_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)',"Mixed (5y)"))
spawners_c_s$plotOM=droplevels(spawners_c_s$plotOM)

spawners_c_s$plotMP=dplyr::recode_factor(factor(spawners_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying",
                                    "Mixed (5y)"="Mixed")



scn_cols_static=c('#045a8d',
  '#41b6c4',
  '#d95f0e',
  '#b30000','#58508d','#bc5090')

scn_cols_rwa=c('#fed98e',
  '#fe9929',
  '#d95f0e',
  '#993404')



## boxplots ##


allpm <- full_join(catch_c_s,bmdat_c_s)

#"m.catch",  "aav.catch",
allpm_long<-allpm |> 
     select("scale.ann.catch", "p.red","p.bm", "plotOM", "plotMP", "iteration","HCR" )|>
     tidyr::pivot_longer(!c("plotOM", "plotMP","iteration","HCR"))

#Figures: trade-off plots in catch vs. conservation risk (% below upper BM)#####


ggplot(allpm_long)+
geom_boxplot(aes(x=plotMP,y=value, fill=plotMP))+
theme_bw(15)+
facet_grid(name~plotOM,scales="free")
#



##trade-off plots
#allvars
allpmselec <-allpm |>  
     select("scale.ann.catch", "p.red", "plotOM", "plotMP", "iteration","HCR")|>
     tidyr::pivot_longer(!c("plotOM", "plotMP","iteration"))


all_comp_cat=expand.grid(mp=levels(hcrdat_all.s$plotMP),scn=levels(hcrdat_all.s$plotOM),
   hcr=levels(as.factor(hcrdat_all.s$HCR)),value=NA,l95=NA,u95=NA)

all_comp_p.red=expand.grid(mp=levels(hcrdat_all.s$plotMP),scn=levels(hcrdat_all.s$plotOM), 
  hcr=levels(as.factor(hcrdat_all.s$HCR)),value=NA,
  l95=NA,u95=NA)

all_comp_p.ub=expand.grid(mp=levels(hcrdat_all.s$plotMP),scn=levels(hcrdat_all.s$plotOM),
  hcr=levels(as.factor(hcrdat_all.s$HCR)),value=NA,l95=NA,u95=NA)

for(i in 1:length(levels(hcrdat_all.s$plotOM))){
  for(u in seq_along(unique(allpm$HCR))){
     
     ests=lm(scale.ann.catch~plotMP-1,data=subset(allpm,plotOM==levels(allpm$plotOM)[i]&HCR==unique(allpm$HCR)[u]))
     confs=confint(ests)

    catests<- allpm%>%subset(plotOM==levels(allpm$plotOM)[i]&HCR==unique(allpm$HCR)[u])%>%
           group_by(plotMP) %>%summarise(avg = median(scale.ann.catch), 
                                        up = quantile(scale.ann.catch,0.975), 
                                        low=quantile(scale.ann.catch,0.025))

 

     
  
  all_comp_cat[match(all_comp_cat$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
               match(all_comp_cat$hcr,unique(allpm$HCR)[u],nomatch =0)==1
                ,4]<-c(catests$avg)

  
  all_comp_cat[match(all_comp_cat$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
               match(all_comp_cat$hcr,unique(allpm$HCR)[u],nomatch =0)==1,5]=c(catests$low)
  
  all_comp_cat[match(all_comp_cat$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
               match(all_comp_cat$hcr,unique(allpm$HCR)[u],nomatch =0)==1,6]=c(catests$up)


    p.redests<- allpm%>%subset(plotOM==levels(allpm$plotOM)[i]&HCR==unique(allpm$HCR)[u])%>%
           group_by(plotMP) %>%summarise(avg = median(p.red), 
                                        up = quantile(p.red,0.975), 
                                        low=quantile(p.red,0.025))


  
  all_comp_p.red[match(all_comp_p.red$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
                 match(all_comp_cat$hcr,unique(allpm$HCR)[u],nomatch =0)==1,4]<-c(p.redests$avg)

  
  all_comp_p.red[match(all_comp_p.red$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
                 match(all_comp_p.red$hcr,unique(allpm$HCR)[u],nomatch =0)==1,5]=c(p.redests$low)


  all_comp_p.red[match(all_comp_p.red$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
                 match(all_comp_p.red$hcr,unique(allpm$HCR)[u],nomatch =0)==1,6]=c(p.redests$up)



p.ubests<- allpm%>%subset(plotOM==levels(allpm$plotOM)[i]&HCR==unique(allpm$HCR)[u])%>%
           group_by(plotMP) %>%summarise(avg = median(p.bm), 
                                        up = quantile(p.bm,0.975), 
                                        low=quantile(p.bm,0.025))
  
  all_comp_p.ub[match(all_comp_p.ub$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
               match(all_comp_p.ub$hcr,unique(allpm$HCR)[u],nomatch =0)==1,4]<-c(p.ubests$avg)

  
  all_comp_p.ub[match(all_comp_p.ub$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
                match(all_comp_p.ub$hcr,unique(allpm$HCR)[u],nomatch =0)==1,5]=c(p.ubests$low)
  all_comp_p.ub[match(all_comp_p.ub$scn,levels(hcrdat_all.s$plotOM)[i],nomatch =0)==1&
                match(all_comp_p.ub$hcr,unique(allpm$HCR)[u],nomatch =0)==1,6]=c(p.ubests$up)


  }

}

all_comp_p.red$variable<-"% above LB"
all_comp_p.ub$variable<-"% above upper benchmark"
all_comp_cat$variable<-"scaled catch"
all_comp_d<-rbind(all_comp_p.red,all_comp_cat)





lb_cat_tradeoff<-ggplot( all_comp_d, aes(y=value, x=variable, colour=mp,group=mp)) + 
    geom_errorbar(aes(ymin=l95, ymax=u95, colour=mp), width=.1) +
    geom_line(aes(y=as.numeric(value), x=variable, colour=mp)) +
    stat_summary(fun=max, geom="line",linewidth=2)+
    theme_bw(15) +
    scale_color_brewer(palette="Dark2") +
    geom_point( size=5)+
    facet_wrap(hcr~scn, scales="free")+
    guides(color=guide_legend(title="Reference point"))+
    theme(legend.position='bottom')


#rescale by scenario

for(i in seq_along(unique(all_comp_d$scn))){

  all_comp_d$value2<-all_comp_d$value
  ub_cat_comp$value2<-ub_cat_comp$value

  tdf<-all_comp_d|>subset(scn==unique(all_comp_d$scn)[i])
  tdfu<-ub_cat_comp|>subset(scn==unique(ub_cat_comp$scn)[i])

  tmp<-tdf$value2[tdf$variable=="% above LB"]
  tdf$value2[tdf$variable=="% above LB"]<-(tmp-min(tmp))/(max(tmp)-min(tmp))


  tmp2<-tdf$value2[tdf$variable=="scaled catch"]
  tdf$value2[tdf$variable=="scaled catch"]<-(tmp2-min(tmp2))/(max(tmp2)-min(tmp2)) *mean(tmp2)

  all_comp_d$value2[all_comp_d$scn==unique(all_comp_d$scn)[i]]<-tdf$value2


  tmp3<-tdfu$value2[tdfu$variable=="% above upper benchmark"]
  tdfu[tdfu=="% above upper benchmark"]<-(tmp3-min(tmp3))/(max(tmp3)-min(tmp3))


  tmp4<-tdfu$value2[tdf$variable=="scaled catch"]
  tdfu$value2[tdfu$variable=="scaled catch"]<-(tmp4-min(tmp4))/(max(tmp4)-min(tmp4))

  all_comp_d$value2[all_comp_d$scn==unique(all_comp_d$scn)[i]]<-tdf$value2
  ub_cat_comp$value2[ub_cat_comp$scn==unique(ub_cat_comp$scn)[i]]<-tdfu$value2

}



#ggsave("figs_HCR/tradeoff_lb_cat_stepwise.png",plot=lb_cat_tradeoff,width=16,height=5)


lb_cat_tradeoff_scaled<-ggplot( all_comp_d, aes(y=value, x=variable, colour=mp,group=mp)) + 
    #geom_errorbar(aes(ymin=l95, ymax=u95, colour=mp), width=.1) +
    #geom_line(aes(y=as.numeric(value), x=variable, colour=mp)) +
    stat_summary(fun=max, geom="line",linewidth=2)+
    theme_bw(15) +
    scale_color_brewer(palette="Dark2") +
    geom_point( size=5)+   
     facet_grid(hcr~scn)+
    #facet_wrap(hcr~scn, scales="free")+
    guides(color=guide_legend(title="Reference point"))+
    theme(legend.position='bottom')





ub_cat_comp<-rbind(all_comp_p.ub,all_comp_cat)

ub_cat_comp$value2<-ub_cat_comp$value


tmp3<-ub_cat_comp$value2[ub_cat_comp$variable=="% above upper benchmark"]
ub_cat_comp$value2[ub_cat_comp$variable=="% above upper benchmark"]<-(tmp3-min(tmp3))/(max(tmp3)-min(tmp3))


tmp4<-ub_cat_comp$value2[ub_cat_comp$variable=="scaled catch"]
ub_cat_comp$value2[ub_cat_comp$variable=="scaled catch"]<-(tmp4-min(tmp4))/(max(tmp4)-min(tmp4))



ub_cat_tradeoff<-ggplot( ub_cat_comp, aes(y=value, x=variable, colour=mp,group=mp)) + 
    #geom_errorbar(aes(ymin=l95, ymax=u95, colour=mp), width=.1) +
    #geom_line(aes(y=as.numeric(value), x=variable, colour=mp)) +
    stat_summary(fun=max, geom="line",linewidth=2)+
    theme_bw(20) +
    scale_color_brewer(palette="Dark2") +
    geom_point( size=5)+
    facet_grid(hcr~scn)+
    guides(color=guide_legend(title="Reference point"))+
    theme(legend.position='bottom')+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16))
ub_cat_tradeoff

ggsave("figs_HCR/tradeoff_ub_cat_stepwise20.png",plot=ub_cat_tradeoff,width=16,height=5)

ub_cat_comp$name<- c("moderate","liberal", "conservative" )[ub_cat_comp$hcr]

ub_cat_comp$name<- factor(ub_cat_comp$name, levels=c("moderate","liberal", "conservative" ))

NCRplot<-ggplot( ub_cat_comp |> subset(scn=="-50%Prod"), aes(y=value, x=variable, colour=mp,group=mp)) + 
    #geom_errorbar(aes(ymin=l95, ymax=u95, colour=mp), width=.1) +
    #geom_line(aes(y=as.numeric(value), x=variable, colour=mp)) +
    stat_summary(fun=max, geom="line",linewidth=2)+
    theme_bw(20) +
    scale_color_brewer(palette="Dark2") +
    geom_point( size=5)+
    facet_wrap(.~name, scales="free")+
    guides(color=guide_legend(title="Reference point"))+
    theme(legend.position='bottom')+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16))
ggsave("figs_HCR/tradeoff_minus50prod_3HCR.png",plot=NCRplot,width=16,height=5)


ub_cat_tradeoff2<-ggplot( ub_cat_comp, aes(y=value2, x=variable, colour=mp,group=mp)) + 
    #geom_errorbar(aes(ymin=l95, ymax=u95, colour=mp), width=.1) +
    #geom_line(aes(y=as.numeric(value), x=variable, colour=mp)) +
    stat_summary(fun=max, geom="line",linewidth=2)+
    theme_bw(20) +
    scale_color_brewer(palette="Dark2") +
    geom_point( size=5)+
    facet_grid(hcr~scn)+
    ylab("value")+
    guides(color=guide_legend(title="Reference point"))+
    theme(legend.position='bottom')+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 16))
ub_cat_tradeoff2


#catch
catch_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP),scn=levels(hcrdat3.s$plotOM),d.catch=NA,l95=NA,u95=NA,
                         m.catch=NA,e95=NA,e95=NA)

for(i in 1:length(levels(hcrdat3.s$plotOM))){

  ests=lm(m.catch~plotMP-1,data=subset(catch_c_s,plotOM==levels(catch_c_s$plotOM)[i]))
  confs=confint(ests)
  
  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,3]<-c(0,
                    ests$coefficients[2]/ests$coefficients[1]-1,
                    ests$coefficients[3]/ests$coefficients[1]-1)

  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,4]=c(0,
                                                                  confs[2,1]/confs[1,1]-1,
                                                                  confs[3,1]/confs[1,1]-1)
  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,5]=c(0,
                                                                   confs[2,2]/confs[1,2]-1,
                                                                   confs[3,2]/confs[1,2]-1)

  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,6]<-c(
                    ests$coefficients[1],
                    ests$coefficients[2],
                    ests$coefficients[3])

  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,7]=c(confs[1,1],
                                                                  confs[2,1],
                                                                  confs[3,1])
  catch_comp_d[match(catch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,8]=c(confs[1,2],
                                                                   confs[2,2],
                                                                   confs[3,2])


  

}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]


## aav catch ###

aavcatch_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP),scn=levels(hcrdat3.s$plotOM),d.aavcatch=NA,
                         l95=NA,u95=NA,m.aavcatch=NA,e95=NA,e95=NA)

for(i in 1:length(levels(hcrdat3.s$plotOM))){

  ests=lm(aav.catch~plotMP-1,data=subset(catch_c_s,plotOM==levels(catch_c_s$plotOM)[i]))
  confs=confint(ests)
  
  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,3]<-c(0,
                    ests$coefficients[2]/ests$coefficients[1]-1,
                    ests$coefficients[3]/ests$coefficients[1]-1)

  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,4]=c(0,
                                                                  confs[2,1]/confs[1,1]-1,
                                                                  confs[3,1]/confs[1,1]-1)
  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,5]=c(0,
                                                                   confs[2,2]/confs[1,2]-1,
                                                                   confs[3,2]/confs[1,2]-1)

  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,6]<-c(
                    ests$coefficients[1],
                    ests$coefficients[2],
                    ests$coefficients[3])

  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,7]=c(confs[1,1],
                                                                  confs[2,1],
                                                                  confs[3,1])
  aavcatch_comp_d[match(aavcatch_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,8]=c(confs[1,2],
                                                                   confs[2,2],
                                                                   confs[3,2])


  

}
aavcatch_comp_d$scn.cols=scn_cols_static[as.numeric(aavcatch_comp_d$scn)]


## no of spawners ####

spawners_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP),scn=levels(hcrdat3.s$plotOM),d.spwnrs=NA,l95=NA,u95=NA)
spawners_comp_d=as.data.frame(spawners_comp_d)
#  data.frame(scn=levels(hcrdat3.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat3.s$plotOM))){
  ests=lm(m.spawners~plotMP-1,data=subset(spawners_c_s,plotOM==levels(spawners_c_s$plotOM)[i]))
  confs=confint(ests)


  spawners_comp_d[match(spawners_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,3]=c(0,
                                                    ests$coefficients[2]/ests$coefficients[1]-1,
                                                    ests$coefficients[3]/ests$coefficients[1]-1)

  spawners_comp_d[match(spawners_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch=0)==1,4]=c(0,
                                                     confs[2,1]/confs[1,1]-1,
                                                     confs[3,1]/confs[1,1]-1)

  spawners_comp_d[match(spawners_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch=0)==1,5]=c(0,
                                                        confs[2,2]/confs[1,2]-1,
                                                        confs[3,2]/confs[1,2]-1)
}



## % of time above lbm  ####

status_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP),scn=levels(hcrdat3.s$plotOM),d.ubm=NA,m.ubm=NA,
                                                      l95=NA,u95=NA,d.red=NA,m.red=NA,l95r=NA,u95r=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat2.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)

for(i in 1:length(levels(hcrdat3.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_c_s,plotOM==levels(bmdat_c_s$plotOM)[i]))
  confs=confint(ests)
  
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,3]<-c(0,
                                                ests$coefficients[2]/ests$coefficients[1]-1,
                                                ests$coefficients[3]/ests$coefficients[1]-1)
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,4]<-c(
                                                ests$coefficients[1],
                                                ests$coefficients[2],
                                                ests$coefficients[3])
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,5]<-c(
                                                confs[1,1],
                                                confs[2,1],
                                                confs[3,1])
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,6]<-c(
                                                confs[1,2],
                                                confs[2,2],
                                                confs[3,2])

  ests2=lm(1-p.red~plotMP-1,data=subset(bmdat_c_s,plotOM==levels(bmdat_c_s$plotOM)[i]))
  confs2=confint(ests2)
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,7]<-c(0,
                                                ests2$coefficients[2]/ests2$coefficients[1]-1,
                                                ests2$coefficients[3]/ests2$coefficients[1]-1)
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,8]<-c(
                                                ests2$coefficients[1],
                                                ests2$coefficients[2],
                                                ests2$coefficients[3])
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,9]<-c(
                                                confs2[1,1],
                                                confs2[2,1],
                                                confs2[3,1])
  status_comp_d[match(status_comp_d$scn,levels(hcrdat3.s$plotOM)[i],nomatch =0)==1,10]<-c(
                                                confs2[1,2],
                                                confs2[2,2],
                                                confs2[3,2])
  
}




trade.off.df<-data.frame(mp=c(spawners_comp_d$mp,catch_comp_d$mp,status_comp_d$mp,status_comp_d$mp,
                              aavcatch_comp_d$mp),
                    scn=c(spawners_comp_d$scn,catch_comp_d$scn,status_comp_d$scn,status_comp_d$scn,
                          aavcatch_comp_d$scn), 
                    valuecomp= c(spawners_comp_d$d.spwnrs,catch_comp_d$d.catch,status_comp_d$d.ubm,
                                 status_comp_d$d.red,aavcatch_comp_d$d.aavcatch),
                    value= c(spawners_comp_d$d.spwnrs,catch_comp_d$d.catch,status_comp_d$m.ubm,
                             status_comp_d$m.red, aavcatch_comp_d$m.aavcatch ),
                    type=c(rep("spawners",nrow(spawners_comp_d)),rep("catch",nrow(catch_comp_d)),
                           rep("pbm",nrow(status_comp_d)),rep("pred",nrow(status_comp_d)),
                           rep("aav",nrow(aavcatch_comp_d))),
                    num=c(rep(1,nrow(spawners_comp_d)),rep(2,nrow(catch_comp_d)),
                           rep(3,nrow(status_comp_d)),rep(4,nrow(status_comp_d)),
                         rep(5,nrow(status_comp_d))))


ggplot(trade.off.df)+
geom_point(aes(x=mp,y=value, col=mp),size=5)+
#geom_line(aes(x=mp,y=valuecomp, col=mp),linewidth=2)+
theme_bw(15)+
facet_grid(type~scn,scales="free")



png(filename='outputs/figs/catch_conservation_tradeoff_UBM_ERtar0.9UMSY_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-20,50)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-20,50))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=50,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=60,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=70,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(spawners_comp_d$d.spwnrs[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  #lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  #lines(c(spawners_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(catch_comp_d$d.catch[i],cex=2,pch=21)
  points(spawners_comp_d$d.uspwnrs[i]*100~1.25+j,data=spawners_comp_d,cex=2,bg=catch_comp_d$scn[i],pch=21)
 }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above UBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_LBM_ERtar0.9UMSY_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-20,70)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-20,50))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=50,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=60,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=70,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(spawners_comp_d$d.red[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95r[i]*100,status_comp_d$u95r[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.red[i]*100~1.25+j,data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
}
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above LBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_spawner_ERtar0.9UMSY_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-20,70)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-20,50))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=50,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=60,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=70,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(spawners_comp_d$d.spawners[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.spwnrs[i]*100~1.25+j,data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
}
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Mean annual spawners',line=1.5)
dev.off()





#Figures: composite plots - status/total status/catch####

#@B - none; F - umsy####
#subset of scenarios:

#@B - assessed; F - assessed####
st_agg<-ggplot(hcrdat3[hcrdat3$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(legend.position = "none", strip.background = element_blank(),
      strip.text.y=element_blank())

st_agg_all<-ggplot(hcrdat3[hcrdat3$year>59&hcrdat3$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('Proportion of years in status')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text.x = element_blank(),strip.text.y=element_blank())


catch_plot_scu<-ggplot(catch_c,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  scale_fill_manual(values=c('blue4','darkgoldenrod','dodgerblue4','goldenrod2')) +
  geom_violin()+
  scale_x_discrete(position = "top") +
  #  geom_boxplot(width=0.1, color="white", alpha=0.5) +
  stat_summary(fun= median,
               geom = "crossbar", 
               width = 0.98,
               position = position_dodge(width = .1),
  )+
  facet_grid(plotOM~.)+
  ylab('Scaled Mean Annual Catch')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none",strip.background = element_blank(),strip.text = element_text(face = "bold"))


stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  ncol=3,greedy=F,align='hv',axis='bt',
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("./outputs/figs/composite_plot_status_catch_ERtar0.9UMSY_0.1BM.png",plot=stclass,width=14,height=8.5)


#@B - fixed; F - assessed####


#Figure: ER trends####

#NO B/assessed F ref pts

#combined - assessed B/F ref pts
er_plot_bmumsy<-ggplot(meanER_c)+
  geom_line(aes(x=year-50,y=m.ER,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_line(aes(x=year-50,y=umsy,group=plotOM),linewidth=1.2,colour='green')+
  geom_ribbon(aes(ymin = l80.er, ymax = u80.er,x=year-50), alpha = 0.2)+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Exploitation Rate')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

er_plot_bmumsy
ggsave("./outputs/figs/er_comp_ERtar0.9UMSY_0.1BM.png",plot=er_plot_bmumsy,width=8,height=6)



#combined - assessed B/F ref pts
sp_plot_bmumsy<-ggplot(spawners_c1)+
  geom_line(aes(x=year-50,y=m.spawners,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_ribbon(aes(ymin = l80.s, ymax = u80.s,x=year-50), alpha = 0.2)+
  geom_line(aes(x=year-50,y=smsy,group=plotOM),linewidth=1.2,colour='green')+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Spawner abundance')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot_bmumsy
ggsave("./outputs/figs/spawners_timeseries_vs_smsy_ERtar0.9UMSY_0.1BM.png",plot=sp_plot_bmumsy,width=8,height=6)


sp_plot2_bmumsy<-ggplot(spawners_c2,aes(x=plotMP,y=log10(m.spawners),fill=plotMP))+
  coord_trans(y = "log10")+
  scale_fill_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  geom_violin()+
  
  stat_summary(fun= median,
               geom = "crossbar", 
               width = 0.98,
               position = position_dodge(width = .1),
  )+
  # geom_boxplot(width=0.1, color="white", alpha=0.5) +
  facet_grid(plotOM~.)+
  ylab('log10[Mean Spawner Abundance]')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot2_bmumsy
ggsave("outputs/figs/spawners_violin_ERtar0.9UMSY_0.1BM.png",plot=sp_plot2_bmumsy,width=8,height=12)
