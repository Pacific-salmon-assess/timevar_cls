#install samsim
#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)

library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)
library(cowplot)

source("R/util_funcs.R")



#gudelines scenarios - load data####
#simPars_um <- read.csv("data/guidelines/SimPars2.0.csv") #simpars for ER tracking umsy, no EG
cuPar <- read.csv("data/guidelines/CUPars2.0.csv") #cu pars
#simPars_b <- read.csv("data/guidelines/SimPars2.1.csv") #simpars for fixed ER, assessed EG
simPars_c <- read.csv("data/guidelines/SimPars2.2.csv") #simpars for ER tracking umsy, assessed EG
#simPars_d <- read.csv("data/guidelines/SimPars2.3.csv") #simpars for ER tracking umsy, fixed EG
#here()

#processing data####
hcrDatalist_c<-list()
srData_c<- list()

for(a in seq_len(nrow(simPars_c))){
  
  hcrDatalist_c[[a]] <- readRDS(paste0("./clusteresult/SamSimOutputs/simData/",
                                       simPars_c$nameOM[a],"/", 
                                       simPars_c$scenario[a],"/",
                                       paste(simPars_c$nameOM[a],"_", simPars_c$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist_c[[a]]$scenario<-simPars_c$scenario[a]
  hcrDatalist_c[[a]]$nameOM<-simPars_c$nameOM[a]
  hcrDatalist_c[[a]]$nameMP<-simPars_c$nameMP[a]
  
  
  srData_c[[a]] <- readRDS(paste0("./clusteresult/SamSimOutputs/simData/", 
                                  simPars_c$nameOM[a],"/",
                                  simPars_c$scenario[a],"/",
                                  paste(simPars_c$nameOM[a],"_", simPars_c$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  
  srData_c[[a]]$scenario<-simPars_c$scenario[a]
  srData_c[[a]]$nameOM<-simPars_c$nameOM[a]
  srData_c[[a]]$nameMP<-simPars_c$nameMP[a]
}

hcrdat_c <- do.call(rbind,hcrDatalist_c)
srdat_c<- do.call(rbind,srData_c)

hcrdat_c<-hcrdat_c[hcrdat_c$year>50,]
srdat_c<-srdat_c[srdat_c$year>50,]



scnwboth<-unique(hcrdat_c$scenario)[grep(pattern="both", unique(hcrdat_c$scenario))]
nameMPboth<-gsub("rwa", "both", hcrdat_c$nameMP[hcrdat_c$scenario%in%scnwboth])
hcrdat_c$nameMP[hcrdat_c$scenario%in%scnwboth]<-nameMPboth

unique(hcrdat_c$nameMP)
hcrdat3=format_dat(hcrdat_c)

names(hcrdat3)

unique(hcrdat3$plotMP)

scnwboth<-unique(hcrdat3$scenario)[grep(pattern="both", unique(hcrdat3$scenario))]
nameMPboth<-gsub("Time-varying", "Mixed", hcrdat3$nameMP[hcrdat3$scenario%in%scnwboth])
hcrdat3$nameMP[hcrdat3$scenario%in%scnwboth]<-nameMPboth

plotMPboth<-gsub("Time-varying", "Mixed", hcrdat3$plotMP[hcrdat3$scenario%in%scnwboth])
hcrdat3$plotMP[hcrdat3$scenario%in%scnwboth]<-plotMPboth
#recode scenario names####



scnwboth<-unique(srdat_c$scenario)[grep(pattern="both", unique(srdat_c$scenario))]

srdat_c$nameMP[srdat_c$scenario%in%scnwboth]
nameMPboth<-gsub("rwa", "both", srdat_c$nameMP[srdat_c$scenario%in%scnwboth])
srdat_c$nameMP[srdat_c$scenario%in%scnwboth]<-nameMPboth

unique(hcrdat3$plotMP)

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


#visual assets
statusCols <- c("#9A3F3F","#DFD98D","#8EB687","gray75","gray25")
statusColsall <- c("#9A3F3F","gray90","gray10","#DFD98D","gray80","gray20","#8EB687","gray70","gray30")

#Figures: true status####

## assessed umsy & abundace benchmarks####
st_agg<-ggplot(hcrdat3[hcrdat3$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
#st_agg
st_agg_all<-ggplot(hcrdat3[hcrdat3$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('')+
  xlab('')+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

stclass<-plot_grid(
  st_agg, st_agg_all,
  align = "h", axis = "bt",
  rel_widths = c(.75,.25)
)

stclass
#ggsave("outputs/figs/true_status_ERtar0.9UMSY_0.1BM.png",plot=stclass,width=14,height=8.5)

#Figures: Guidance doc outputs with a limited set of scenarios####


unique(srdat_c$nameMP)


#Figures: true vs obs benchmarks - abundance scenarios####

#Fishery metric calcs - total catch, CV catch,####

srdat_c$plotOM<-dplyr::recode_factor(factor(srdat_c$nameOM),
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


srdat_c$plotMP<-dplyr::recode_factor(factor(srdat_c$nameMP),
                                     "10yr_autocorr_adaptER" = "Stationary (10y)",
                                     "10yr_rwa_adaptER" = "Time-varying (10y)",
                                     "10yr_both_adaptER" = "Mixed (10y)", 
                                     "5yr_autocorr_adaptER" = "Stationary (5y)",
                                     "5yr_rwa_adaptER" = "Time-varying (5y)",
                                     "5yr_both_adaptER" = "Mixed (5y)")

#total catch and variance in annual catch - from 10y after burnin to end####

catch1_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_c$scale.ann.catch=catch_c$m.catch/catch1_c$max.catch[match(catch_c$plotOM,catch1_c$plotOM)]


#proportion above benchmarks####
bmdat_c= hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())

bmdat_tot_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())

#ER series####

meanER_c=srdat_c[srdat_c$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_c$scnmp.t=paste(meanER_c$plotOM,meanER_c$plotMP,meanER_c$year,sep="_")
meanER_c=distinct(meanER_c,scnmp.t,.keep_all=T)


#spawners through time####
spawners_c1=srdat_c[srdat_c$year>59&srdat_c$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_c1$mpy=paste(spawners_c1$plotOM,spawners_c1$plotMP,spawners_c1$year)
spawners_c1=distinct(spawners_c1,mpy,.keep_all=TRUE)

#total spawners####
spawners_c2=srdat_c[srdat_c$year>59&srdat_c$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_c2$mpy=paste(spawners_c2$plotOM,spawners_c2$plotMP,spawners_c2$iteration)
spawners_c2=distinct(spawners_c2,mpy,.keep_all=TRUE)

  

#Figures: catch violin plots####

catch_plot_scc<-ggplot(catch_c,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  scale_fill_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
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
ggsave("outputs/figs/catch_violins_ERtar0.9UMSY_0.1BM.png",plot=catch_plot_scc,width=8,height=12)


#Scenario subsets####

#total catch and variance in annual catch
hcrdat3.s=subset(hcrdat3,plotOM %in%  c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))
hcrdat3.s$plotOM=droplevels(hcrdat3.s$plotOM)
hcrdat3.s=subset(hcrdat3.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)"))
hcrdat3.s$plotMP=droplevels(hcrdat3.s$plotMP)
hcrdat3.s$plotMP=dplyr::recode_factor(factor(hcrdat3.s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")


#subset of catch data

catch_c_s=subset(catch_c,plotOM %in% levels(hcrdat3.s$plotOM))
catch_c_s=subset(catch_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
catch_c_s$plotOM=droplevels(catch_c_s$plotOM)
catch_c_s$plotOM=droplevels(catch_c_s$plotOM)
catch_c_s$plotMP=dplyr::recode_factor(factor(catch_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")


#subset of benchmark data

bmdat_c_s=subset(bmdat_c,plotOM %in% levels(hcrdat3.s$plotOM))
bmdat_c_s=subset(bmdat_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
bmdat_c_s$plotOM=droplevels(bmdat_c_s$plotOM)
bmdat_c_s$plotOM=droplevels(bmdat_c_s$plotOM)
bmdat_c_s$plotMP=dplyr::recode_factor(factor(bmdat_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")


#subset of spawner data

spawners_c_s=subset(spawners_c2,plotOM %in% levels(hcrdat3.s$plotOM))
spawners_c_s=subset(spawners_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
spawners_c_s$plotOM=droplevels(spawners_c_s$plotOM)
spawners_c_s$plotOM=droplevels(spawners_c_s$plotOM)
spawners_c_s$plotMP=dplyr::recode_factor(factor(spawners_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")


scn_cols_static=c('#045a8d',
  '#41b6c4',
  '#d95f0e',
  '#b30000','#58508d','#bc5090')

scn_cols_rwa=c('#fed98e',
  '#fe9929',
  '#d95f0e',
  '#993404')


#Figures: trade-off plots in catch vs. conservation risk (% below upper BM)#####
##UMSY track - no EG####



## Umsy and bm####

catch_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP)[-1],scn=levels(hcrdat3.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat3.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat3.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_c_s,plotOM==levels(catch_c_s$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat3.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat3.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat3.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  
}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP)[-1],scn=levels(hcrdat3.s$plotOM),d.ubm=NA,l95=NA,u95=NA,d.red=NA,l95r=NA,u95r=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat3.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat3.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_c_s,plotOM==levels(bmdat_c_s$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  ests2=lm(1-p.red~plotMP-1,data=subset(bmdat_c,plotOM==levels(bmdat_c$plotOM)[i]))
  confs2=confint(ests2)
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),6]=ests2$coefficients[2]/ests2$coefficients[1]-1
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),7]=confs2[2,1]/confs2[1,2]-1
  status_comp_d[match(levels(hcrdat3.s$plotOM)[i],status_comp_d$scn),8]=confs2[2,2]/confs2[1,1]-1
}

spawners_comp_d=expand.grid(mp=levels(hcrdat3.s$plotMP)[-1],scn=levels(hcrdat3.s$plotOM),d.spwnrs=NA,l95=NA,u95=NA)
spawners_comp_d=as.data.frame(spawners_comp_d)
#  data.frame(scn=levels(hcrdat3.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat3.s$plotOM))){
  ests=lm(m.spawners~plotMP-1,data=subset(spawners_c_s,plotOM==levels(spawners_c_s$plotOM)[i]))
  confs=confint(ests)
  spawners_comp_d[match(levels(hcrdat3.s$plotOM)[i],spawners_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  spawners_comp_d[match(levels(hcrdat3.s$plotOM)[i],spawners_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  spawners_comp_d[match(levels(hcrdat3.s$plotOM)[i],spawners_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
}


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
  lines(c(spawners_comp_d$d.ubm[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.ubm[i]*100~1.25+j,data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
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
