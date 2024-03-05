#=================================================
#run closed loop simulations
#=================================================

#setwd(paste0("C:\\Users\\worc\\Documents\\timevar\\samSim\\R","/.."))
#devtools::document()
#devtools::load_all()

#install samsim
#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)


library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)
library(cowplot)

source("R/util_funcs.R")



#gudelines scenarios example
simPars_um <- read.csv("data/guidelines/SimPars2.0.csv")
cuPar_um <- read.csv("data/guidelines/CUPars2.0.csv")
#here()

hcrDatalist_um<-list()
srData_um<- list()

for(a in seq_len(nrow(simPars_um))){

   hcrDatalist_um[[a]] <- readRDS(paste0("umsy_track/SamSimOutputs/simData/",
   simPars_um$nameOM[a],"/", 
   simPars_um$scenario[a],"/",
   paste(simPars_um$nameOM[a],"_", simPars_um$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist_um[[a]]$scenario<-simPars_um$scenario[a]
   hcrDatalist_um[[a]]$nameOM<-simPars_um$nameOM[a]
   hcrDatalist_um[[a]]$nameMP<-simPars_um$nameMP[a]


   srData_um[[a]] <- readRDS(paste0("umsy_track/SamSimOutputs/simData/", 
                         simPars_um$nameOM[a],"/",
                         simPars_um$scenario[a],"/",
                         paste(simPars_um$nameOM[a],"_", simPars_um$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout


  srData_um[[a]]$scenario<-simPars_um$scenario[a]
  srData_um[[a]]$nameOM<-simPars_um$nameOM[a]
  srData_um[[a]]$nameMP<-simPars_um$nameMP[a]
}

hcrdat_um <- do.call(rbind,hcrDatalist_um)
srdat_um<- do.call(rbind,srData_um)

hcrdat_um<-hcrdat_um[hcrdat_um$year>50,]
srdat_um<-srdat_um[srdat_um$year>50,]

hcrdat1=format_dat(hcrdat_um)

simPars_b <- read.csv("data/guidelines/SimPars2.1.csv")

hcrDatalist_b<-list()
srData_b<- list()

for(a in seq_len(nrow(simPars_b))){
  
  hcrDatalist_b[[a]] <- readRDS(paste0("bm_track/SamSimOutputs/simData/",
                                        simPars_b$nameOM[a],"/", 
                                        simPars_b$scenario[a],"/",
                                        paste(simPars_b$nameOM[a],"_", simPars_b$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist_b[[a]]$scenario<-simPars_b$scenario[a]
  hcrDatalist_b[[a]]$nameOM<-simPars_b$nameOM[a]
  hcrDatalist_b[[a]]$nameMP<-simPars_b$nameMP[a]
  
  
  srData_b[[a]] <- readRDS(paste0("bm_track/SamSimOutputs/simData/", 
                                   simPars_b$nameOM[a],"/",
                                   simPars_b$scenario[a],"/",
                                   paste(simPars_b$nameOM[a],"_", simPars_b$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  
  srData_b[[a]]$scenario<-simPars_b$scenario[a]
  srData_b[[a]]$nameOM<-simPars_b$nameOM[a]
  srData_b[[a]]$nameMP<-simPars_b$nameMP[a]
}

hcrdat_b <- do.call(rbind,hcrDatalist_b)
srdat_b<- do.call(rbind,srData_b)

hcrdat_b<-hcrdat_b[hcrdat_b$year>50,]
srdat_b<-srdat_b[srdat_b$year>50,]

hcrdat2=format_dat(hcrdat_b)

##true status - ER feedback scenarios

statusCols <- c("#9A3F3F","#DFD98D","#8EB687","gray75","gray25")
statusColsall <- c("#9A3F3F","gray90","gray10","#DFD98D","gray80","gray20","#8EB687","gray70","gray30")

#umsy benchmark
st_agg<-ggplot(hcrdat1[hcrdat1$year<111,])+
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
st_agg_all<-ggplot(hcrdat1[hcrdat1$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('')+
  xlab('')+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),)

stclass<-plot_grid(
  st_agg, st_agg_all,
  align = "h", axis = "bt",
  rel_widths = c(.75,.25)
)

stclass
ggsave("figs/true_status_ERtar_0.9UMSY.png",plot=stclass,width=14,height=8.5)

#abundance benchmark
st_agg<-ggplot(hcrdat2[hcrdat2$year<111,])+
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
st_agg_all<-ggplot(hcrdat2[hcrdat2$year<111,])+
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
ggsave("figs/true_status_ERtar0.65_0.1BM.png",plot=stclass,width=14,height=8.5)

##estimated status - ER feedback scenarios
st_agg<-ggplot(hcrdat2[hcrdat2$year<111,])+
  geom_bar(aes(x=year-50, fill=status_agg),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
#st_agg
st_agg_all<-ggplot(hcrdat2[hcrdat2$year<111,])+
  geom_bar(aes(x=plotMP, fill=status_agg),position = "fill")+
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
ggsave("figs/status_classification_ERtar0.65_0.1BM.png",plot=stclass,width=14,height=8.5)


#==================================================
#true vs obs benchmarks and umsy


BMcomp2<-data.frame(year=hcrdat2$year,
  iteration=hcrdat2$iteration,
  typeest=as.factor(rep(c("truth","obs"),each=length(hcrdat2$year))),
  lowerBM=c(hcrdat2$lowerBM,hcrdat2$lowerObsBM),
  upperBM=c(hcrdat2$upperBM,hcrdat2$upperObsBM),
  nameOM=hcrdat2$nameOM,
  nameMP=hcrdat2$nameMP,
  plotOM=hcrdat2$plotOM
  )


upperBMtrue_hcr<-ggplot(BMcomp2)+
geom_line(aes(x=year-50,y=upperBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d(begin=.1, end=.8) +
facet_grid(plotOM~nameMP)+
theme_bw(14)
upperBMtrue_hcr
ggsave("figs/upperBM_true_hcr_ERtar0.65_0.1BM.png",plot=upperBMtrue_hcr,width=14,height=8.5)

lowerBMtrue_hcr<-ggplot(BMcomp2)+
geom_line(aes(x=year-50,y=lowerBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d(begin=.1, end=.8) +
facet_grid(plotOM~nameMP)+
coord_cartesian(ylim=c(0,70000))+
theme_bw(14)
lowerBMtrue_hcr
ggsave("figs/lowerBM_true_hcr_ERtar0.65_0.1BM.png",plot=lowerBMtrue_hcr,width=14,height=8.5)

BMcompt=subset(BMcomp,typeest=='truth')

bm_plot1<-ggplot(BMcompt)+
  geom_line(aes(x=year-50,y=upperBM,color=plotOM,group=plotOM),linewidth=1.2)+
  scale_color_manual(values=c('#a1dab4','#a1dab4','#d95f0e','#fed98e','#253494','#41b6c4','darkred'),name='') +
  ylab('Upper benchmark: 80% Smsy')+
  xlab(
    ''
  )+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))

bm_plot2<-ggplot(BMcompt)+
  geom_line(aes(x=year-50,y=lowerBM,color=plotOM,group=plotOM),linewidth=1.2)+
  scale_color_manual(values=c('#a1dab4','#a1dab4','#d95f0e','#fed98e','#253494','#41b6c4','darkred')) +
  ylab('Lower benchmark: Sgen')+
  xlab(
    'Year of simulation'
  )+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))

legend = cowplot::get_legend(bm_plot1)
bm_trends1<-plot_grid(
  bm_plot1+ theme(legend.position="none"),
bm_plot2+ theme(legend.position="none"),
ncol=1,nrow=2
)
bm_trends2=plot_grid(bm_trends1,legend,rel_widths=c(.75,.25))
ggsave("figs/trueBM_comp.png",plot=bm_trends2,width=8,height=6)


##Fishery Risk - total catch, CV catch,
srdat_um$plotOM<-dplyr::recode_factor(factor(srdat_um$nameOM),"stationaryhAR1" = "AR1 high",
                                    "stationarylAR1" = "AR1 low",  
                                    "decLinearcap0.5"="cap * 0.5", 
                                    "decLinearcap0.25"="cap * 0.25",
                                    "incLinearcap1.5"='cap * 1.5',
                                    "decLinearProd2to1.5"="prod 2 -> 1.5",     
                                    "decLinearProd2to1"="prod 2 -> 1",
                                    "decLinearProd2to0.5"="prod 2 -> 0.5",
                                    "incLinearProd2to3"="prod 2 -> 3")

srdat_um$plotMP<-dplyr::recode_factor(factor(srdat_um$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                    "10yr_rwa_constER" = "Time-varying (10y)",  
                                    "5yr_autocorr_constER" = "Stationary (5y)",
                                    "5yr_rwa_constER" = "Time-varying (5y)",
                                    "10yr_autocorr_adaptER" = "Stationary (10y)",
                                    "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                    "5yr_autocorr_adaptER" = "Stationary (5y)",
                                    "5yr_rwa_adaptER" = "Time-varying (5y)")

srdat_b$plotOM<-dplyr::recode_factor(factor(srdat_um$nameOM),"stationaryhAR1" = "AR1 high",
                                      "stationarylAR1" = "AR1 low",  
                                      "decLinearcap0.5"="cap * 0.5", 
                                      "decLinearcap0.25"="cap * 0.25",
                                      "incLinearcap1.5"='cap * 1.5',
                                      "decLinearProd2to1.5"="prod 2 -> 1.5",     
                                      "decLinearProd2to1"="prod 2 -> 1",
                                      "decLinearProd2to0.5"="prod 2 -> 0.5",
                                      "incLinearProd2to3"="prod 2 -> 3")

srdat_b$plotMP<-dplyr::recode_factor(factor(srdat_um$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                      "10yr_rwa_constER" = "Time-varying (10y)",  
                                      "5yr_autocorr_constER" = "Stationary (5y)",
                                      "5yr_rwa_constER" = "Time-varying (5y)",
                                      "10yr_autocorr_adaptER" = "Stationary (10y)",
                                      "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                      "5yr_autocorr_adaptER" = "Stationary (5y)",
                                      "5yr_rwa_adaptER" = "Time-varying (5y)")



#total catch and variance in annual catch
catch1_u=hcrdat1 %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_u=hcrdat1 %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(totalCatch),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_u$scale.ann.catch=catch_u$m.catch/catch1_u$max.catch[match(catch_u$plotOM,catch1_u$plotOM)]


catch1_b=hcrdat2 %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_b=hcrdat2 %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(totalCatch),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_b$scale.ann.catch=catch_b$m.catch/catch1_b$max.catch[match(catch_b$plotOM,catch1_b$plotOM)]

#proportion above benchmarks
bmdat_u= hcrdat1 %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_b= hcrdat2 %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())

#ER series
meanER_u=srdat_um %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_u$scnmp.t=paste(meanER_u$plotOM,meanER_u$plotMP,meanER_u$year,sep="_")
meanER_u=distinct(meanER_u,scnmp.t,.keep_all=T)
meanER_b=srdat_b %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_b$scnmp.t=paste(meanER_b$plotOM,meanER_b$plotMP,meanER_b$year,sep="_")
meanER_b=distinct(meanER_b,scnmp.t,.keep_all=T)

#total spawners
spawners_u=srdat_um %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_u$mpy=paste(spawners_u$plotOM,spawners_u$plotMP,spawners_u$year)
spawners_u=distinct(spawners_u,mpy,.keep_all=TRUE)
spawners_b=srdat_b %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_b$mpy=paste(spawners_b$plotOM,spawners_b$plotMP,spawners_b$year)
spawners_b=distinct(spawners_b,mpy,.keep_all=TRUE)

#green status spawners
spawners_u_g=srdat_um[hcrdat1$status=='green',] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_u_g$mpy=paste(spawners_u_g$plotOM,spawners_u_g$plotMP,spawners_u_g$year)
spawners_u_g=distinct(spawners_u_g,mpy,.keep_all=TRUE)
spawners_b_g=srdat_b[hcrdat2$status=='green',] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_b_g$mpy=paste(spawners_b_g$plotOM,spawners_b_g$plotMP,spawners_b_g$year)
spawners_b_g=distinct(spawners_b_g,mpy,.keep_all=TRUE)


  

#comparison of catch:
catch_plot_scu<-ggplot(catch_u,aes(x=plotMP,y=cv.catch,fill=plotMP))+
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

catch_plot_scu
ggsave("figs/catch_violins_ERtar_0.9Umsy.png",plot=catch_plot_scu,width=8,height=12)

catch_plot_scb<-ggplot(catch_b,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
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

catch_plot_scb
ggsave("figs/catch_violins_ERtar0.65_0.1BM.png",plot=catch_plot_scb,width=8,height=12)


#Summary plot: proportional change in catch vs. conservation risk (% below upper BM)
#scenarios
hcrdat1.s=subset(hcrdat1,plotOM %in% c('AR1 high','cap * 0.25','prod 2 -> 1','prod 2 -> 0.5'))
hcrdat1.s$plotOM=droplevels(hcrdat1.s$plotOM)
hcrdat1.s=subset(hcrdat1.s,plotMP %in% c("Stationary (10y)" ,"Time-varying (10y)"))
hcrdat1.s$plotMP=droplevels(hcrdat1.s$plotMP)

#total catch and variance in annual catch
catch1_u=hcrdat1.s %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_u=hcrdat1.s %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(totalCatch),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_u$scale.ann.catch=catch_u$m.catch/catch1_u$max.catch[match(catch_u$plotOM,catch1_u$plotOM)]

hcrdat2.s=subset(hcrdat2,plotOM %in% c('AR1 high','cap * 0.25','prod 2 -> 1','prod 2 -> 0.5'))
hcrdat2.s$plotOM=droplevels(hcrdat2.s$plotOM)
hcrdat2.s=subset(hcrdat2.s,plotMP %in% c("Stationary (10y)" ,"Time-varying (10y)"))
hcrdat2.s$plotMP=droplevels(hcrdat2.s$plotMP)

#total catch and variance in annual catch
catch2_u=hcrdat2.s %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_u2=hcrdat2.s %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(totalCatch),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_u2$scale.ann.catch=catch_u2$m.catch/catch2_u$max.catch[match(catch_u2$plotOM,catch2_u$plotOM)]


##To do: finish summary plot by just comparing 10-y static to 10-y rwa

scn_cols_static=c('#bdd7e7',
  '#6baed6',
  '#3182bd',
  '#08519c')

scn_cols_rwa=c('#fed98e',
  '#fe9929',
  '#d95f0e',
  '#993404')

catch_comp_d=expand.grid(mp=levels(hcrdat1.s$plotMP)[-1],scn=levels(hcrdat1.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat1.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat1.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_u,plotOM==levels(catch_u$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1

}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat1.s$plotMP)[-1],scn=levels(hcrdat1.s$plotOM),d.status=NA,l95=NA,u95=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat1.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat1.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_u,plotOM==levels(bmdat_u$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
}

png(filename='figs/catch_conservation_tradeoff_ERtar_0.9Umsy.png',width=8,height=6,units='in',res=600)
plot(c(-50,50)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-40,40))
abline(h=0)
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  lines(c(catch_comp_d$d.catch[i]*100,status_comp_d$d.status[i]*100)~c(1.25,1.75),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.25,2))
  lines(c(status_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.75,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
}
points(d.catch*100~rep(1.25,nrow(catch_comp_d)),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols,pch=21)
points(d.status*100~rep(1.75,nrow(status_comp_d)),data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols,pch=21)
mtext(side=1,at=1.25,'Mean Annual Catch',line=1.5)
mtext(side=1,at=1.75,'Proportion of time above EG/UBM',line=1.5)
dev.off()


catch_comp_d=expand.grid(mp=levels(hcrdat2.s$plotMP)[-1],scn=levels(hcrdat2.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat2.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat2.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_u2,plotOM==levels(catch_u2$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  
}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat2.s$plotMP)[-1],scn=levels(hcrdat2.s$plotOM),d.status=NA,l95=NA,u95=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat2.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat2.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_b,plotOM==levels(bmdat_b$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
}

png(filename='figs/catch_conservation_tradeoff_ERtar0.65_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-50,50)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-40,40))
abline(h=0)
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  lines(c(catch_comp_d$d.catch[i]*100,status_comp_d$d.status[i]*100)~c(1.25,1.75),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.25,2))
  lines(c(status_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.75,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
}
points(d.catch*100~rep(1.25,nrow(catch_comp_d)),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols,pch=21)
points(d.status*100~rep(1.75,nrow(status_comp_d)),data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols,pch=21)
mtext(side=1,at=1.25,'Mean Annual Catch',line=1.5)
mtext(side=1,at=1.75,'Proportion of time above EG/UBM',line=1.5)
dev.off()


#copy for 0.1BM


##ggplot wont plot error bars - making this effort USELESS and wasting my damn time. ggplot is trash.

#catch_plot2<-ggplot(as.data.frame(catch_comp_d),aes(x=factor(mp),y=d.catch*100))+
 # geom_abline(aes(intercept=0,slope=0))+
#  geom_bar(stat="identity", colour="black",position="dodge") + 
 # geom_errorbar(aes(ymin=l95, ymax =u95), width = 0.5,position=position_dodge(.9))+
#  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
#  geom_point(size=3)+
#  geom_pointrange(aes(ymin = catch_comp_d$l95*100, ymax = catch_comp_d$u95*100,x=factor(mp)))+
#  facet_grid(scn~.)+
 # ylab('% Difference in mean annual catch relative to 10y stationary')+
 # xlab('')+
 # theme_bw(12)+
 # theme(legend.position="none")

#catch_plot2
#ggsave("figs/catch_violins_ERtar_0.9Umsy.png",plot=catch_plot,width=8,height=6)

#er trends
er_plot_umsy<-ggplot(meanER_u)+
  geom_line(aes(x=year-50,y=m.ER,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_line(aes(x=year-50,y=umsy,group=plotOM),linewidth=1.2,colour='green')+
  geom_ribbon(aes(ymin = l80.er, ymax = u80.er,x=year-50), alpha = 0.2)+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Exploitation Rate')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

er_plot_umsy
ggsave("figs/er_comp_ERtar_0.9UMSY.png",plot=er_plot_umsy,width=8,height=6)

er_plot_bm<-ggplot(meanER_b)+
  geom_line(aes(x=year-50,y=m.ER,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_line(aes(x=year-50,y=umsy,group=plotOM),linewidth=1.2,colour='green')+
  geom_ribbon(aes(ymin = l80.er, ymax = u80.er,x=year-50), alpha = 0.2)+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Exploitation Rate')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

er_plot_bm
ggsave("figs/er_comp_ERtar0.65_0.1BM.png",plot=er_plot_bm,width=8,height=6)

#spawner trends
sp_plot_umsy<-ggplot(spawners_u)+
  geom_line(aes(x=year-50,y=m.spawners,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_ribbon(aes(ymin = l80.s, ymax = u80.s,x=year-50), alpha = 0.2)+
  geom_line(aes(x=year-50,y=smsy,group=plotOM),linewidth=1.2,colour='green')+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Spawner abundance')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot_umsy
ggsave("figs/spawners_timeseries_vs_smsy_ERtar_0.9UMSY.png",plot=sp_plot_umsy,width=8,height=6)

sp_plot_bm<-ggplot(spawners_b)+
  geom_line(aes(x=year-50,y=m.spawners,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_ribbon(aes(ymin = l80.s, ymax = u80.s,x=year-50), alpha = 0.2)+
  geom_line(aes(x=year-50,y=smsy,group=plotOM),linewidth=1.2,colour='green')+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Spawner abundance')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot_bm
ggsave("figs/spawners_timeseries_vs_smsy_ERtar0.65_0.1BM.png",plot=sp_plot_bm,width=8,height=6)

#spawner violins
sp_plot2_umsy<-ggplot(spawners_u_g,aes(x=plotMP,y=m.spawners,fill=plotMP))+
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
  ylab('log10[Spawner Abundance (in Green Status)]')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot2_umsy
ggsave("figs/spawners_violin_ERtar_0.9UMSY.png",plot=sp_plot2_umsy,width=8,height=12)

sp_plot2_bm<-ggplot(spawners_b_g,aes(x=plotMP,y=m.spawners,fill=plotMP))+
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
  ylab('log10[Spawner Abundance (in Green Status)]')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot2_bm
ggsave("figs/spawners_violin_ERtar0.65_0.1BM.png",plot=sp_plot2_bm,width=8,height=12)


#composite plots####
#subset of scenarios:
st_agg<-ggplot(hcrdat1.s[hcrdat1.s$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text.x = element_blank(),strip.text.y=element_blank())
#st_agg
st_agg_all<-ggplot(hcrdat1.s[hcrdat1.s$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('Proportion of years in status')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text.x = element_blank(),strip.text.y=element_blank())

catch_plot_scu<-ggplot(catch_u,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
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


stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  align = "h", axis = "l",ncol=3,
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("figs/composite_plot_status_catch_ERtar_0.9UMSY.png",plot=stclass,width=14,height=8.5)


st_agg<-ggplot(hcrdat2.s[hcrdat2.s$year<111,])+
  geom_bar(aes(x=year-50, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(legend.position = "none")
#st_agg
st_agg_all<-ggplot(hcrdat2.s[hcrdat2.s$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  facet_grid(plotOM~.)+
  ylab('Proportion of years in status')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position = "none")

catch_plot_scu<-ggplot(catch_u2,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
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


stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  align = "h", axis = "l",ncol=3,
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("figs/composite_plot_status_catch_ERtar0.65_0.1BM.png",plot=stclass,width=14,height=8.5)












###


scen=subset(hcrdat,scenario==unique(hcrdat1.5$scenario)[1])
cols=c('navy','dodgerblue','darkred','goldenrod')
par(mfrow=c(1,1))
plot(c(0,1),ylim=c(min(sb$totalCatch),max(sb$totalCatch)),xlim=c(1,36),bty='l',type='n')

for(s in 9:12){
  sb=subset(hcrdat1.5,scenario==unique(hcrdat1.5$scenario)[s])
  for(i in 1:max(sb$iteration)){
    i_sb=subset(sb,iteration==i)
    lines(i_sb$totalCatch~seq(1,36),col=adjustcolor(cols[s-8],alpha.f=0.5))
  }
  
}

er=hcrdat1.5 %>% group_by(scenario,iteration) %>% summarize(m.er=mean(ER),cv.catch=sd(ER)/mean(ER),)


scen=subset(hcrdat1.5,scenario==unique(hcrdat1.5$scenario)[1])
cols=c('navy','dodgerblue','darkred','goldenrod')
par(mfrow=c(2,2))

for(s in 9:12){
  sb=subset(srdat1.5,scenario==unique(hcrdat1.5$scenario)[s])
  plot(c(0,1),ylim=c(0,1),xlim=c(1,36),bty='l',type='n')
for(i in 1:max(sb$iteration)){
    i_sb=subset(sb,iteration==i)
    lines(i_sb$ER~seq(1,36),col=adjustcolor(cols[s-8],alpha.f=0.5))
  }
  
}

#just assessing benchmarks####
simPars1.5 <- read.csv("data/guidelines/SimPars1.5.csv")
cuPar1.5 <- read.csv("data/guidelines/CUPars1.5.csv")

hcrDatalist1.5<-list()
srData1.5<- list()

#1.5 initial alpha scenarios
for(a in seq_len(nrow(simPars1.5))){
  #a=3
  
  hcrDatalist1.5[[a]] <- readRDS(paste0("gdlout/SamSimOutputs/simData/",
                                        simPars1.5$nameOM[a],"/", 
                                        simPars1.5$scenario[a],"/",
                                        paste(simPars1.5$nameOM[a],"_", simPars1.5$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist1.5[[a]]$scenario<-simPars1.5$scenario[a]
  hcrDatalist1.5[[a]]$nameOM<-simPars1.5$nameOM[a]
  hcrDatalist1.5[[a]]$nameMP<-simPars1.5$nameMP[a]
  
  
  srData1.5[[a]] <- readRDS(paste0("gdlout/SamSimOutputs/simData/", 
                                   simPars1.5$nameOM[a],"/",
                                   simPars1.5$scenario[a],"/",
                                   paste(simPars1.5$nameOM[a],"_", simPars1.5$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  
  srData1.5[[a]]$scenario<-simPars1.5$scenario[a]
  srData1.5[[a]]$nameOM<-simPars1.5$nameOM[a]
  srData1.5[[a]]$nameMP<-simPars1.5$nameMP[a]
}


hcrdat1.5 <- do.call(rbind,hcrDatalist1.5)
srdat1.5<- do.call(rbind,srData1.5)
hcrdat1.5<-hcrdat1.5[hcrdat1.5$year>54,]
srdat1.5<-srdat1.5[srdat1.5$year>54,]

hcrdat1.5$aboveLowerassess<-"correct below lower BM"
hcrdat1.5$aboveLowerassess[hcrdat1.5$aboveLowerBM==1& hcrdat1.5$aboveLowerObsBM==1]<-"correct above lower BM"
hcrdat1.5$aboveLowerassess[hcrdat1.5$aboveLowerBM==0& hcrdat1.5$aboveLowerObsBM==1]<-"wrong optimistic"
hcrdat1.5$aboveLowerassess[hcrdat1.5$aboveLowerBM==1& hcrdat1.5$aboveLowerObsBM==0]<-"wrong pessimistic"
hcrdat1.5$aboveUpperassess<-"correct below upper BM"
hcrdat1.5$aboveUpperassess[hcrdat1.5$aboveUpperBM==1& hcrdat1.5$aboveUpperObsBM==1]<-"correct above upper BM"
hcrdat1.5$aboveUpperassess[hcrdat1.5$aboveUpperBM==0& hcrdat1.5$aboveUpperObsBM==1]<-"wrong optimistic"
hcrdat1.5$aboveUpperassess[hcrdat1.5$aboveUpperBM==1& hcrdat1.5$aboveUpperObsBM==0]<-"wrong pessimistic"

head(hcrdat1.5)

abovelower<-ggplot(hcrdat1.5)+
  geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
  scale_color_viridis_d(begin=.1, end=.8) +
  scale_fill_viridis_d(begin=.1, end=.8) +
  facet_grid(nameOM~nameMP)+
  ggtitle("proportion above lower BM")+
  theme_bw(14)
abovelower
ggsave("figs/LB_assess_flatER.png",plot=abovelower)


abovelower<-ggplot(hcrdat2.0)+
  geom_bar(aes(x=year, fill=factor(aboveLowerassess)),position = "fill")+
  scale_color_viridis_d(begin=.1, end=.8) +
  scale_fill_viridis_d(begin=.1, end=.8) +
  facet_grid(nameOM~nameMP)+
  ggtitle("proportion above lower BM")+
  theme_bw(14)
abovelower
ggsave("figs/LB_assess.png",plot=abovelower)

hcrdat1.5$wsp.status<-NA
hcrdat1.5$wsp.status[hcrdat1.5$aboveLowerBM==1&
                    hcrdat1.5$aboveUpperBM==1]<-"green"
hcrdat1.5$wsp.status[hcrdat1.5$aboveLowerBM==1&
                    hcrdat1.5$aboveUpperBM==0]<-"amber"
hcrdat1.5$wsp.status[hcrdat1.5$aboveLowerBM==0&
                    hcrdat1.5$aboveUpperBM==0]<-"red"

hcrdat1.5$status<-NA
hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==1&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==1]<-"green"
hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==0]<-"amber"
hcrdat1.5$status[hcrdat1.5$aboveLowerBM==0&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==0&
                hcrdat1.5$aboveUpperObsBM==0]<-"red"


hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==1&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==0]<-"pessimistic green-> amber"


hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==1]<-"optimistic amber -> green"


hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==0&
                hcrdat1.5$aboveUpperObsBM==0]<-"pessimistic amber -> red"

hcrdat1.5$status[hcrdat1.5$aboveLowerBM==0&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==0]<-"optimistic red -> amber"

hcrdat1.5$status[hcrdat1.5$aboveLowerBM==1&
                hcrdat1.5$aboveUpperBM==1&
                hcrdat1.5$aboveLowerObsBM==0&
                hcrdat1.5$aboveUpperObsBM==0]<-"pessimistic green -> red"

hcrdat1.5$status[hcrdat1.5$aboveLowerBM==0&
                hcrdat1.5$aboveUpperBM==0&
                hcrdat1.5$aboveLowerObsBM==1&
                hcrdat1.5$aboveUpperObsBM==1]<-"optimistic red -> green"


unique(hcrdat1.5$status)

hcrdat1.5$status_agg<-hcrdat1.5$status
hcrdat1.5$status_agg[hcrdat1.5$status %in% c( "pessimistic green-> amber",                       
                                        "pessimistic amber -> red",  
                                        "pessimistic green -> red" )]<- "pessimistic"
hcrdat1.5$status_agg[hcrdat1.5$status %in% c( "optimistic amber -> green", 
                                        "optimistic red -> amber",  
                                        "optimistic red -> green")]<- "optimistic"

hcrdat1.5$wsp.status<- factor(hcrdat1.5$wsp.status,levels=c('red','amber','green'))
hcrdat1.5$status_agg<-factor(hcrdat1.5$status_agg,levels=c("red",
                                                     "amber",
                                                     "green",
                                                     "pessimistic",
                                                     "optimistic"))



#these palletes are color bling friendlish. -- second one is a bit harder to see in black and white. Change with caution.
statusCols <- c("#9A3F3F","#DFD98D","#8EB687","gray75","gray25")
statusColsall <- c("#9A3F3F","gray90","gray10","#DFD98D","gray80","gray20","#8EB687","gray70","gray30")
hcrdat1.5$status<-factor(hcrdat1.5$status,levels=c("red",
                                             "optimistic red -> amber",
                                             "pessimistic amber -> red",                                                      
                                             "amber",
                                             "optimistic amber -> green",
                                             "pessimistic green-> amber",
                                             "green",
                                             "optimistic red -> green",
                                             "pessimistic green -> red"
))

#reduced by one scenario - include prod 2-0.5 and 2-1


#hcrdat1.52=subset(hcrdat1.5,nameOM!='decLinearProd1.5to0.5')

#unique(hcrdat1.52$nameOM)
hcrdat1.5$plotOM<-dplyr::recode_factor(factor(hcrdat1.5$nameOM),"stationaryhAR1" = "AR1 high",
                                    "stationarylAR1" = "AR1 low",  
                                    "decLinearcap0.5"="cap * 0.5", 
                                    "decLinearcap0.25"="cap * 0.25",
                                    "decLinearProd2to1.5"="prod 2 -> 1.5",     
                                    "decLinearProd2to1"="prod 2 -> 1",
                                    "decLinearProd2to0.6"="prod 2 -> 0.5")

hcrdat1.5$plotMP<-dplyr::recode_factor(factor(hcrdat1.5$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                    "10yr_rwa_constER" = "Time-varying (10y)",  
                                    "5yr_autocorr_constER" = "Stationary (5y)",
                                    "5yr_rwa_constER" = "Time-varying (5y)",
                                    "10yr_autocorr_adaptER" = "Stationary (10y)",
                                    "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                    "5yr_autocorr_adaptER" = "Stationary (5y)",
                                    "5yr_rwa_adaptER" = "Time-varying (5y)")


st_agg<-ggplot(hcrdat1.5)+
  geom_bar(aes(x=year-50, fill=status_agg),position = "fill")+
  scale_fill_manual(values = 
                      statusCols)+
  facet_grid(plotOM~plotMP)+
  ylab('Proportion of simulations')+
  xlab('Years of simulation')+
  #  ggtitle("status")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
#st_agg
st_agg_all<-ggplot(hcrdat1.5[hcrdat1.5$year>69,])+
  geom_bar(aes(x=plotMP, fill=status_agg),position = "fill")+
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
ggsave("figs/status_classification_nofeedback.png",plot=stclass,width=14,height=8.5)


##true status - no feedback scenarios
st_agg<-ggplot(hcrdat1.5)+
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
st_agg_all<-ggplot(hcrdat1.5[hcrdat1.5$year>69,])+
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
ggsave("figs/true_status_nofeedback.png",plot=stclass,width=14,height=8.5)


