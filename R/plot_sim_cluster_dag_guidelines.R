
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



#gudelines scenarios example####
simPars_um <- read.csv("data/guidelines/SimPars2.0.csv") #simpars for ER tracking umsy, no EG
cuPar <- read.csv("data/guidelines/CUPars2.0.csv") #cu pars
simPars_b <- read.csv("data/guidelines/SimPars2.1.csv") #simpars for fixed ER, assessed EG
simPars_c <- read.csv("data/guidelines/SimPars2.2.csv") #simpars for ER tracking umsy, assessed EG
simPars_d <- read.csv("data/guidelines/SimPars2.3.csv") #simpars for ER tracking umsy, fixed EG
#here()

hcrDatalist_um<-list()
srData_um<- list()

for(a in seq_len(nrow(simPars_um))){

   hcrDatalist_um[[a]] <- readRDS(paste0("./outputs/umsy_track/SamSimOutputs/simData/",
   simPars_um$nameOM[a],"/", 
   simPars_um$scenario[a],"/",
   paste(simPars_um$nameOM[a],"_", simPars_um$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout

   hcrDatalist_um[[a]]$scenario<-simPars_um$scenario[a]
   hcrDatalist_um[[a]]$nameOM<-simPars_um$nameOM[a]
   hcrDatalist_um[[a]]$nameMP<-simPars_um$nameMP[a]


   srData_um[[a]] <- readRDS(paste0("./outputs/umsy_track/SamSimOutputs/simData/", 
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
  
  hcrDatalist_b[[a]] <- readRDS(paste0("./outputs/bm_track/SamSimOutputs/simData/",
                                        simPars_b$nameOM[a],"/", 
                                        simPars_b$scenario[a],"/",
                                        paste(simPars_b$nameOM[a],"_", simPars_b$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist_b[[a]]$scenario<-simPars_b$scenario[a]
  hcrDatalist_b[[a]]$nameOM<-simPars_b$nameOM[a]
  hcrDatalist_b[[a]]$nameMP<-simPars_b$nameMP[a]
  
  
  srData_b[[a]] <- readRDS(paste0("./outputs/bm_track/SamSimOutputs/simData/", 
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

simPars_c <- read.csv("data/guidelines/SimPars2.2.csv")

hcrDatalist_c<-list()
srData_c<- list()

for(a in seq_len(nrow(simPars_c))){
  
  hcrDatalist_c[[a]] <- readRDS(paste0("./outputs/umsy_bm_track/SamSimOutputs/simData/",
                                       simPars_c$nameOM[a],"/", 
                                       simPars_c$scenario[a],"/",
                                       paste(simPars_c$nameOM[a],"_", simPars_c$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist_c[[a]]$scenario<-simPars_c$scenario[a]
  hcrDatalist_c[[a]]$nameOM<-simPars_c$nameOM[a]
  hcrDatalist_c[[a]]$nameMP<-simPars_c$nameMP[a]
  
  
  srData_c[[a]] <- readRDS(paste0("./outputs/umsy_bm_track/SamSimOutputs/simData/", 
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

hcrdat3=format_dat(hcrdat_c)


hcrDatalist_d<-list()
srData_d<- list()

for(a in seq_len(nrow(simPars_d))){
  
  hcrDatalist_d[[a]] <- readRDS(paste0("./outputs/umsy_fixedbm/SamSimOutputs/simData/",
                                       simPars_d$nameOM[a],"/", 
                                       simPars_d$scenario[a],"/",
                                       paste(simPars_d$nameOM[a],"_", simPars_d$nameMP[a], "_", "CU_HCR_PM.RData",sep="")))$hcrDatout
  
  hcrDatalist_d[[a]]$scenario<-simPars_d$scenario[a]
  hcrDatalist_d[[a]]$nameOM<-simPars_d$nameOM[a]
  hcrDatalist_d[[a]]$nameMP<-simPars_d$nameMP[a]
  
  
  srData_d[[a]] <- readRDS(paste0("./outputs/umsy_fixedbm/SamSimOutputs/simData/", 
                                  simPars_d$nameOM[a],"/",
                                  simPars_d$scenario[a],"/",
                                  paste(simPars_d$nameOM[a],"_", simPars_d$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  
  srData_d[[a]]$scenario<-simPars_d$scenario[a]
  srData_d[[a]]$nameOM<-simPars_d$nameOM[a]
  srData_d[[a]]$nameMP<-simPars_d$nameMP[a]
}

hcrdat_d <- do.call(rbind,hcrDatalist_d)
srdat_d<- do.call(rbind,srData_d)

hcrdat_d<-hcrdat_d[hcrdat_d$year>50,]
srdat_d<-srdat_d[srdat_d$year>50,]

hcrdat4=format_dat(hcrdat_d)

#recode scenarios
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

hcrdat4$plotOM<-dplyr::recode_factor(factor(hcrdat4$plotOM),
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


#true status - ER feedback scenarios####

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
ggsave("./outputs/figs/true_status_ERtar_1.0UMSY.png",plot=stclass,width=14,height=8.5)

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
st_agg_all<-ggplot(hcrdat2[hcrdat2$year>59&hcrdat2$year<111,])+
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
ggsave("./outputs/figs/true_status_ERtar0.75_0.5BM.png",plot=stclass,width=14,height=8.5)


#abundance&umsy benchmarks
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
ggsave("outputs/figs/true_status_ERtar0.9UMSY_0.1BM.png",plot=stclass,width=14,height=8.5)

#umsy & fixed benchmark
st_agg<-ggplot(hcrdat4[hcrdat4$year<111,])+
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
st_agg_all<-ggplot(hcrdat4[hcrdat4$year<111,])+
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
ggsave("outputs/figs/true_status_ERtar0.9UMSY_0.1fixedBM.png",plot=stclass,width=14,height=8.5)

#Guidance doc outputs with a limited set of scenarios####

#true status subset - abundance benchmark
hcrdat2x=subset(hcrdat2,plotOM %in% c('Static - 0.4AR1','Static - 0.8AR1','-50%Cap','-75%Cap','-50%Prod','-75%Prod'))
hcrdat2x$plotOM=droplevels(hcrdat2x$plotOM)

st_agg<-ggplot(hcrdat2x[hcrdat2x$year<111,])+
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
st_agg_all<-ggplot(hcrdat2x[hcrdat2x$year>59&hcrdat2x$year<111,])+
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
ggsave("./outputs/figs/GDsub_true_status_ERtar0.75_0.3BM.png",plot=stclass,width=14,height=8.5)

#true status subset - F benchmark
hcrdat4x=subset(hcrdat4,plotOM %in% c('Static - 0.4AR1','Static - 0.8AR1','-50%Cap','-75%Cap','-50%Prod','-75%Prod'))
hcrdat4x$plotOM=droplevels(hcrdat4x$plotOM)

st_agg<-ggplot(hcrdat4x[hcrdat4x$year<111,])+
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
st_agg_all<-ggplot(hcrdat4x[hcrdat4x$year>59&hcrdat4x$year<111,])+
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
ggsave("./outputs/figs/GDsub_true_status_ERtarUmsy_fixEG.png",plot=stclass,width=14,height=8.5)


##estimated status - ER feedback scenarios####


st_agg<-ggplot(hcrdat2x[hcrdat2x$year<111,])+
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
st_agg_all<-ggplot(hcrdat2x[hcrdat2x$year>59&hcrdat2x$year<111,])+
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
ggsave("outputs/figs/status_classification_ERtar0.75_0.3BM.png",plot=stclass,width=14,height=8.5)






##true vs obs benchmarks and umsy####
BMcomp2<-data.frame(year=hcrdat2x$year,
  iteration=hcrdat2x$iteration,
  typeest=as.factor(rep(c("True","Est."),each=length(hcrdat2x$year))),
  lowerBM=c(hcrdat2x$lowerBM,hcrdat2x$lowerObsBM),
  upperBM=c(hcrdat2x$upperBM,hcrdat2x$upperObsBM),
  nameOM=hcrdat2x$nameOM,
  nameMP=hcrdat2x$plotMP,
  plotOM=hcrdat2x$plotOM
  )


upperBMtrue_hcr<-ggplot(BMcomp2)+
geom_line(aes(x=year-50,y=upperBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d('',begin=.1, end=.8) +
facet_grid(plotOM~nameMP)+
ylab('Upper Benchmark')+
theme_bw(14)
upperBMtrue_hcr
ggsave("outputs/figs/upperBM_true_hcr_ERtar0.75_0.5BM.png",plot=upperBMtrue_hcr,width=14,height=8.5)

lowerBMtrue_hcr<-ggplot(BMcomp2)+
geom_line(aes(x=year-50,y=lowerBM,color=typeest,group=interaction(typeest, iteration)),linewidth=1.2, alpha=.4)+
scale_color_viridis_d('',begin=.1, end=.8) +
facet_grid(plotOM~nameMP)+
coord_cartesian(ylim=c(0,70000))+
ylab('Lower Benchmark')+
theme_bw(14)
lowerBMtrue_hcr
ggsave("outputs/figs/lowerBM_true_hcr_ERtar0.75_0.5BM.png",plot=lowerBMtrue_hcr,width=14,height=8.5)

BMcompt=subset(BMcomp2,typeest=='truth')

bm_plot1<-ggplot(BMcompt)+
  geom_line(aes(x=year-50,y=upperBM,color=plotOM,group=plotOM),linewidth=1.2)+
  scale_color_manual(values=c('#a1dab4','#d95f0e','#fed98e','#253494','#41b6c4'),name='') +
  annotate('text',x=71,y=unique(BMcompt$upperBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[1]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[1],size=4,colour='#a1dab4',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$upperBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[2]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[2],size=4,colour='#d95f0e',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$upperBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[3]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[3],size=4,colour='#fed98e',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$upperBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[4]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[4],size=4,colour='#253494',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$upperBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[5]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[5],size=4,colour='#41b6c4',hjust = 0)+
  coord_cartesian(clip = 'off')+
  ylab('Upper benchmark: 80% Smsy')+
  xlab(
    ''
  )+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),plot.margin = margin(0.5,3,0.5,0.5, "cm"))

bm_plot2<-ggplot(BMcompt)+
  geom_line(aes(x=year-50,y=lowerBM,color=plotOM,group=plotOM),linewidth=1.2)+
  scale_color_manual(values=c('#a1dab4','#d95f0e','#fed98e','#253494','#41b6c4')) +
  annotate('text',x=71,y=unique(BMcompt$lowerBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[1]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[1],size=4,colour='#a1dab4',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$lowerBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[2]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[2],size=4,colour='#d95f0e',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$lowerBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[3]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[3],size=4,colour='#fed98e',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$lowerBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[4]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[4],size=4,colour='#253494',hjust = 0)+
  annotate('text',x=71,y=unique(BMcompt$lowerBM[BMcompt$plotOM==levels(factor(BMcompt$plotOM))[5]&BMcompt$year==120]),label=levels(factor(BMcompt$plotOM))[5],size=4,colour='#41b6c4',hjust = 0)+
  coord_cartesian(xlim = c(0, 70),clip = 'off')+
   ylab('Lower benchmark: Sgen')+
  xlab(
    'Year of simulation'
  )+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),plot.margin = margin(0.5,3,0.5,0.5, "cm"))

bm_trends1<-plot_grid(
  bm_plot1+ theme(legend.position="none"),
bm_plot2+ theme(legend.position="none"),
ncol=1,nrow=2
)
bm_trends1
ggsave("outputs/figs/trueBM_comp.png",plot=bm_trends1,width=6,height=8)


#Fishery Risk - total catch, CV catch,####
srdat_um$plotOM<-dplyr::recode_factor(factor(srdat_um$nameOM),
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


srdat_um$plotMP<-dplyr::recode_factor(factor(srdat_um$nameMP),
                                    "10yr_autocorr_constER" = "Stationary (10y)",
                                    "10yr_rwa_constER" = "Time-varying (10y)",  
                                    "5yr_autocorr_constER" = "Stationary (5y)",
                                    "5yr_rwa_constER" = "Time-varying (5y)",
                                    "10yr_autocorr_adaptER" = "Stationary (10y)",
                                    "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                    "5yr_autocorr_adaptER" = "Stationary (5y)",
                                    "5yr_rwa_adaptER" = "Time-varying (5y)")

srdat_b$plotOM<-dplyr::recode_factor(factor(srdat_um$nameOM),
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

srdat_b$plotMP<-dplyr::recode_factor(factor(srdat_b$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                      "10yr_rwa_constER" = "Time-varying (10y)",  
                                      "5yr_autocorr_constER" = "Stationary (5y)",
                                      "5yr_rwa_constER" = "Time-varying (5y)",
                                      "10yr_autocorr_adaptER" = "Stationary (10y)",
                                      "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                      "5yr_autocorr_adaptER" = "Stationary (5y)",
                                      "5yr_rwa_adaptER" = "Time-varying (5y)")

srdat_c$plotOM<-dplyr::recode_factor(factor(srdat_um$nameOM),
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


srdat_c$plotMP<-dplyr::recode_factor(factor(srdat_c$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                     "10yr_rwa_constER" = "Time-varying (10y)",  
                                     "5yr_autocorr_constER" = "Stationary (5y)",
                                     "5yr_rwa_constER" = "Time-varying (5y)",
                                     "10yr_autocorr_adaptER" = "Stationary (10y)",
                                     "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                     "5yr_autocorr_adaptER" = "Stationary (5y)",
                                     "5yr_rwa_adaptER" = "Time-varying (5y)")

srdat_d$plotOM<-dplyr::recode_factor(factor(srdat_d$nameOM),
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


srdat_d$plotMP<-dplyr::recode_factor(factor(srdat_d$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                     "10yr_rwa_constER" = "Time-varying (10y)",  
                                     "5yr_autocorr_constER" = "Stationary (5y)",
                                     "5yr_rwa_constER" = "Time-varying (5y)",
                                     "10yr_autocorr_adaptER" = "Stationary (10y)",
                                     "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                     "5yr_autocorr_adaptER" = "Stationary (5y)",
                                     "5yr_rwa_adaptER" = "Time-varying (5y)")


#total catch and variance in annual catch
catch1_u=hcrdat1[hcrdat1$year>59&hcrdat1$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_u=hcrdat1[hcrdat1$year>59&hcrdat1$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_u$scale.ann.catch=catch_u$m.catch/catch1_u$max.catch[match(catch_u$plotOM,catch1_u$plotOM)]

catch1_b=hcrdat2[hcrdat2$year>59&hcrdat2$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_b=hcrdat2[hcrdat2$year>59&hcrdat2$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_b$scale.ann.catch=catch_b$m.catch/catch1_b$max.catch[match(catch_b$plotOM,catch1_b$plotOM)]

catch1_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_c$scale.ann.catch=catch_c$m.catch/catch1_c$max.catch[match(catch_c$plotOM,catch1_c$plotOM)]

catch1_d=hcrdat4[hcrdat4$year>59&hcrdat4$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(m.catch=exp(mean(log(totalCatch)))) %>% group_by(plotOM) %>% summarize(max.catch=max(m.catch))
catch_d=hcrdat4[hcrdat4$year>59&hcrdat4$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% summarize(total.catch=sum(log(totalCatch)),m.catch=exp(mean(log(totalCatch))),cv.catch=sd(totalCatch)/exp(mean(log(totalCatch))))
catch_d$scale.ann.catch=catch_d$m.catch/catch1_d$max.catch[match(catch_d$plotOM,catch1_d$plotOM)]

#proportion above benchmarks
bmdat_u= hcrdat1[hcrdat1$year>59&hcrdat1$year<111,] %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_b= hcrdat2[hcrdat2$year>59&hcrdat2$year<111,] %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_c= hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_d= hcrdat4[hcrdat4$year>59&hcrdat4$year<111,] %>% group_by(plotOM,plotMP,iteration)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())

bmdat_tot_u= hcrdat1[hcrdat1$year>59&hcrdat1$year<111,] %>% group_by(plotOM,plotMP)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_tot_b=hcrdat2[hcrdat2$year>59&hcrdat2$year>59&hcrdat2$year<111,] %>% group_by(plotOM,plotMP)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_tot_c=hcrdat3[hcrdat3$year>59&hcrdat3$year<111,] %>% group_by(plotOM,plotMP)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())
bmdat_tot_d= hcrdat4[hcrdat4$year>59&hcrdat4$year<111,] %>% group_by(plotOM,plotMP)  %>% summarize(p.red=1-sum(aboveLowerBM)/n(),p.bm=sum(aboveUpperBM)/n())

#ER series
meanER_u=srdat_um[srdat_um$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_u$scnmp.t=paste(meanER_u$plotOM,meanER_u$plotMP,meanER_u$year,sep="_")
meanER_u=distinct(meanER_u,scnmp.t,.keep_all=T)
meanER_b=srdat_b[srdat_b$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_b$scnmp.t=paste(meanER_b$plotOM,meanER_b$plotMP,meanER_b$year,sep="_")
meanER_b=distinct(meanER_b,scnmp.t,.keep_all=T)

meanER_c=srdat_c[srdat_c$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_c$scnmp.t=paste(meanER_c$plotOM,meanER_c$plotMP,meanER_c$year,sep="_")
meanER_c=distinct(meanER_c,scnmp.t,.keep_all=T)

meanER_d=srdat_d[srdat_d$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.ER=mean(ER),sd.er=sd(ER),l80.er=quantile(ER,0.1),u80.er=quantile(ER,0.9),umsy=uMSY)
meanER_d$scnmp.t=paste(meanER_d$plotOM,meanER_d$plotMP,meanER_d$year,sep="_")
meanER_d=distinct(meanER_d,scnmp.t,.keep_all=T)

#spawners through time
spawners_u1=srdat_um[srdat_um$year>59&srdat_um$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_u1$mpy=paste(spawners_u1$plotOM,spawners_u1$plotMP,spawners_u1$year)
spawners_u1=distinct(spawners_u1,mpy,.keep_all=TRUE)
spawners_b1=srdat_b[srdat_b$year>59&srdat_b$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_b1$mpy=paste(spawners_b1$plotOM,spawners_b1$plotMP,spawners_b1$year)
spawners_b1=distinct(spawners_b1,mpy,.keep_all=TRUE)
spawners_c1=srdat_c[srdat_c$year>59&srdat_c$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_c1$mpy=paste(spawners_c1$plotOM,spawners_c1$plotMP,spawners_c1$year)
spawners_c1=distinct(spawners_c1,mpy,.keep_all=TRUE)
spawners_d1=srdat_d[srdat_d$year>59&srdat_d$year<111,] %>% group_by(plotOM,plotMP,year) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen)
spawners_d1$mpy=paste(spawners_d1$plotOM,spawners_d1$plotMP,spawners_d1$year)
spawners_d1=distinct(spawners_d1,mpy,.keep_all=TRUE)

#total spawners
spawners_u2=srdat_um[srdat_um$year>59&srdat_um$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_u2$mpy=paste(spawners_u2$plotOM,spawners_u2$plotMP,spawners_u2$iteration)
spawners_u2=distinct(spawners_u2,mpy,.keep_all=TRUE)
spawners_b2=srdat_b[srdat_b$year>59&srdat_b$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_b2$mpy=paste(spawners_b2$plotOM,spawners_b2$plotMP,spawners_b2$iteration)
spawners_b2=distinct(spawners_b2,mpy,.keep_all=TRUE)
spawners_c2=srdat_c[srdat_c$year>59&srdat_c$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_c2$mpy=paste(spawners_c2$plotOM,spawners_c2$plotMP,spawners_c2$iteration)
spawners_c2=distinct(spawners_c2,mpy,.keep_all=TRUE)
spawners_d2=srdat_d[srdat_d$year>59&srdat_d$year<111,] %>% group_by(plotOM,plotMP,iteration) %>% reframe(m.spawners=exp(mean(log(spawners))),l80.s=quantile(spawners,0.1),u80.s=quantile(spawners,0.9),smsy=sMSY,sgen=sGen,cv.spawners=sd(spawners)/exp(mean(log(spawners))))
spawners_d2$mpy=paste(spawners_d2$plotOM,spawners_d2$plotMP,spawners_d2$iteration)
spawners_d2=distinct(spawners_d2,mpy,.keep_all=TRUE)


  

##comparison of catch:####
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

catch_plot_scu
ggsave("outputs/figs/catch_violins_ERtar_0.9Umsy.png",plot=catch_plot_scu,width=8,height=12)

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
ggsave("outputs/figs/catch_violins_ERtar0.65_0.1BM.png",plot=catch_plot_scb,width=8,height=12)

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

catch_plot_scd<-ggplot(catch_d,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
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

catch_plot_scd
ggsave("outputs/figs/catch_violins_ERtar0.9UMSY_0.1fixedBM.png",plot=catch_plot_scd,width=8,height=12)


#Summary plot: proportional change in catch vs. conservation risk (% below upper BM)#####
#scenarios
hcrdat1.s=subset(hcrdat1,plotOM %in% c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))
hcrdat1.s$plotOM=droplevels(hcrdat1.s$plotOM)
hcrdat1.s=subset(hcrdat1.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)"))
hcrdat1.s$plotMP=droplevels(hcrdat1.s$plotMP)
hcrdat1.s$plotMP=dplyr::recode_factor(factor(hcrdat1.s$plotMP),"Stationary (5y)" = "Stationary",
                     "Time-varying (5y)" = "Time-varying")

#total catch and variance in annual catch
hcrdat2.s=subset(hcrdat2,plotOM %in% c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))
hcrdat2.s$plotOM=droplevels(hcrdat2.s$plotOM)
hcrdat2.s=subset(hcrdat2.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)"))
hcrdat2.s$plotMP=droplevels(hcrdat2.s$plotMP)
hcrdat2.s$plotMP=dplyr::recode_factor(factor(hcrdat2.s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

#total catch and variance in annual catch
hcrdat3.s=subset(hcrdat3,plotOM %in%  c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))
hcrdat3.s$plotOM=droplevels(hcrdat3.s$plotOM)
hcrdat3.s=subset(hcrdat3.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)"))
hcrdat3.s$plotMP=droplevels(hcrdat3.s$plotMP)
hcrdat3.s$plotMP=dplyr::recode_factor(factor(hcrdat3.s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

#total catch and variance in annual catch
hcrdat4.s=subset(hcrdat4,plotOM %in%  c('Static - 0.4AR1','-75%Cap','-50%Prod','-75%Prod'))
hcrdat4.s$plotOM=droplevels(hcrdat4.s$plotOM)
hcrdat4.s=subset(hcrdat4.s,plotMP %in% c("Stationary (5y)" ,"Time-varying (5y)"))
hcrdat4.s$plotMP=droplevels(hcrdat4.s$plotMP)
hcrdat4.s$plotMP=dplyr::recode_factor(factor(hcrdat4.s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

#subset of catch data
catch_u_s=subset(catch_u,plotOM %in% levels(hcrdat1.s$plotOM))
catch_u_s=subset(catch_u_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
catch_u_s$plotOM=droplevels(catch_u_s$plotOM)
catch_u_s$plotOM=droplevels(catch_u_s$plotOM)
catch_u_s$plotMP=dplyr::recode_factor(factor(catch_u_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")
catch_b_s=subset(catch_b,plotOM %in% levels(hcrdat2.s$plotOM))
catch_b_s=subset(catch_b_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
catch_b_s$plotOM=droplevels(catch_b_s$plotOM)
catch_b_s$plotOM=droplevels(catch_b_s$plotOM)
catch_b_s$plotMP=dplyr::recode_factor(factor(catch_b_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

catch_c_s=subset(catch_c,plotOM %in% levels(hcrdat3.s$plotOM))
catch_c_s=subset(catch_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
catch_c_s$plotOM=droplevels(catch_c_s$plotOM)
catch_c_s$plotOM=droplevels(catch_c_s$plotOM)
catch_c_s$plotMP=dplyr::recode_factor(factor(catch_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

catch_d_s=subset(catch_c,plotOM %in% levels(hcrdat4.s$plotOM))
catch_d_s=subset(catch_d_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
catch_d_s$plotOM=droplevels(catch_d_s$plotOM)
catch_d_s$plotOM=droplevels(catch_d_s$plotOM)
catch_d_s$plotMP=dplyr::recode_factor(factor(catch_d_s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

#subset of benchmark data
bmdat_u_s=subset(bmdat_u,plotOM %in% levels(hcrdat1.s$plotOM))
bmdat_u_s=subset(bmdat_u_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
bmdat_u_s$plotOM=droplevels(bmdat_u_s$plotOM)
bmdat_u_s$plotOM=droplevels(bmdat_u_s$plotOM)
bmdat_u_s$plotMP=dplyr::recode_factor(factor(bmdat_u_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")
bmdat_b_s=subset(bmdat_b,plotOM %in% levels(hcrdat2.s$plotOM))
bmdat_b_s=subset(bmdat_b_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
bmdat_b_s$plotOM=droplevels(bmdat_b_s$plotOM)
bmdat_b_s$plotOM=droplevels(bmdat_b_s$plotOM)
bmdat_b_s$plotMP=dplyr::recode_factor(factor(bmdat_b_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

bmdat_c_s=subset(bmdat_c,plotOM %in% levels(hcrdat3.s$plotOM))
bmdat_c_s=subset(bmdat_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
bmdat_c_s$plotOM=droplevels(bmdat_c_s$plotOM)
bmdat_c_s$plotOM=droplevels(bmdat_c_s$plotOM)
bmdat_c_s$plotMP=dplyr::recode_factor(factor(bmdat_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

bmdat_d_s=subset(bmdat_c,plotOM %in% levels(hcrdat4.s$plotOM))
bmdat_d_s=subset(bmdat_d_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
bmdat_d_s$plotOM=droplevels(bmdat_d_s$plotOM)
bmdat_d_s$plotOM=droplevels(bmdat_d_s$plotOM)
bmdat_d_s$plotMP=dplyr::recode_factor(factor(bmdat_d_s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

#subset of spawner data
spawners_u_s=subset(spawners_u2,plotOM %in% levels(hcrdat1.s$plotOM))
spawners_u_s=subset(spawners_u_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
spawners_u_s$plotOM=droplevels(spawners_u_s$plotOM)
spawners_u_s$plotOM=droplevels(spawners_u_s$plotOM)
spawners_u_s$plotMP=dplyr::recode_factor(factor(spawners_u_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")
spawners_b_s=subset(spawners_b2,plotOM %in% levels(hcrdat2.s$plotOM))
spawners_b_s=subset(spawners_b_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
spawners_b_s$plotOM=droplevels(spawners_b_s$plotOM)
spawners_b_s$plotOM=droplevels(spawners_b_s$plotOM)
spawners_b_s$plotMP=dplyr::recode_factor(factor(spawners_b_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

spawners_c_s=subset(spawners_c2,plotOM %in% levels(hcrdat3.s$plotOM))
spawners_c_s=subset(spawners_c_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
spawners_c_s$plotOM=droplevels(spawners_c_s$plotOM)
spawners_c_s$plotOM=droplevels(spawners_c_s$plotOM)
spawners_c_s$plotMP=dplyr::recode_factor(factor(spawners_c_s$plotMP),"Stationary (5y)" = "Stationary",
                                      "Time-varying (5y)" = "Time-varying")

spawners_d_s=subset(spawners_d2,plotOM %in% levels(hcrdat4.s$plotOM))
spawners_d_s=subset(spawners_d_s, plotMP %in% c("Stationary (5y)",'Time-varying (5y)'))
spawners_d_s$plotOM=droplevels(spawners_d_s$plotOM)
spawners_d_s$plotOM=droplevels(spawners_d_s$plotOM)
spawners_d_s$plotMP=dplyr::recode_factor(factor(spawners_d_s$plotMP),"Stationary (5y)" = "Stationary",
                                         "Time-varying (5y)" = "Time-varying")

##To do: finish summary plot by just comparing 10-y static to 10-y rwa

scn_cols_static=c('#045a8d',
  '#41b6c4',
  '#d95f0e',
  '#b30000','#58508d','#bc5090')

scn_cols_rwa=c('#fed98e',
  '#fe9929',
  '#d95f0e',
  '#993404')



##0.9UMSY track####

catch_comp_d=expand.grid(mp=levels(hcrdat1.s$plotMP)[-1],scn=levels(hcrdat1.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat1.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat1.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_u_s,plotOM==levels(catch_u_s$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat1.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1

}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat1.s$plotMP)[-1],scn=levels(hcrdat1.s$plotOM),d.ubm=NA,l95=NA,u95=NA,d.red=NA,l95r=NA,u95r=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat1.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat1.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_u_s,plotOM==levels(bmdat_u_s$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  ests2=lm(1-p.red~plotMP-1,data=subset(bmdat_u,plotOM==levels(bmdat_u$plotOM)[i]))
  confs2=confint(ests2)
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),6]=ests2$coefficients[2]/ests2$coefficients[1]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),7]=confs2[2,1]/confs2[1,2]-1
  status_comp_d[match(levels(hcrdat1.s$plotOM)[i],status_comp_d$scn),8]=confs2[2,2]/confs2[1,1]-1
}

spawners_comp_d=expand.grid(mp=levels(hcrdat1.s$plotMP)[-1],scn=levels(hcrdat1.s$plotOM),d.spwnrs=NA,l95=NA,u95=NA)
spawners_comp_d=as.data.frame(spawners_comp_d)
#  data.frame(scn=levels(hcrdat1.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat1.s$plotOM))){
  ests=lm(m.spawners~plotMP-1,data=subset(spawners_u_s,plotOM==levels(spawners_u_s$plotOM)[i]))
  confs=confint(ests)
  spawners_comp_d[match(levels(hcrdat1.s$plotOM)[i],spawners_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  spawners_comp_d[match(levels(hcrdat1.s$plotOM)[i],spawners_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  spawners_comp_d[match(levels(hcrdat1.s$plotOM)[i],spawners_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  }

png(filename='outputs/figs/catch_conservation_tradeoff_UBM_ERtar_0.9Umsy.png',width=8,height=6,units='in',res=600)
plot(c(-20,120)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-20,120))
abline(h=0)
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=60,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=80,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=100,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=120,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(status_comp_d$d.ubm[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(status_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.94,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.ubm[i]*100~1.25+j,data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above UBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_LBM_ERtar_0.9Umsy.png',width=8,height=6,units='in',res=600)
plot(c(-10,70)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-10,70))
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
  lines(c(status_comp_d$d.red[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(status_comp_d$l95r[i]*100,status_comp_d$u95r[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.red[i]*100~1.25+j,data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above LBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_spawners_ERtar_0.9Umsy.png',width=8,height=6,units='in',res=600)
plot(c(-10,210)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-10,210))
abline(h=0)
abline(h=50,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=100,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=150,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=200,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(spawners_comp_d$d.spwnrs[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95[i]*100,spawners_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~1.75+j,2,data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.spwnrs[i]*100~1.25+j,2,data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Mean annual spawners',line=1.5)
dev.off()


##assessed Abundance benchmark ####
catch_comp_d=expand.grid(mp=levels(hcrdat2.s$plotMP)[-1],scn=levels(hcrdat2.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat2.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat2.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_b_s,plotOM==levels(catch_b_s$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat2.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  
}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat2.s$plotMP)[-1],scn=levels(hcrdat2.s$plotOM),d.ubm=NA,l95=NA,u95=NA,d.red=NA,l95r=NA,u95r=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat2.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat2.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_b_s,plotOM==levels(bmdat_b_s$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  ests2=lm(1-p.red~plotMP-1,data=subset(bmdat_b_s,plotOM==levels(bmdat_b_s$plotOM)[i]))
  confs2=confint(ests2)
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),6]=ests2$coefficients[2]/ests2$coefficients[1]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),7]=confs2[2,1]/confs2[1,2]-1
  status_comp_d[match(levels(hcrdat2.s$plotOM)[i],status_comp_d$scn),8]=confs2[2,2]/confs2[1,1]-1
}

spawners_comp_d=expand.grid(mp=levels(hcrdat2.s$plotMP)[-1],scn=levels(hcrdat2.s$plotOM),d.spwnrs=NA,l95=NA,u95=NA)
spawners_comp_d=as.data.frame(spawners_comp_d)
#  data.frame(scn=levels(hcrdat2.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat2.s$plotOM))){
  ests=lm(m.spawners~plotMP-1,data=subset(spawners_b_s,plotOM==levels(spawners_b_s$plotOM)[i]))
  confs=confint(ests)
  spawners_comp_d[match(levels(hcrdat2.s$plotOM)[i],spawners_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  spawners_comp_d[match(levels(hcrdat2.s$plotOM)[i],spawners_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  spawners_comp_d[match(levels(hcrdat2.s$plotOM)[i],spawners_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
}


png(filename='outputs/figs/catch_conservation_tradeoff_UBM_ERtar0.65_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-50,40)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-50,30))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-50,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(status_comp_d$d.ubm[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(status_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.12*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~c(1.75+j),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.ubm[i]*100~c(1.25+j),data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above UBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_LBM_ERtar0.65_0.1BM.png',width=8,height=6,units='in',res=600)
plot(c(-30,30)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-30,30))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(status_comp_d$d.red[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(status_comp_d$l95r[i]*100,status_comp_d$u95r[i]*100)~rep(1.25+j,1))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~c(1.75+j),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.red[i]*100~c(1.25+j),data=status_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Proportion of time above LBM',line=1.5)
dev.off()

png(filename='outputs/figs/catch_conservation_tradeoff_spawner_ERtar0.75_0.3BM.png',width=8,height=6,units='in',res=600)
plot(c(-40,20)~c(1,2),bty='l',type='n',xaxt='n',ylab='% Difference from stationary',xlab='',ylim=c(-40,20))
abline(h=0)
abline(h=10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-10,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-20,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-30,lty=5,col=adjustcolor('black',alpha.f = 0.5))
abline(h=-40,lty=5,col=adjustcolor('black',alpha.f = 0.5))
for(i in 1:nrow(status_comp_d)){
  j=jitter(0,amount=0.025)
  lines(c(spawners_comp_d$d.spwnrs[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.16*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~c(1.75+j),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.spwnrs[i]*100~c(1.25+j),data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  }
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Mean annual spawners',line=1.5)
dev.off()


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


## uMSY & fixed BM####
catch_comp_d=expand.grid(mp=levels(hcrdat4.s$plotMP)[-1],scn=levels(hcrdat4.s$plotOM),d.catch=NA,l95=NA,u95=NA)
#  data.frame(scn=levels(hcrdat4.s$plotOM),catch.10y.s=NA,catch.10y.rwa=NA,catch.10y.rwa.l95=NA,catch.10y.rwa.u95=NA,catch.5y.s=NA,catch.5y.s.l95=NA,catch.5y.s.u95=NA,catch.5y.rwa=NA,catch.5y.rwa.l95=NA,catch.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat4.s$plotOM))){
  ests=lm(m.catch~plotMP-1,data=subset(catch_c_s,plotOM==levels(catch_c_s$plotOM)[i]))
  confs=confint(ests)
  catch_comp_d[match(levels(hcrdat4.s$plotOM)[i],catch_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  catch_comp_d[match(levels(hcrdat4.s$plotOM)[i],catch_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  catch_comp_d[match(levels(hcrdat4.s$plotOM)[i],catch_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  
}
catch_comp_d$scn.cols=scn_cols_static[as.numeric(catch_comp_d$scn)]

status_comp_d=expand.grid(mp=levels(hcrdat4.s$plotMP)[-1],scn=levels(hcrdat4.s$plotOM),d.ubm=NA,l95=NA,u95=NA,d.red=NA,l95r=NA,u95r=NA)
status_comp_d=as.data.frame(status_comp_d)
#  data.frame(scn=levels(hcrdat4.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat4.s$plotOM))){
  ests=lm(p.bm~plotMP-1,data=subset(bmdat_c_s,plotOM==levels(bmdat_c_s$plotOM)[i]))
  confs=confint(ests)
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
  ests2=lm(1-p.red~plotMP-1,data=subset(bmdat_c,plotOM==levels(bmdat_c$plotOM)[i]))
  confs2=confint(ests2)
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),6]=ests2$coefficients[2]/ests2$coefficients[1]-1
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),7]=confs2[2,1]/confs2[1,2]-1
  status_comp_d[match(levels(hcrdat4.s$plotOM)[i],status_comp_d$scn),8]=confs2[2,2]/confs2[1,1]-1
}

spawners_comp_d=expand.grid(mp=levels(hcrdat4.s$plotMP)[-1],scn=levels(hcrdat4.s$plotOM),d.spwnrs=NA,l95=NA,u95=NA)
spawners_comp_d=as.data.frame(spawners_comp_d)
#  data.frame(scn=levels(hcrdat4.s$plotOM),status.10y.s=NA,status.10y.rwa=NA,status.10y.rwa.l95=NA,status.10y.rwa.u95=NA,status.5y.s=NA,status.5y.s.l95=NA,status.5y.s.u95=NA,status.5y.rwa=NA,status.5y.rwa.l95=NA,status.5y.rwa.u95=NA)
for(i in 1:length(levels(hcrdat4.s$plotOM))){
  ests=lm(m.spawners~plotMP-1,data=subset(spawners_c_s,plotOM==levels(spawners_c_s$plotOM)[i]))
  confs=confint(ests)
  spawners_comp_d[match(levels(hcrdat4.s$plotOM)[i],spawners_comp_d$scn),3]=ests$coefficients[2]/ests$coefficients[1]-1
  spawners_comp_d[match(levels(hcrdat4.s$plotOM)[i],spawners_comp_d$scn),4]=confs[2,1]/confs[1,2]-1
  spawners_comp_d[match(levels(hcrdat4.s$plotOM)[i],spawners_comp_d$scn),5]=confs[2,2]/confs[1,1]-1
}


png(filename='outputs/figs/catch_conservation_tradeoff_UBM_ERtar0.9UMSY_0.1fixedBM.png',width=8,height=6,units='in',res=600)
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

png(filename='outputs/figs/catch_conservation_tradeoff_LBM_ERtar0.9UMSY_0.1fixedBM.png',width=8,height=6,units='in',res=600)
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

png(filename='outputs/figs/catch_conservation_tradeoff_spawner_ERtar1.0UMSY_0.3fixedBM.png',width=8,height=6,units='in',res=600)
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
  lines(c(spawners_comp_d$d.spwnrs[i]*100,catch_comp_d$d.catch[i]*100)~c(1.25+j,1.75+j),col=catch_comp_d$scn.cols[i],lwd=2)
  lines(c(catch_comp_d$l95[i]*100,catch_comp_d$u95[i]*100)~rep(1.75+j,2))
  lines(c(spawners_comp_d$l95[i]*100,status_comp_d$u95[i]*100)~rep(1.25+j,2))
  text(catch_comp_d$scn[i],x=par('usr')[2]*0.95,y=par('usr')[4]-par('usr')[4]*(0.13*i),col=catch_comp_d$scn.cols[i])
  points(d.catch[i]*100~c(1.75+j),data=catch_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
  points(d.spwnrs[i]*100~c(1.25+j),data=spawners_comp_d,cex=2,bg=catch_comp_d$scn.cols[i],pch=21)
}
mtext(side=1,at=1.75,'Mean annual catch',line=1.5)
mtext(side=1,at=1.25,'Mean annual spawners',line=1.5)
dev.off()





#composite plots####




#subset of scenarios:
st_agg<-ggplot(hcrdat1[hcrdat1$year<111,])+
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
#st_agg


st_agg_all<-ggplot(hcrdat1[hcrdat1$year>59&hcrdat1$year<111,])+
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

catch_plot_scu<-ggplot(catch_u,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  scale_fill_manual(values=c('blue4','darkgoldenrod','dodgerblue4','goldenrod2')) +
  scale_x_discrete(position = "top") +
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
  theme(legend.position="none",strip.background = element_blank(),strip.text = element_text(face = "bold"))


stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  ncol=3,greedy=F,align='hv',axis='bt',
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("./outputs/figs/composite_plot_status_catch_ERtar_0.9UMSY.png",plot=stclass,width=14,height=8.5)


st_agg<-ggplot(hcrdat2.s[hcrdat2.s$year<111,])+
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

#st_agg
st_agg_all<-ggplot(hcrdat2.s[hcrdat2.s$year>59&hcrdat2.s$year<111,])+
  geom_bar(aes(x=plotMP, fill=wsp.status),position = "fill")+
  scale_fill_manual(values = 
                      statusCols,name='Status')+
  scale_x_discrete(position = "top") +
  facet_grid(plotOM~.)+
  ylab('Proportion of years in status')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text.x = element_blank(),strip.text.y=element_blank())


catch_plot_scu<-ggplot(catch_b_s,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  scale_fill_manual(values=c('dodgerblue4','goldenrod2')) +
  scale_x_discrete(position = "top") +
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
  theme(legend.position="none",strip.background = element_blank(),strip.text = element_text(face = "bold"))



stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  align='hv',axis='bt',ncol=3,
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("./outputs/figs/composite_plot_status_catch_ERtar0.75_0.3BM.png",plot=stclass,width=14,height=8.5)

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



st_agg<-ggplot(hcrdat4[hcrdat4$year<111,])+
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

st_agg_all<-ggplot(hcrdat4[hcrdat4$year<111,])+
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


catch_plot_scu<-ggplot(catch_d,aes(x=plotMP,y=scale.ann.catch,fill=plotMP))+
  scale_fill_manual(values=c('blue4','darkgoldenrod','dodgerblue4','goldenrod2')) +
  geom_violin()+
  #  geom_boxplot(width=0.1, color="white", alpha=0.5) +
  stat_summary(fun= median,
               geom = "crossbar", 
               width = 0.98,
               position = position_dodge(width = .1),
  )+
  scale_x_discrete(position = "top") +
  facet_grid(plotOM~.)+
  ylab('Scaled Mean Annual Catch')+
  xlab('')+
  theme_bw(12)+
  theme(legend.position="none",strip.background = element_blank(),strip.text = element_text(face = "bold"))


stclass<-plot_grid(
  st_agg, st_agg_all+ theme(legend.position="none"),catch_plot_scu,
  align = "hv", axis = "bt",ncol=3,
  rel_widths = c(.5,.25,.25)
)

stclass
ggsave("./outputs/figs/composite_plot_status_catch_ERtar0.9UMSY_0.1fixedBM.png",plot=stclass,width=14,height=8.5)


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

#er trends####
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
ggsave("./outputs/figs/er_comp_ERtar_0.9UMSY.png",plot=er_plot_umsy,width=8,height=6)

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
ggsave("./outputs/figs/er_comp_ERtar0.75_0.1BM.png",plot=er_plot_bm,width=8,height=6)

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

er_plot_bmumsyf<-ggplot(meanER_d)+
  geom_line(aes(x=year-50,y=m.ER,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_line(aes(x=year-50,y=umsy,group=plotOM),linewidth=1.2,colour='green')+
  geom_ribbon(aes(ymin = l80.er, ymax = u80.er,x=year-50), alpha = 0.2)+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Exploitation Rate')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

er_plot_bmumsyf
ggsave("./outputs/figs/er_comp_ERtar0.9UMSY_0.1fixedBM.png",plot=er_plot_bmumsyf,width=8,height=6)

#spawner trends
sp_plot_umsy<-ggplot(spawners_u1)+
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
ggsave("./outputs/figs/spawners_timeseries_vs_smsy_ERtar_0.9UMSY.png",plot=sp_plot_umsy,width=8,height=6)

sp_plot_bm<-ggplot(spawners_b1)+
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
ggsave("./outputs/figs/spawners_timeseries_vs_smsy_ERtar0.65_0.1BM.png",plot=sp_plot_bm,width=8,height=6)

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

sp_plot_bmumsyf<-ggplot(spawners_d1)+
  geom_line(aes(x=year-50,y=m.spawners,color=plotMP,group=plotOM),linewidth=1.2)+
  geom_ribbon(aes(ymin = l80.s, ymax = u80.s,x=year-50), alpha = 0.2)+
  geom_line(aes(x=year-50,y=smsy,group=plotOM),linewidth=1.2,colour='green')+
  scale_color_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  facet_grid(plotOM~plotMP)+
  ylab('Spawner abundance')+
  xlab('Year of simulation')+
  theme_bw(12)+
  theme(legend.position="none")

sp_plot_bmumsyf
ggsave("./outputs/figs/spawners_timeseries_vs_smsy_ERtar0.9UMSY_0.1fixedBM.png",plot=sp_plot_bmumsyf,width=8,height=6)

#spawner violins####
sp_plot2_umsy<-ggplot(spawners_u2,aes(x=plotMP,y=log10(m.spawners),fill=plotMP))+
   coord_trans(y = "log10")+
  scale_fill_manual(values=c('dodgerblue4','goldenrod2','dodgerblue1','goldenrod1')) +
  ylim(4,5.5)+
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

sp_plot2_umsy
ggsave("./outputs/figs/spawners_violin_ERtar_0.9UMSY.png",plot=sp_plot2_umsy,width=8,height=12)

sp_plot2_bm<-ggplot(spawners_b2,aes(x=plotMP,y=log10(m.spawners),fill=plotMP))+
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

sp_plot2_bm
ggsave("outputs/figs/spawners_violin_ERtar0.65_0.1BM.png",plot=sp_plot2_bm,width=8,height=12)


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

sp_plot2_bmumsyf<-ggplot(spawners_d2,aes(x=plotMP,y=log10(m.spawners),fill=plotMP))+
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

sp_plot2_bmumsyf
ggsave("outputs/figs/spawners_violin_ERtar0.9UMSY_0.1fixedBM.png",plot=sp_plot2_bmumsyf,width=8,height=12)












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



