format_dat=function(hcrdat){

  hcrdat$aboveLowerassess<-"correct below lower BM"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==1]<-"correct above lower BM"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==0& hcrdat$aboveLowerObsBM==1]<-"wrong optimistic"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==0]<-"wrong pessimistic"
  hcrdat$aboveUpperassess<-"correct below upper BM"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==1]<-"correct above upper BM"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==0& hcrdat$aboveUpperObsBM==1]<-"wrong optimistic"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==0]<-"wrong pessimistic"
  hcrdat$HCRtype='UMSY'
  
  #==================================================
  #red amber green
  
  hcrdat$wsp.status<-NA
  hcrdat$wsp.status[hcrdat$aboveLowerBM==1&
                      hcrdat$aboveUpperBM==1]<-"green"
  hcrdat$wsp.status[hcrdat$aboveLowerBM==1&
                      hcrdat$aboveUpperBM==0]<-"amber"
  hcrdat$wsp.status[hcrdat$aboveLowerBM==0&
                      hcrdat$aboveUpperBM==0]<-"red"
  
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
  
  hcrdat$wsp.status<- factor(hcrdat$wsp.status,levels=c('red','amber','green'))
  hcrdat$status_agg<-factor(hcrdat$status_agg,levels=c("red",
                                                       "amber",
                                                       "green",
                                                       "pessimistic",
                                                       "optimistic"))
  
  
  
  #these palletes are color bling friendlish. -- second one is a bit harder to see in black and white. Change with caution.
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
  
  #reduced by one scenario - include prod 2-0.5 and 2-1
  
  
  #hcrdat2=subset(hcrdat,nameOM!='decLinearProd1.5to0.5')
  
  #unique(hcrdat2$nameOM)
  hcrdat$plotOM<-dplyr::recode_factor(factor(hcrdat$nameOM),"stationaryhAR1" = "AR1 high",
                                      "stationarylAR1" = "AR1 low",  
                                      "decLinearcap0.50"="cap * 0.5", 
                                      "decLinearcap0.25"="cap * 0.25",
                                      "incLinearcap1.5"='cap * 1.5',
                                      "decLinearProd2to1.5"="prod 2 -> 1.5",     
                                      "decLinearProd2to1"="prod 2 -> 1",
                                      "decLinearProd2to0.5"="prod 2 -> 0.5",
                                      "incLinearProd2to3"="prod 2 -> 3")
  
  hcrdat$plotMP<-dplyr::recode_factor(factor(hcrdat$nameMP),"10yr_autocorr_constER" = "Stationary (10y)",
                                      "10yr_rwa_constER" = "Time-varying (10y)",  
                                      "5yr_autocorr_constER" = "Stationary (5y)",
                                      "5yr_rwa_constER" = "Time-varying (5y)",
                                      "10yr_autocorr_adaptER" = "Stationary (10y)",
                                      "10yr_rwa_adaptER" = "Time-varying (10y)",  
                                      "10yr_both_adaptER" = "Mixed (10y)",
                                      "5yr_autocorr_adaptER" = "Stationary (5y)",
                                      "5yr_rwa_adaptER" = "Time-varying (5y)",
                                      "5yr_both_adaptER" = "Mixed (5y)")
  
  
  return(hcrdat)
}

format_dat2=function(hcrdat){
  
  hcrdat$aboveLowerassess<-"correct below lower BM"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==1]<-"correct above lower BM"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==0& hcrdat$aboveLowerObsBM==1]<-"wrong optimistic"
  hcrdat$aboveLowerassess[hcrdat$aboveLowerBM==1& hcrdat$aboveLowerObsBM==0]<-"wrong pessimistic"
  hcrdat$aboveUpperassess<-"correct below upper BM"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==1]<-"correct above upper BM"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==0& hcrdat$aboveUpperObsBM==1]<-"wrong optimistic"
  hcrdat$aboveUpperassess[hcrdat$aboveUpperBM==1& hcrdat$aboveUpperObsBM==0]<-"wrong pessimistic"
  hcrdat$HCRtype='UMSY'
  
  #==================================================
  #red amber green
  
  hcrdat$wsp.status<-NA
  hcrdat$wsp.status[hcrdat$aboveLowerBM==1&
                      hcrdat$aboveUpperBM==1]<-"green"
  hcrdat$wsp.status[hcrdat$aboveLowerBM==1&
                      hcrdat$aboveUpperBM==0]<-"amber"
  hcrdat$wsp.status[hcrdat$aboveLowerBM==0&
                      hcrdat$aboveUpperBM==0]<-"red"
  
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
  
  hcrdat$wsp.status<- factor(hcrdat$wsp.status,levels=c('red','amber','green'))
  hcrdat$status_agg<-factor(hcrdat$status_agg,levels=c("red",
                                                       "amber",
                                                       "green",
                                                       "pessimistic",
                                                       "optimistic"))
  
  
  
  #these palletes are color bling friendlish. -- second one is a bit harder to see in black and white. Change with caution.
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

  
  
  return(hcrdat)
}