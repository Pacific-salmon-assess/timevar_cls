#edit simPars
library(stringr)

simpar<-read.csv("data/cls/SimPars.csv")
#change HCR numbers so that they are oredere  from most conservartionis to most aggressive

head(simpar)



#maphcr <- c(HCR1 = "HCR3", HCR2 = "HCR4", HCR3 = "HCR1", HCR4 = "HCR2")
#simpar$scenario <- str_replace(simpar$scenario, "HCR[1-4]", function(x) maphcr[x])
#simpar$nameMP <- str_replace(simpar$nameMP, "HCR[1-4]", function(x) maphcr[x])


#mapscn<-c(stationarylAR1="stationary0.4AR1", 
#        stationaryhAR1 = "stationary0.8AR1",
#        decLinearcap0.25 ="decLinearcap0.25x",                   
#        decLinearProd2to0.5 ="decLinearProd0.25x",
#        incLinearProd2to3="incLinearProd1.5x" , 
#        incLinearcap2= "incLinearcap2x",
#        regProd2to1 ="regProd0.5x",
#        regProd2to0.5= "regProd0.25x",
#        regCap0.25="regCap0.25x")
 
#keys_sorted <- names(mapscn)[order(-nchar(names(mapscn)))]
#pattern <- paste0("^(", paste(str_replace_all(keys_sorted, "([.])", "\\\\\\1"), collapse = "|"), ")")

#simpar <- simpar %>%
#  mutate(scenario = str_replace(scenario, pattern, function(x) mapscn[x]))

#simpar$nameOM <- mapscn[simpar$nameOM]

#write.csv(simpar,file="data/cls/SimPars.csv", row.names=FALSE)


head(simpar)

simpar$scenario[simpar$singleHCR=="retro"]<-paste0(simpar$scenario[simpar$singleHCR=="retro"],"_retro")
simpar$nameMP[simpar$singleHCR=="retro"]<-paste0(simpar$nameMP[simpar$singleHCR=="retro"],"_retro")


simpar$scenario[simpar$singleHCR=="forecast"]<-paste0(simpar$scenario[simpar$singleHCR=="forecast"],"_forecast")
simpar$nameMP[simpar$singleHCR=="forecast"]<-paste0(simpar$nameMP[simpar$singleHCR=="forecast"],"_forecast")



write.csv(simpar, file ="data/cls/SimPars.csv", row.names = FALSE)
head(simpar)

unique(simpar$singleHCR)