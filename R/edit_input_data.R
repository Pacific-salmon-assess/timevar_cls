#edit simPars


simpar<-read.csv("data/cls/SimPars.csv")


simpar$scenario[simpar$singleHCR=="retro"]<-paste0(simpar$scenario[simpar$singleHCR=="retro"],"_retro")
simpar$nameMP[simpar$singleHCR=="retro"]<-paste0(simpar$nameMP[simpar$singleHCR=="retro"],"_retro")


simpar$scenario[simpar$singleHCR=="forecast"]<-paste0(simpar$scenario[simpar$singleHCR=="forecast"],"_forecast")
simpar$nameMP[simpar$singleHCR=="forecast"]<-paste0(simpar$nameMP[simpar$singleHCR=="forecast"],"_forecast")



write.csv(simpar, file ="data/cls/SimPars.csv", row.names = FALSE)
head(simpar)

unique(simpar$singleHCR)