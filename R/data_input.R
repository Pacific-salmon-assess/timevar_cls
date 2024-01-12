#=================================================
#Data population
#=================================================

library(samEst)
library(samSim)
library(ggplot2)
?cuParexample

#cvER = 0.228 -- need to check if this is in line with onbserved variability


?simParexample


prop<-c(.1,.4,.3,.2)
truedf<-data.frame(Var2=1:4,prop=prop)
dfage<-reshape2::melt(t(replicate(1000,ppnAgeErr(prop,.5,runif(4, 0.0001, 0.9999)))))


ggplot(dfage)+
geom_line(aes(x=Var2,y=value,group=Var1),alpha=.3)+
geom_line(data=truedf,aes(x=Var2,y=prop),color="blue",linewidth=2)


1.49-0.65


0.65/1.49


1.5*0.44