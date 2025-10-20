#install samsim
#remotes::install_github("Pacific-salmon-assess/samEst", force=TRUE)
#remotes::install_github("Pacific-salmon-assess/samSim", ref="sbccnk", force=TRUE)

library(samEst)
#library(samSim)
library(ggplot2)
library(dplyr)
library(cowplot)

png(filename='figs_AFS/stepwiseHCR1.png',width=6,height=4,units='in',res=600)
plot(1:12,1:12,bty='l',type='n',xaxt='n',yaxt='n',xlab='',
     ylab='',ylim=c(1,13), cex.lab=1.4)
title(ylab="Exploitation Rate", xlab="Spawner Abundance", line=1, cex.lab=1.5, family="Calibri Light")
lines(x=c(4,4),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(x=c(8,8),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(c(0,4),c(1,1),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5),)
lines(c(4,8),c(4,4),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
lines(c(8,13),c(8,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
lines(c(4,4),c(1,4),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
lines(c(8,8),c(4,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
text(2,1,labels="10% ER",col="black",adj=c(.5,-1), cex=1.4, font.lab=2)
text(6,4,labels=expression(paste("50% of ",U[MSY])),col="black",adj=c(.5,-1), font.lab=2, cex=1.4)
text(10,8,labels=expression(paste(U[MSY])),col="black",adj=c(.5,-1), cex=1.4)
text(2,10,labels="RED",col="red",adj=c(.5,-1), cex=1.4, font.lab=2)
text(6,10,labels="AMBER",col="darkgoldenrod2",adj=c(.5,-1), cex=1.4, font.lab=2)
text(10,10,labels="GREEN",col="chartreuse4",adj=c(.5,-1), cex=1.4, font.lab=3)
text(4,12,labels=expression(paste(S[gen])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
text(8,12,labels=expression(paste(S[MSY])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
dev.off()



png(filename='figs_AFS/highERamberHCR2.png',width=6,height=4,units='in',res=600)
plot(1:12,1:12,bty='l',type='n',xaxt='n',yaxt='n',xlab='',
     ylab='',ylim=c(1,13), cex.lab=1.4)
title(ylab="Exploitation Rate", xlab="Spawner Abundance", line=1, cex.lab=1.5, family="Calibri Light")
lines(x=c(4,4),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(x=c(8,8),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(c(0,4),c(1,1),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5),)
lines(c(4,13),c(8,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
lines(c(4,4),c(1,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
text(2,1,labels="10% ER",col="black",adj=c(.5,-1), cex=1.4, font.lab=2)
text(10,8,labels=expression(paste(U[MSY])),col="black",adj=c(.5,-1), cex=1.4)
text(2,10,labels="RED",col="red",adj=c(.5,-1), cex=1.4, font.lab=2)
text(6,10,labels="AMBER",col="darkgoldenrod2",adj=c(.5,-1), cex=1.4, font.lab=2)
text(10,10,labels="GREEN",col="chartreuse4",adj=c(.5,-1), cex=1.4, font.lab=3)
text(4,12,labels=expression(paste(S[gen])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
text(8,12,labels=expression(paste(S[MSY])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
dev.off()





png(filename='figs_AFS/lowERamberHCR3.png',width=6,height=4,units='in',res=600)
plot(1:12,1:12,bty='l',type='n',xaxt='n',yaxt='n',xlab='',
     ylab='',ylim=c(1,13), cex.lab=1.4)
title(ylab="Exploitation Rate", xlab="Spawner Abundance", line=1, cex.lab=1.5, family="Calibri Light")
lines(x=c(4,4),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(x=c(8,8),y=c(0,12),lty=3,col=adjustcolor('gray30',alpha.f = 0.5),lwd=2)
lines(c(0,8),c(1,1),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5),)
lines(c(8,13),c(8,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
lines(c(8,8),c(1,8),lty=2, lwd=4,col=adjustcolor('darkred',alpha.f = 0.5))
text(2,1,labels="10% ER",col="black",adj=c(.5,-1), cex=1.4, font.lab=2)
text(10,8,labels=expression(paste(U[MSY])),col="black",adj=c(.5,-1), cex=1.4)
text(2,10,labels="RED",col="red",adj=c(.5,-1), cex=1.4, font.lab=2)
text(6,10,labels="AMBER",col="darkgoldenrod2",adj=c(.5,-1), cex=1.4, font.lab=2)
text(10,10,labels="GREEN",col="chartreuse4",adj=c(.5,-1), cex=1.4, font.lab=3)
text(4,12,labels=expression(paste(S[gen])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
text(8,12,labels=expression(paste(S[MSY])),col="black",adj=c(.5,-0), cex=1.4, font.lab=3)
dev.off()