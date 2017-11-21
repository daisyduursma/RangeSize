
rm(list = ls())

dat<-read.csv('/Users/daisy/Google Drive/PhD/BreedingRanges/tables/Breeding_alphahull_area_2015-03-03.csv')

#Plots for sensitivty analysis
  #Goal: Understand the effect of alpha on alpha-hull area

#rearange data for plotting in mutliple colours
DF<-data.frame(Area_km2=c(
          dat$AreaKM2.E00_10kmres,
          dat$AreaKM2.alpha_30_100km_25cellper100_AND_observations,
          dat$Aalpha_30_100km25cellper100_OR_smallArea_5cellper100,
          dat$Aalpha_25_100km25cellper100_OR_smallArea_5cellper100),
          Obs_100km_Areakm2=rep(dat$speciesOccurence100km,4),alpha=(rep(c(6:3),each=nrow(dat))))

DF$Area_km2<-DF$Area_km2/10000
DF$Obs_100km_Areakm2<-DF$Obs_100km_Areakm2/10000

DF$PercentChangeAOO_EOO<-(DF$Area_km2-DF$Obs_100km_Areakm2)/DF$Area_km2*100


# Plot
plot(Area_km2 ~ Obs_100km_Areakm2, pch=16, col=alpha, data = DF,cex=.7)
legend('topleft',c('e00','a30_obs', 'a30', 'a25'),col=6:3,pch=16, cex=.7)
abline(a = 0, b = 1, col = "gray60")



plot (Area_km2~log10(PercentChangeAOO_EOO), pch=16, col=alpha, data = DF,cex=.7)
legend('topleft',c('e00','a30_obs', 'a30', 'a25'),col=6:3,pch=16, cex=.7)






DF2<-subset(DF,alpha==5 | alpha == 6 | alpha == 7)
plot (PercentChangeAOO_EOO~Areakm, pch=16, col=alpha, data = DF2,cex=.5)
legend('topright',paste('Alpha',c(2,2.5,3),sep=" "),col=7:5,pch=8,cex=1)

#fit line for alpha = 2
a20<-subset(DF2, alpha==7)
loess_fit <- loess(PercentChangeAOO_EOO ~ Areakm, a20)
j <- order(a20$Areakm)
lines(a20$Areakm[j],loess_fit$fitted[j], col = 7,lwd=7)

#fit line for alpha = 2.5
a25<-subset(DF2, alpha==6)
loess_fit <- loess(PercentChangeAOO_EOO ~ Areakm, a25)
j <- order(a25$Areakm)
lines(a25$Areakm[j],loess_fit$fitted[j], col = 6,lwd=7)

#fit line for alpha = 3
a30<-subset(DF2, alpha==5)
loess_fit <- loess(PercentChangeAOO_EOO ~ Areakm, a30)
j <- order(a30$Areakm)
lines(a30$Areakm[j],loess_fit$fitted[j], col = 5,lwd=7)




