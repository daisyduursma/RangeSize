
#.libPaths('~/MyRlibs')

rm(list = ls())

#get source code
source('~/Google Drive/R_scripts/PhD/rangesize/fromStu/globals.R')

library(raster)

#setwd('~/Work/Range size/data')

#if (length(commandArgs(TRUE)) == 0) {

#load species presence
pres <- read.csv('~/Google Drive/PhD/BreedingRanges/tables/Breeding_XYalphahull_locations_2015-03-18.csv')
pres[pres==TRUE] <- 1
pres[prese==FALSE] <- 0

#loadRange area
data <- read.csv('~/Google Drive/PhD/BreedingRanges/tables/Breeding_alphahull_area_2015-03-18.csv')
data$species<-gsub(' ', '.', data$Taxon.scientific.name)

#load trait data

traits<-read.csv('~/Google Drive/PhD/birdTraits/Garnett paper/NSD-Data Descriptor/Australian_Bird_Data_Version_1_0.csv')
traits$species<-gsub(' ', '.', traits$Taxon.scientific.name)



  #output.filename <- '~/Google Drive/PhD/BreedingRanges/tables/Breeding_mean_range_size.csv'
  #method <- 'AreaKM2.25alpha_100kmres'
  #} else { 
  
#   # be sure to include relative paths with filenames  
#   presence.filename <- commandArgs(TRUE)[1]
#   range.size.filename <- commandArgs(TRUE)[2]
#   output.filename <- commandArgs(TRUE)[3]
#   method <- commandArgs(TRUE)[4]
#   
# }

# test with subset if you like
#m.presence <- matrix(m.presence[2001:3000])

# create empty data frame for storing computed range size data (mean, etc)
#mean.range.size <- data.frame(matrix(ncol = 7, nrow = nrow(presence.filename)))
#r_names <- c('n', 'min', 'max', 'sum', 'mean', 'sd', 'sd.n')

# iterate through species presence
# calculating range size stats for each cell
# note that R calculates stdev using (n-1), so also calculate using n

mean.range.size <- list()

for (i in 1:nrow(pres)) {
	message(i)
	if (sum(pres[i,7:ncol(pres)]) == 0) {
    next
	}else{
		# there are some species present in this cell
		# so do range size calcs
		
		# do this by creating temp variable x, which is a list of same length as the list of
		# species present in the cell, but containing range size data for each species
		# rather than species name
    summary<-data.frame(pres[i,c('x','y')])
    celldat<-pres[i,]
    celldat<-data.frame(t(celldat[7:ncol(celldat)]))
    colnames(celldat)<-'values'
    celldat$species<-rownames(celldat)
    celldat<-subset(celldat,values >0)
    x<-merge(celldat,data,by='species')
    x<-merge(x,traits,by='species')
		#x <- sapply(presence.filename[[i]], function(i) data[data$species == i, method])
    #x<-x[!sapply(x,is.null)]
    #message(x)
		
		# compute whatever you want from the data in x
		summary$n <- nrow(x)
		summary$RANGEmin <- min(x$AreaOfOccurance)
		summary$RANGEmax <- max(x$AreaOfOccurance)
		summary$RANGEsum <- sum(x$AreaOfOccurance)
		summary$RANGEmean <- mean(x$AreaOfOccurance)
		summary$RANGEsd <- sd(x$AreaOfOccurance)
		summary$RANGEsd.n <- sd(x$AreaOfOccurance) / sqrt(length(x$AreaOfOccurance)  / (length(x$AreaOfOccurance) - 1))
		
    lm(x$AreaOfOccurance,)
	
	
	}
	mean.range.size[[i]]<-summary

}

ENDdat<-do.call("rbind",mean.range.size)
ENDdat$log10<-log(ENDdat$mean)

aus<-raster('~/Google Drive/PhD/Data/Spatial/Climate/biovars/bio_1.asc',
            crs = "+proj=longlat +datum=WGS84")
aus<-trim(aggregate(aus,fact=10))
aus <- projectRaster(from = aus, crs = aea, res = 10000, method = 'bilinear')
aus100<-aggregate(aus,fact=10)#100km resolution ascii and points


dat<-rasterize(cbind(ENDdat$x,ENDdat$y),aus100,ENDdat$mean)
qqnorm(ENDdat$mean); qqline(ENDdat$mean, col = 2)
qqplot(ENDdat$mean, rt(300, df = 5))
#dat<-rasterize(cbind(presence.filename$x,presence.filename$y),aus100,presence.filename$Acanthagenys.rufogularis)



write.csv(ENDmean.range.size, output.filename, row.names = FALSE, na = '')  
  



