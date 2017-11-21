rm(list = ls())

#read in libraries
library(rgdal)
library(alphahull)
library(raster)
library(maptools)
library(rgeos)

#get source code
source('~/Google Drive/R_scripts/PhD/rangesize/fromStu/globals.R')


#turn on and off plots, plots of alphahulls and grid cells
#that are included in 100km grid cell presence
SAVE_PLOTS <- FALSE
xlim <- c(-2e+06, 2e+06)#set long limitis used in plotting
ylim <- c(-5e+06, -1e+06)#set lat limitis used in plotting

#read in an equal area mask of Australia,transform to equal area
shape <- readShapeSpatial('~/Google Drive/PhD/Data/Spatial/AustraliaMask.shp')
proj4string(shape) <- CRS('+init=epsg:4326')
aus_aea <- spTransform(shape, CRS(aea))

#10km resolution raster, points of Australia land,projected to equal area, 
#10km and 100km grid cell size, small values to work with alphahull
aus<-raster('~/Google Drive/PhD/Data/Spatial/Climate/biovars/bio_1.asc',
            crs = "+proj=longlat +datum=WGS84")
aus<-trim(aggregate(aus,fact=10))
aus <- projectRaster(from = aus, crs = aea, res = 10000, method = 'bilinear')
ausP<- rasterToPoints(aus)[,1:2]
ausPsmall<-ausP/1000000 #values small so that they work with alphahull
aus100<-aggregate(aus,fact=10)#100km resolution ascii and points
aus100mask<-aus100/aus100
aus100[] <- 1:ncell(aus100) #reasign values so every cell has own value
aus100<-aus100/aus100mask
aus100Points<-data.frame(rasterToPoints(aus100))

#read in observations
species.list <- read.csv('~/Google Drive/PhD/BreedingRanges/tables/BreedingSummary_2015-02-04.csv',
                stringsAsFactors = FALSE)
#species.list<-subset(species.list, UniqueLocsCount >= 10)
obs <- read.csv('~/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/UniqueLatLongYearMonthSourchEggobs_2015-03-17.csv',
                stringsAsFactors = FALSE)

#make empty list to feed results into
sub_obs <- list()
sp_pres <- list()

for (i in 1:nrow(species.list)) {
  #get species information
  species_name <- species.list[i, 'Taxon.scientific.name']
  message(paste(i, 'of', nrow(species.list), ': Species is', species_name))
  end.dat<-subset(species.list,Taxon.scientific.name==species_name)[,1:2]
  
  #get lat and lon and change to equal area
  species.data <- subset(obs,Scientific.Name == species_name)
  if(nrow(species.data)<=5)
    next
  lonlat <- SpatialPoints(data.frame(species.data[,c('lon', 'lat')]),
                          CRS("+init=epsg:4326"))
  lonlat <- spTransform(lonlat, CRS(aea))
  
  #use the center of the 10km gridcells to make the actual alphahulls
  # rasterize to grid of resolution as specified above in metres 
  obs.raster <- mask(rasterize(lonlat, aus, fbackground = 0,field=1),aus)
  obs.raster100 <- mask(rasterize(lonlat, aus100,  background = 0,field=1),aus100)
  p100<-rasterToPoints(obs.raster100, fun=function(x) {x == 1})[,c('x','y')]
  cellvalue100<-extract(aus100,p100,method='simple')
  
  #count of occupancy at differnt scales
  end.dat$observationCount<-nrow(species.data)
  end.dat$observationCount10km<-cellStats(obs.raster,stat='sum')
  end.dat$observationCount100km<-cellStats(obs.raster100,stat='sum')
  #area of occupancy100km
  end.dat$AreaSpeciesObservations100km<- end.dat$observationCount100km*10000
  
   
  #find locations and convert to use in alpahull
  p <- rasterToPoints(obs.raster, fun=function(x) {x == 1})
  p <- p[, 1:2]
  multiplier <- 1000000
  p <- p/multiplier
  p <- p + matrix(rnorm(length(p), sd = 10^-4), nc = 2)  
  
  #calc alpha hull area(km2) based on centroids
  a25 <- ahull(p, alpha = (0.25 * (1000000 / multiplier)))
  a30 <- ahull(p, alpha = (0.3 * (1000000 / multiplier)))
  a40 <- ahull(p, alpha = (0.4 * (1000000 / multiplier)))
  a60 <- ahull(p, alpha = (0.6 * (1000000 / multiplier)))
  a80 <- ahull(p, alpha = (0.8 * (1000000 / multiplier)))
  EOO <- ahull(p, alpha = 100)
  
  #find out what 10km resolution grid cells the species alpha hulls are in
  TFpres.25<-apply(ausPsmall, 1, function(x) inahull(a25, x))
  TFpres.30<-apply(ausPsmall, 1, function(x) inahull(a30, x))
  TFpres.40<-apply(ausPsmall, 1, function(x) inahull(a40, x))
  TFpres.60<-apply(ausPsmall, 1, function(x) inahull(a60, x))
  TFpres.80<-apply(ausPsmall, 1, function(x) inahull(a80, x))
  TFpres.E00<-apply(ausPsmall, 1, function(x) inahull(EOO, x))
  
  #area of 10km cells 
  end.dat$AreaKM2.25alpha_10kmres <- length(TFpres.25[TFpres.25 %in% TRUE])*100
  end.dat$AreaKM2.30alpha_10kmres <- length(TFpres.30[TFpres.30 %in% TRUE])*100
  end.dat$AreaKM2.40alpha_10kmres <- length(TFpres.40[TFpres.40 %in% TRUE])*100
  end.dat$AreaKM2.60alpha_10kmres <- length(TFpres.60[TFpres.60 %in% TRUE])*100
  end.dat$AreaKM2.80alpha_10kmres <- length(TFpres.80[TFpres.80 %in% TRUE])*100
  end.dat$AreaKM2.E00_10kmres <- length(TFpres.E00[TFpres.E00 %in% TRUE])*100
  
  
  #area of 100km cells,find which 100km grid cells 10km grid cells belong to
  cellvalue<-extract(aus100,ausP,method='simple')#get 100km cell values
  
  dat25<-as.data.frame(subset(cbind(ausP,cellvalue,TFpres.25),TFpres.25==1))#combine with TRUE FAlSE presence data
  cells25<-unique(as.numeric(dat25$cellvalue))#find unique 100km cell values
  dat30<-as.data.frame(subset(cbind(ausP,cellvalue,TFpres.30),TFpres.30==1))#combine with TRUE FAlSE presence data
  cells30<-unique(as.numeric(dat30$cellvalue))#find unique 100km cell values
  
  #set threshold for level of occurance within 10000km2 grid cells
  threshold <- 30
 
#alpha = .25
 n.cells <- table(dat25$cellvalue)
 pres25<-names(n.cells)[n.cells >= threshold]
 #calculate area
 ifelse(length(pres25) == 0,
      end.dat$AreaKM2.25alpha_100kmres<-0,
      end.dat$AreaKM2.25alpha_100kmres<-length(pres25) *10000)

#alpha = .30
 n.cells <- table(dat30$cellvalue)
 pres30<-names(n.cells)[n.cells >= threshold]
 ifelse(length(pres30) == 0,
        end.dat$AreaKM2.30alpha_100kmres<-0,
        end.dat$AreaKM2.30alpha_100kmres<-length(pres30) *10000)


if(end.dat$observationCount<=3){
      end.dat$AreaOfOccurance<-end.dat$AreaSpeciesObservations100km
      end.dat$AreaOfOccuranceDerived<-'AreaSpeciesObservations100km'
  }else{ 
      if (end.dat$AreaSpeciesObservations100km>end.dat$AreaKM2.30alpha_100kmres){
            end.dat$AreaOfOccurance<-end.dat$AreaSpeciesObservations100km
            end.dat$AreaOfOccuranceDerived<-'AreaSpeciesObservations100km(largerArea)'
      }else{
        end.dat$AreaOfOccurance<-end.dat$AreaKM2.30alpha_100kmres
        end.dat$AreaOfOccuranceDerived<-'AreaKM2.30alpha_100kmres'
      }
} 
  #make TRUE and FALSE values of species presence for every grid cell
  cellPres<-data.frame(aus100Points$layer %in% pres30) #make column of TF values if species occur in 
  colnames(cellPres)<-species_name


  
  
  
  #loop to make jpgs of alpha-hulls and cells that are included in 100km raster
  if (SAVE_PLOTS) {
    #for plotting
    smallP <- rasterToPoints(obs.raster, fun=function(x) {x == 1})[,1:2]
    multiplier <- 1
    smallP <- smallP/multiplier
    smallP <- smallP + matrix(rnorm(length(p), sd = 10^-4), nc = 2) 
    ahull25 <- ahull(smallP, alpha = (0.25 * (1000000 / multiplier)))
    ahull30 <- ahull(smallP, alpha = (0.3 * (1000000 / multiplier)))
    ahullEOO <- ashape(smallP, alpha = (100 * (1000000 / multiplier)))
    
    jpeg(filename = paste0('/Users/daisy/Google Drive/PhD/BreedingRanges/figures/alpha_hull_plots/',
                          tolower(gsub(' ', '_', species_name)), '.jpg'),
                          width = 1000,
                          height = 600, 
                          units = "px", 
                          quality = 100)
    par( mfrow = c( 2, 3 ), mar=c(0,0,0,0),oma = c( 0, 0, 2, 0 ) )
    
    #extent of occurrence
    plot(aus_aea, xlim = xlim, ylim = ylim, border = "grey")
    plot(ahullEOO, wpoints = FALSE, xlim = xlim, ylim = ylim, 
         col = "purple", xlab = '', ylab = '',add=TRUE)
    points(smallP, pch=20,col="black",cex=.5)
    mtext("Extent of Occurrence",
          cex= 1,line = -3, font= 1)
    
    #area of occurrence
    plot(aus_aea, xlim = xlim, ylim = ylim, border = "grey")
    plot(ahull30, wpoints = FALSE, xlim = xlim, ylim = ylim, 
        col = "blue", xlab = '', ylab = '',add=TRUE)
    plot(ahull25, wpoints = FALSE, xlim = xlim, ylim = ylim, 
        col = "orange", xlab = '', ylab = '',add=TRUE)
    points(p, pch=20,col="black",cex=.5)
    mtext("orange: alpha=.25, blue: alpha=.30",
        cex= 1,line = -3, font= 1)
   
    #plot 100km cells where observations are
    plot(trim(obs.raster100),axes=FALSE,box=FALSE,legend=FALSE)
    mtext("Obs - 100km res.",
        cex= 1,line = -3, font= 1)
    
    #plot (alpha25 with size)
    r<-aus100
    r[r%in%c(pres25)] <- 1
    r[!r%in%c(1)] <- 0
    r<-mask(r,aus100)
    plot(trim(r),axes=FALSE,box=FALSE,legend=FALSE)
    mtext("alpha .25",cex= 1,line = -3, font= 1)
    
    #plot (alpha30 with size)
    r<-aus100
    r[r%in%c(pres30)] <- 1
    r[!r%in%c(1)] <- 0
    r<-mask(r,aus100)
    plot(trim(r),axes=FALSE,box=FALSE,legend=FALSE)
    mtext("alpha .30",cex= 1,line = -3, font= 1) 
    
   title(main = paste(species_name),outer=TRUE,cex.main= 2, font.main = 3)
    dev.off()      
  }
  
  #put information to keep into lists
  sub_obs[[i]]<-end.dat #summary data for species
  sp_pres [[i]]<-cellPres #alpha hull presence at 100 km
  
}

#put data togther in dataframe
breeding_ahull<-do.call("rbind",sub_obs)
spPres100km<-do.call("cbind",sp_pres)

spPres100km<-cbind(aus100Points,spPres100km)
spPres100km[spPres100km==TRUE] <- 1





# z <- sub_obs[[1]]
# for(i in 2:300){
#   x <- sub_obs[[i]]
#   names(x) <- names(z)
#   z <- rbind(z, x)
#   message(i)
# }

#make map of one of the suitable areas
#xy<-spPres100km[,c('x','y')]
#vals<-spPres100km$"Acanthagenys rufogularis"
#r0 <- rasterize(xy,aus100,vals)
#plot(r0)

#write out file
write.csv(z, paste0('~/Google Drive/PhD/BreedingRanges/tables/Breeding_alphahull_area_',
                              as.Date(Sys.time()),'.csv'),row.names=FALSE)
write.csv(spPres100km, paste0('~/Google Drive/PhD/BreedingRanges/tables/Breeding_XYalphahull_locations_',
                              as.Date(Sys.time()),'.csv'),row.names=FALSE)


