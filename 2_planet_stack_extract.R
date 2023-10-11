#Planet Stack and Extract 
#K Tiedeman
#adapted from BCI repository: script 2a_planet_timeseries_extract 
#11 October 2023

#Planet extract...

library(terra)
library(raster)
library(dismo)
library(stringr)
library(geodata)

###Function###

addVI <- function(x){
  x$NDVI <- (x$nir - x$red) / (x$nir + x$red)
  x$GNDVI <- (x$nir - x$green) / (x$nir + x$green)
  x$GRVI <- (x$green - x$red) / (x$green + x$red)
  x$GBVI <- (x$blue - x$green) / (x$blue + x$green)  
  x$EVI <- 2.5 * ((x$nir - x$red) / (x$nir + 6 * x$red - 7.5 * x$blue + 1))
  x$GCVI <- (x$nir/x$green)-1
  return(x)
}



# This can be commented out, useful for multiple users/computers 
this <- system('hostname', TRUE)
if (this == "kn-eas-c003") {
  dp <- "C:/Users/Hyaena/Documents/planet"
} else {
  dp <- "//10.126.19.90/EAS_ind/ktiedeman/data/bci" #
}


setwd(dp) ### set working directory

#The directory with your processed planet files - same as planet mask script
outMo <-'data/planet_processed/'


rastlist <- list.files(path = outMo, pattern='.tif$', 
                       all.files=TRUE, full.names=TRUE)

#import all raster files in folder using lapply - thank you https://stackoverflow.com/questions/52746936/how-to-efficiently-import-multiple-raster-tif-files-into-r
d <- lapply(rastlist, rast)

d


####

crown <- vect('data/Crowns/Crowns_2020_08_01_FullyMergedWithPlotData.shp')
plot(crown) #not sure what's happening some artifact or point outside the norm. - could clip this later to bci for ease 
crs(crown)
# 
# panama<-getData('GADM', country='PAN', level=0)
# panama <- geodata::gadm(country="PAN", level=0,path=tempdir()) #down until october 12 2023 

# 



crown <- project(crown, crs(d[[1]]))

plotRGB(d[[1]], r=3, g=2, b=1, stretch = "lin")
plot(crown, add=T)



bci_boundary <- vect("data/Barro_Colorado_Nature_Monument_Boundaries-shp/Barro_Colorado_Nature_Monument_Boundaries.shp")
bci_boundary <- project(bci_boundary, crs(d[[1]]))

crown <- terra::crop(crown, ext(bci_boundary))
plot(crown) # cool, removed whatever what happening outside of BCI for ease of plotting


#Decision point - 1. a centroid of each polygon or 2. all values inside of the polygon 


##Start with 1. Centroids 

y <- terra::centroids(crown) #,TRUE) #not a "true centroid", but will be inside of the polygon 
#(the points returned are guaranteed to be inside the polygons or on the lines, but they are not the true centroids.
#True centroids may be outside a polygon, for example when a polygon is "bean shaped")


extFund <- function(listrast, pts){
  
  al <- extract(listrast, pts,  xy = T, method = "simple", bind=T)
  x <- cbind(pts, al)
  
  
  dr <- sources(listrast)
  j <- str_extract(dr, "[^/]*$")[1] #unneccesary, but if they're in a subfolder helpful 
  #Extract just the date itself 
  x <- as.data.frame(x)
  x$date_Pl <- str_extract(j, "[^_]+") 
 
  return(x)
}



fgg <- lapply(d, extFund, pts = y) #fgg <- lapply(d, extFund, allpt = allpt)
fgrg <- do.call("rbind", fgg)


#Cleaning
fgrg <- fgrg[complete.cases(fgrg$blue),]

fgrg$day <-as.Date(as.character(fgrg$date_Pl), "%Y%m%d")

fgrg <- addVI(fgrg)


startday <- min(fgrg$date_Pl)
endday <- max(fgrg$date_Pl)

write.csv(fgrg, paste0("data/planet/planet_timeseries_nondipt_pts_", startday, "_",endday, ".csv" ))


#https://www.rdocumentation.org/packages/RStoolbox/versions/0.2.6/topics/spectralIndices


##################Polygons ################







ac <- extract(d[[1]], crown,  xy = T, method = "simple", bind=T)












#okay so what if we extract the polygons, maybe that's why we're not seeing a big signal??
dip <- shapefile("data/BCI_Dipteryx_maps/BCI_Dipteryx_maps/BCI_Dipteryx_Master_final.shp")
dip2 <- vect(dip)
dip2$pkuid <- as.numeric(dip2$pkuid)
#dex <- d[[2]]

ExtrPolyFun <- function(rex, dip2){
  
  polyid <- unique(unlist(as.numeric(dip2$pkuid)))
  polylist <- list() 
  
  for (s in 1:length(polyid)){
    lip <- dip2[dip2$pkuid == polyid[s],]
    sf <- data.frame(extract(rex, lip))
    names(sf) <- c("ID", "blue", "green", "red", "nir") 
    if (dim(sf)[1] != 0){
      sf$pkuid <- lip$pkuid
      sf$pixid <- paste0(sf$pkuid,1:nrow(sf))
      polylist[[s]] <- sf}
    
    
  } 
  
  plist <- do.call(rbind, polylist)
  
  dr <- sources(rex)
  j <- str_extract(dr, "[^/]*$")[1] 
  plist$date <- str_extract(j, "[^_]+")
  return(plist)  
  
}


ppe <- lapply(d, ExtrPolyFun, dip2=dip2)

ppe2 <- do.call(rbind, ppe)



ppe2$day <-as.Date(as.character(ppe2$date), "%Y%m%d")
ppe2$pixid <- as.numeric(as.character(ppe2$pixid))

#This should be moved to a different script, but for now, it stays
ppe2$NDVI <- (ppe2$nir - ppe2$red) / (ppe2$nir + ppe2$red)
ppe2$GNDVI <- (ppe2$nir - ppe2$green) / (ppe2$nir + ppe2$green)
ppe2$GRVI <- (ppe2$green - ppe2$red) / (ppe2$green + ppe2$red)
ppe2$GBVI <- (ppe2$blue - ppe2$green) / (ppe2$blue + ppe2$green)  
ppe2$EVI <- 2.5 * ((ppe2$nir - ppe2$red) / (ppe2$nir + 6 * ppe2$red - 7.5 * ppe2$blue + 1))
ppe2$GCVI <- (ppe2$nir/ppe2$green)-1

startday <- min(ppe2$date)
endday <- max(ppe2$date)
write.csv(ppe2, paste0("data/planet/planet_timeseries_dip_polygons_", startday, "_",endday, ".csv" ))


m <- read.csv("data/planet/planet_timeseries_dip_polygons_20180301_20200420.csv")
mh <- aggregate(cbind(NDVI,GNDVI,GRVI,GBVI,EVI, GCVI, blue, green,  red,  nir )~pkuid+day, ppe2, median)

mh$dipteryx <- 1

andere <- rbind(mh, subset(fgrg, select = -c(date, ID)))


write.csv(andere, paste0("data/planet/planet_timeseries_dip_polygonMedian_", startday, "_",endday, ".csv" ))


sg <- table(as.data.frame(andere$day, andere$pkuid))
