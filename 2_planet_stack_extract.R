#Planet Stack and Extract 

#Planet extract...
library(terra)
library(dismo)
library(stringr)

this <- system('hostname', TRUE)
if (this == "DESKTOP-J9EEJ0L") {
  dp <- "C:/Users/Kate/Dropbox/R_data/bci"
} else {
  dp <- "//10.126.19.90/EAS_ind/ktiedeman/data/bci" 
}

setwd(dp) ### set directory

allpt <- vect( "data/intermediate_crowns_points/DiptPresAbse1921pts_inBoth.shp")
#read in bci imagery 

s <- vect("//10.126.19.90/EAS_ind/ktiedeman/data/bci/data/Tree2006AstTab.shp")


##Let's move away from being a hack. 

rastlist <- list.files(path = "data/imagery/planet_processed/", pattern='.tif$', 
                       all.files=TRUE, full.names=TRUE)

#import all raster files in folder using lapply - thank you https://stackoverflow.com/questions/52746936/how-to-efficiently-import-multiple-raster-tif-files-into-r
d <- lapply(rastlist, rast)


extFund <- function(listrast, allpt){
  
  al <- extract(listrast, allpt) 
  al <- data.frame(pkuid = allpt$pkuid, dipteryx= allpt$dipteryx, al)
  
  dr <- sources(listrast)
  j <- str_extract(dr, "[^/]*$")[1] #unneccesary, but if they're in a subfolder helpful 
  #Extract just the date itself 
  #al <- as.data.frame(al)
  al$date <- str_extract(j, "[^_]+")
  names(al) <- c("pkuid","dipteryx","ID", "blue", "green", "red", "nir", "date")
  
  return(al)
}



extFundND <- function(listrast, allpt){
  
  al <- extract(listrast, allpt) 
  al <- data.frame(Species = allpt$Species, ID= allpt$id_1, al)
  
  dr <- sources(listrast)
  j <- str_extract(dr, "[^/]*$")[1] #unneccesary, but if they're in a subfolder helpful 
  #Extract just the date itself 
  #al <- as.data.frame(al)
  al$date <- str_extract(j, "[^_]+")
  names(al) <- c("Species","ID", "n", "blue", "green", "red", "nir", "date")
  
  return(al)
}



fgg <- lapply(d, extFundND, allpt = s) #fgg <- lapply(d, extFund, allpt = allpt)


fgrg <- do.call("rbind", fgg)

fgrg <- fgrg[complete.cases(fgrg$blue),]

fgrg$day <-as.Date(as.character(fgrg$date), "%Y%m%d")
fgrg$pkuid <- as.numeric(as.character(fgrg$pkuid))

#This should be moved to a different script, but for now, it stays
fgrg$NDVI <- (fgrg$nir - fgrg$red) / (fgrg$nir + fgrg$red)
fgrg$GNDVI <- (fgrg$nir - fgrg$green) / (fgrg$nir + fgrg$green)
fgrg$GRVI <- (fgrg$green - fgrg$red) / (fgrg$green + fgrg$red)
fgrg$GBVI <- (fgrg$blue - fgrg$green) / (fgrg$blue + fgrg$green)  
fgrg$EVI <- 2.5 * ((fgrg$nir - fgrg$red) / (fgrg$nir + 6 * fgrg$red - 7.5 * fgrg$blue + 1))
fgrg$GCVI <- (fgrg$nir/fgrg$green)-1

startday <- min(fgrg$date)
endday <- max(fgrg$date)
write.csv(fgrg, paste0("data/planet/planet_timeseries_nondipt_pts_", startday, "_",endday, ".csv" ))


#https://www.rdocumentation.org/packages/RStoolbox/versions/0.2.6/topics/spectralIndices


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
