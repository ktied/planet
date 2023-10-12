#Planet Stack and Extract 
#K Tiedeman
#adapted from BCI repository: script 2a_planet_timeseries_extract 
#11 October 2023

#Planet extract: this script stacks raster files into a list, and has two helper functions to extract the data either from centroids or points or from the polygons of the crowns 

library(terra)
library(raster)
library(dismo)
library(stringr)
library(geodata)

###Function to add Vegetation indicies### Can add more as desired 
#https://www.rdocumentation.org/packages/RStoolbox/versions/0.2.6/topics/spectralIndices
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
# Keeping as a list allows different extents hopefully 
d


####Read crown data #####

crown <- vect('data/Crowns/Crowns_2020_08_01_FullyMergedWithPlotData.shp')
plot(crown) #not sure what's happening some artifact or point outside the norm. - could clip this later to bci for ease 
crs(crown)

crown <- project(crown, crs(d[[1]]))

plotRGB(d[[1]], r=3, g=2, b=1, stretch = "lin")
plot(crown, add=T)

#Clipping to BCI boundary 
bci_boundary <- vect("data/Barro_Colorado_Nature_Monument_Boundaries-shp/Barro_Colorado_Nature_Monument_Boundaries.shp")
bci_boundary <- project(bci_boundary, crs(d[[1]]))

crown <- terra::crop(crown, ext(bci_boundary))
plot(crown) # cool, removed whatever what happening outside of BCI for ease of plotting

#############################################
#############################################
#Decision point - 1. a centroid of each polygon or 2. all values inside of the polygon 


###### Centroids ##############

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


fgg <- lapply(d, extFund, pts = y) 
fgrg <- do.call("rbind", fgg)


#Cleaning, in case of cloud cover 
fgrg <- fgrg[complete.cases(fgrg$blue),]

#add date of planet image from the string
fgrg$day <-as.Date(as.character(fgrg$date_Pl), "%Y%m%d")

#add vegetation indicies
fgrg <- addVI(fgrg)

#Write file to csv (could also be rds)
startday <- min(fgrg$date_Pl)
endday <- max(fgrg$date_Pl)

write.csv(fgrg, paste0("data/extracted_data/planet_timeseries_50HA_pts_", startday, "_",endday, ".csv" ))



##########################################
##################Polygons ################




#If you have a single crown or want to define the function 
ac <- extract(d[[1]], crown, fun=median, xy = T, cells=T, method = "simple", bind=T)

#Id field, I'm choosing stemID 

# 
# 
# 
# 
# 
# 
# #okay so what if we extract the polygons, maybe that's why we're not seeing a big signal??
# dip <- shapefile("data/BCI_Dipteryx_maps/BCI_Dipteryx_maps/BCI_Dipteryx_Master_final.shp")
# dip2 <- vect(dip)
# dip2$pkuid <- as.numeric(dip2$pkuid)
# #dex <- d[[2]]
# 
# dip2 <- crown
# id <- "stemID"
# rex = d[[1]]

ExtrPolyFun <- function(rex, dip2, id){
  
  polyid <- unique(unlist(dip2[[id]]))
  polylist <- list() 
  
  for (s in 1:length(polyid)){
    lip <- dip2[dip2[[id]] == polyid[s],]
    sfg <- data.frame(extract(rex, lip, xy=T, cells = T, method="simple"))

    if (dim(sfg)[1] != 0){
      sfg[[id]] <- lip[[id]][1,]
      sfg <- merge(sfg, lip, by = id) #This isn't particularly elegant, and I'm sure a function is out there that does exactly this already 
      polylist[[s]] <- sfg }
    
    
  } 
  
  plist <- do.call(rbind, polylist)
  
  dr <- sources(rex)
  j <- str_extract(dr, "[^/]*$")[1] 
  plist$date_Pl <- str_extract(j, "[^_]+")
  return(plist)  
  
}

#run the function and rbind results 
ppe <- lapply(d, ExtrPolyFun, dip2=crown, id="stemID")
ppe2 <- do.call(rbind, ppe)


#add date of planet image from the string
ppe2$day <-as.Date(as.character(ppe2$date_Pl), "%Y%m%d")

#clean, remove rows (pixels) that were NA due to cloud cover 
ppe2 <- ppe2[complete.cases(ppe2$blue),]

#add vegetation indicies 
ppe2 <- addVI(ppe2)

#get the start and end of your raster list and write file to csv 
startday <- min(ppe2$date)
endday <- max(ppe2$date)
write.csv(ppe2, paste0("data/planet/planet_timeseries_50HA_polygons_", startday, "_",endday, ".csv" ))


