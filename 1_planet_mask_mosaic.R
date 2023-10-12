#K Tiedeman
#adapted from BCI repository: 0_planet_masking 
#4 October 2023


#The following code uses the Surface Reflectance (SR) Planet Images,
#Masks the images with the udm2 clip, and then mosaics if there are multiple images for a single day 
# 

# Set up: 
#input folder
indir <- paste0("data/input/")


#output folder
outdir <-'data/planet_masked/'
#mosaic folder (taking multiple images per day and mosaicking them)\

#install.packages("gdalUtils", repos="http://R-Forge.R-project.org")

#packages 
library(terra)
library(dismo)
library(randomForest)
#library(gdalUtils)
library(stringr)


# This can be commented out, useful for multiple users/computers 
this <- system('hostname', TRUE)
if (this == "kn-eas-c003") {
  dp <- "C:/Users/Hyaena/Documents/planet"
} else {
  dp <- "//10.126.19.90/EAS_ind/ktiedeman/data/bci" #
}

setwd(dp) ### set directory
#https://developers.planet.com/docs/data/psscene4band/#available-asset-type

###FIRST, A Single image 

sep <- rast("data/input/PSScene/20230911_153434_97_2484_3B_AnalyticMS_SR_clip.tif")
#1 blue, 2 green, 3 red , 4 NIR 
udm2 <- rast("data/input/PSScene/20230911_153434_97_2484_3B_udm2_clip.tif")
plot(udm2)

#Udm2 is the quality control - cloud, could shadow and the haze 
#https://medium.com/planet-stories/finding-the-usable-data-in-planet-imagery-64582d21430e
#8 layers https://developers.planet.com/docs/data/udm-2/

# Band 1	Clear map	[0, 1]	0: not clear, 1: clear
# Band 2	Snow map	[0, 1]	0: no snow or ice, 1: snow or ice
# Band 3	Shadow map	[0, 1]	0: no shadow, 1: shadow
# Band 4	Light haze map	[0, 1]	0: no light haze, 1: light haze
# Band 5	Heavy haze map	[0, 1]	0: no heavy haze, 1: heavy haze
# Band 6	Cloud map	[0, 1]	0: no cloud, 1: cloud
# Band 7	Confidence map	[0-100]	percentage value: per-pixel algorithmic confidence in classification
# Band 8	Unusable pixels	--	Equivalent to the UDM asset: see Planet's Imagery Specification for complete details

#From here, the idea 
idealmask <- udm2[[1]] == 1 & udm2[[2]] == 0 & udm2[[3]]==0 & udm2[[4]]==0 & udm2[[5]]==0 & udm2[[6]] == 0& udm2[[7]] > 70
idealmask[idealmask==0] <- NA
#This space for testing. For loops that should be functions below:::##
SepMasked <- mask(sep, idealmask)
plot(SepMasked)

plotRGB(SepMasked, r=3, g=2, b=1, stretch = "lin")


#####For use with multiple images per day, or multiple downloads#####
####Mask with planet's udm2 mask in a for loop bc of creative ending....#####

#Repeats from above 
#input folder
indir <- paste0("data/input/")

#output folder
outdir <-'data/planet_masked/'

#Find all the files
dn <- list.files(path = indir, pattern =  c("*AnalyticMS_clip.tif$|*AnalyticMS_SR_clip.tif$"), recursive = T) # 

dn




MaskFun <- function(x, indir , outdir, CF = 0 ){ #CF is the confidence 0-100
  
  j <- str_extract(x, "[^/]*$")
  outfile <- sub("_3B_AnalyticMS_SR_clip.tif.*", "", j)   
  #Paths 
  pat <- sub("_3B_AnalyticMS_SR_clip.tif.*", "", x)    
  
  
  outname <- paste0(outdir, outfile,'_ps4_3B_Analytic_MS_masked', '.tif')
  print(outname)
  
#  if(!file.exists(outname)){ #This is to avoid rerunning 
    rawrast <- rast(paste0(indir, x))
    print(paste0(indir,x))
    
    if(file.exists(paste0(indir, pat, "_3B_udm2_clip.tif"))){
      maskj <- rast(paste0(indir, pat, "_3B_udm2_clip.tif"))
      
      print(paste0(indir, pat, "_3B_udm2_clip.tif"))
      
      ###Should be a function, but here you can define your confidence level 
      idealmask <- maskj[[1]] == 1 & maskj[[2]] == 0 & maskj[[3]]==0 & maskj[[4]]==0 & maskj[[5]]==0 & maskj[[6]] == 0 & maskj[[7]] > CF
      
    } else {
      
      maskj <- rast(paste0(indir, pat, "_3B_AnalyticMS_DN_udm_clip.tif"))
      
      print(paste0(indir, pat, "_3B_AnalyticMS_DN_udm_clip.tif"))
      
      idealmask <- maskj < 5 
      
    }
    
    
    maskit <- mask(rawrast, idealmask, maskvalue=0, filename = outname, overwrite=T)
    return(maskit)
  #}
  
}

#One individual one 
dn1 <- MaskFun(dn[1], indir, outdir, CF = 80)


#On the list as a whole 
s <- lapply(dn, MaskFun, indir, outdir, CF = 80)

##################################################################################

#Now you have cloud masked rasters, and for each date you may want to mosaic 

#THis worked  well to mosaic the rasters - however GdalUtils has depreciated,
#I'm now using vrt as suggested by Robert: #https://stackoverflow.com/questions/67169266/error-in-mosaic-of-rasters-from-different-extent-using-terra-package-in-r

outMo <-'data/planet_processed/'


dn <- list.files(path = outdir, pattern = "*Analytic_MS_masked.tif$", recursive = T)
j <- str_extract(dn, "[^/]*$") #unneccesary, but if they're in a subfolder helpful 

#Extract just the date itself 
pa <- str_extract(j, "[^_]+")
# find the unique dates
#dt <- unlist(lapply(strsplit(dn, "_"),"[[",1))
dt <- unlist(unique(pa))




for (i in (1:length(dt))){
  
  outname <- paste0(outMo, dt[i],'_ps4_3B_Analytic_MS_processed', '.tif')
  print(outname)
  #if(!file.exists(outname)){  #This piece can be uncommented if you're rerunning and don't want to redo work 

    rr <- dn[grep(dt[i], dn)]
    rr <- paste0(outdir, rr)
   # v <- vrt(c(rr)) #https://stackoverflow.com/questions/67169266/error-in-mosaic-of-rasters-from-different-extent-using-terra-package-in-r
    gg <- lapply(rr, rast)
    rsrc <- sprc(gg)
    m <- mosaic(rsrc)
    writeRaster(m, outname, overwrite =T)
    #mosaic_rasters(gdalfile = rr, dst_dataset = outname, of="GTiff", verbose=TRUE)#, co = list(list("COMPRESS" = "DEFLATE")))
  #}
  
}


