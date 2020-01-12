####PRELUDE
#Step 1. Download and Install all required packages.
list.of.packages <- c("RColorBrewer","classInt", "sp", "rgeos",
                      "sf","rgdal","geojsonio","tidyverse","raster",
                      "rasterVis","RStoolbox","plotly",
                      "htmlwidgets","magrittr","viridis", "tmap", "tmaptools", "shiny", "shinyjs")

#Checks whether any packages needed are not already downloaded&installed.
#Then installs all new packages.
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])] 
if(length(new.packages)) install.packages(new.packages)

#LIBRARY ALL PACKAGES
library(sf) #simplefeatures
library(sp)
library(rgeos)
library(rgdal)
library(geojsonio)
library(raster)
library(rasterVis)
library(RStoolbox)

library(plotly)
library(tidyverse)  #prerequisite for lots of packages
library(htmlwidgets)
library(RColorBrewer) #color scales
library(classInt)

library(magrittr) #Pipes - used for plotting
library(viridis)  #color scales
library(tmap)
library(tmaptools)
library(shiny)
library(shinyjs)
####BARCELONA
##STEP 1.
#Store GeoJSON of Barcelona file into simple features object. This is used later for cropping the raster.
#
barcelonaSF <- st_read("Barcelona/Barcelona.geojson")
#Change projection to match the projection of raster file.
newProjection <- "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs"
barcelonaSF <- st_transform(barcelonaSF, newProjection)
#
##STEP 2: RASTER PRE-PROCESSING
#Read in the meta file, and stack the raster bands using the meta file. Faster than reading in traditionally.
#Atmospheric processing, applying DOS2 correction from chavez (1996).
barcMTL <- readMeta("Barcelona/May_10_19/LC08_L1TP_198031_20190510_20190521_01_T1_MTL.txt")
espRS <- stackMeta("Barcelona/May_10_19/LC08_L1TP_198031_20190510_20190521_01_T1_MTL.txt")
#Atmospheric processing, applying DOS2 correction from chavez (1996).
#determine the SHV on the blue band, with a dark prop of 1% 
hazeDN <- estimateHaze(espRS, hazeBands = "B2_dn", darkProp = 0.01) 
#correct atmospherical scattering to deliver surface reflectance.
atmCorrectedRS <- radCor(espRS, metaData = barcMTL,
                         method = "costz",hazeValues = hazeDN, hazeBands = "B2_dn")

##Crop the Corrected Raster to the Barcelona extent and border. Write these out as rasters for storage purposes.
#Crop stack to barcelona shape.
barcelonaBound <- as(barcelonaSF, Class = "Spatial")
barcelonaExtent <- extent(barcelonaBound)
barcelonaRSCrop <- crop(atmCorrectedRS, barcelonaExtent)
barcelonaRSMask <- mask(barcelonaRSCrop, barcelonaBound)
#
##STEP 3: PROCESSING - SAVI, IBI, LST, UHI, UHS, UTFVI
#Create Index Function.
indexFunction <- function(band1, band2) {
  Index <- (band1 - band2) / (band1 + band2)
  return(Index)
}

##SAVI
#used to detect vegetation over NDVI because more sensitive for areas where vegetation density is low (i.e cities)
L = 0.5
barcelonaSAVI <- ((barcelonaRSMask$B5_sre - barcelonaRSMask$B4_sre)*(1+L))/
                    (barcelonaRSMask$B5_sre + barcelonaRSMask$B4_sre + L)

##NDBI
barcelonaNDBI <- indexFunction(barcelonaRSMask$B6_sre, barcelonaRSMask$B5_sre)

##MNDWI (Xu 2005)
barcelonaMNDWI <- indexFunction(barcelonaRSMask$B3_sre,barcelonaRSMask$B6_sre)

#We can now use these indexes to find specific types of land cover.

#Employ use of IBI to show built up land XU 2008
NDBI1 <- barcelonaNDBI + 1
SAVI1 <- barcelonaSAVI + 1
MNDWI1 <- barcelonaMNDWI + 1
barcelonaIBI <- ((NDBI1-((SAVI1+MNDWI1)/2))/
                   (NDBI1+((SAVI1+MNDWI1)/2)))


###LST
#As we have performed atmospheric correction, we already have the Brightness temperature values in bands 10.
#So no need to recalculate.
barcBT1 <- barcelonaRSMask$B10_bt

#calc fractional vegetation of each pixel in raster, 
barcFracV=(barcelonaSAVI-0.2/0.5-0.2)^2
emiss1=0.004*barcFracV+0.986

#p calculation
Boltzmann=1.38*10e-23
Plank=6.626*10e-34
c=2.998*10e8
p=Plank*(c/Boltzmann)

#define remaining varaibles
lambda=1.09e-5

#run the LST calculation
barcLST <-barcBT1/(1 +(lambda*barcBT1/p)*log(emiss1))
# convert LST from kelvin to celsius
barcLST<-barcLST-273.15


##DETERMINING UHI AND UHS (urban heat islands and urban heat spots)
sd1 <- cellStats(barcLST, stat = "sd")
mean1 <- cellStats(barcLST, stat = "mean")
barcUHI <- barcLST > (mean1 + (0.5*sd1))     #(Ma,Kuang & Huang 2010)
barcNonUHI <- barcLST <= (mean1 + (0.5*sd1)) #(Guha,Govil & Mukherjee 2017)
barcUHS <- barcLST > (mean1 + (2*sd1))

##UTFVI urban thermal field variance index for mapping heat stress (Zhang, 2006).
barcUTFVI <- (barcLST - mean1)/mean1



###LOS ANGELES - Same methodology as barcelona.
##STEP 1.
#Store GeoJSON of Los Angeles file into simple features object. This is used later for cropping the raster.
#
laSF <- st_read("LA/LosAngeles.geojson")

#
##STEP 2: RASTER PRE-PROCESSING
#Read in the meta file, and stack the raster bands using the meta file. Faster than reading in traditionally.

laMTL <- readMeta("LA/May_30_19/LC08_L1TP_041036_20190530_20190605_01_T1_MTL.txt")
laRS <- stackMeta("LA/May_30_19/LC08_L1TP_041036_20190530_20190605_01_T1_MTL.txt")
crs(laRS)
#Change projection of LA GeoJSON to match raster.
newProjection <- "+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs"
laSF <- st_transform(laSF, newProjection)
#Atmospheric processing, applying DOS2 correction from chavez (1996).
#determine the SHV on the blue band, with a dark prop of 1% 
hazeDN1 <- estimateHaze(laRS, hazeBands = "B2_dn", darkProp = 0.01) 
#correct atmospherical scattering to deliver surface reflectance.
laCorRS <- radCor(laRS, metaData = laMTL,
                         method = "costz",hazeValues = hazeDN1, hazeBands = "B2_dn")

##Crop the Corrected Raster to the LA extent and border. Write these out as rasters for storage purposes.
#Crop stack to LA shape.
laBound <- as(laSF, Class = "Spatial")
laExtent <- extent(laBound)
laRSCrop <- crop(laCorRS, laExtent)
laRSMask <- mask(laRSCrop, laBound)

#writeRaster(...)


##STEP 3: PROCESSING - SAVI, IBI, LST, UHI, UHS, UTFVI

##SAVI
L = 0.5 #Constant depends on amount of vegetation,
#0 = high veg density, 1 = low veg density.
#Generally taken as 0.5.

laSAVI <- ((laRSMask$B5_sre - laRSMask$B4_sre)*(1+L))/
                    (laRSMask$B5_sre + laRSMask$B4_sre + L)

##NDBI
laNDBI <- indexFunction(laRSMask$B6_sre, laRSMask$B5_sre)

##MNDWI - modified water index(Xu 2005)
laMNDWI <- indexFunction(laRSMask$B3_sre,laRSMask$B6_sre)

#We can now use these indexes to find specific types of land cover.

#Employ use of IBI to show built up land.
NDBI2 <- laNDBI + 1
SAVI2 <- laSAVI + 1
MNDWI2 <- laMNDWI + 1 #plus one to make all values positive

laIBI <- (NDBI2-(SAVI2+MNDWI2)/2)/
            (NDBI2+(SAVI2+MNDWI2)/2)

##LST
#As we have performed atmospheric correction, we already have the Brightness temperature values in bands 10.
#So no need to recalculate.
laBT <- laRSMask$B10_bt

#calc fractional vegetation of each pixel in raster, 
laFracV=(laSAVI-0.2/0.5-0.2)^2
emiss2=0.004*laFracV+0.986

#LST equation
#p calc - already stored variables but declare again to read easier.
Boltzmann=1.38*10e-23
Plank=6.626*10e-34
c=2.998*10e8
p=Plank*(c/Boltzmann)
lambda=1.09e-5
#run the LST calculation
laLST <-laBT/(1 +(lambda*laBT/p)*log(emiss2))
# convert LST from kelvin to celsius
laLST<-laLST-273.15

##DETERMINING UHI AND UHS
sd2 <- cellStats(laLST, stat = "sd")
mean2 <- cellStats(laLST, stat = "mean")
laUHI <- laLST > (mean2 + (0.5*sd2))      
laNonUHI <- laLST <= (mean2 + (0.5*sd2))  
laUHS <- laLST > (mean2 + (2*sd2))        

##UTFVI
laUTFVI <- (laLST - mean2)/mean2


####PLOTTING AND GRAPHS
#Regression of IBI/SAVI as my two variates with LST. 
#Summary statistics etc.

#Create map function so easier plot creation.
createMap <- function(ras, geoshp, colour,br, legTitle, Mtitle){
  tm_shape(ras) + 
    tm_raster(style = "fixed", title = Mtitle,
              
              palette = colour, midpoint = NA, breaks = br)+
    tm_legend(outside = FALSE, legend.position = c("right", "bottom"))+
    tm_shape(geoshp) + 
    tm_borders(alpha =  0.4, col = "black") + 
    tm_layout(inner.margins = c(0.1, 0.1, 0.1, 0.1),
              title = legTitle,
              title.size = 1.1,
              title.position = c("center", "top")) +
    tm_compass(type = "4star", size = 3, text.size = 1,position = c("left", "top")) +
    tm_scale_bar(position = c("left","bottom"))
  
}
createMap2 <- function(ras, geoshp, colour,br, legTitle, Mtitle){
  tm_shape(ras) + 
    tm_raster(style = "fixed", title = Mtitle,
              
              palette = colour, midpoint = NA, breaks = br)+
    tm_legend(outside = FALSE, legend.position = c("left", "bottom"))+
    tm_shape(geoshp) + 
    tm_borders(alpha =  0.4, col = "black") + 
    tm_layout(inner.margins = c(0.08, 0.1, 0.1, 0.1),
              title = legTitle,
              title.size = 1.1,
              title.position = c("center", "top")) +
    tm_compass(type = "4star", size = 3, text.size = 1,position = c("left", "top")) +
    tm_scale_bar(position = c("right","bottom"))
  
}

#Colour palettes
RedGreen <- c("#A93226","#E67E22","#F1C40F","#52BE80","#196F3D")
GreenRed <- c("#196F3D","#52BE80","#F1C40F","#E67E22","#A93226")
#Break points for SAVI threshold intervals
cuts = c(-1,0,0.1, 0.2, 0.3,1)

SAVImap1 <- createMap(barcelonaSAVI, barcelonaSF, RedGreen, cuts, "Barcelona - SAVI","SAVI")
SAVImap2 <- createMap2(laSAVI, laSF, RedGreen, cuts, "Los Angeles - SAVI", "SAVI")

SAVIplot <- tmap_arrange(SAVImap1, SAVImap2, ncol = 2)

IBImap1 <- createMap(barcelonaIBI, barcelonaSF, GreenRed, c(-1,-0.2,-0.1,0,0.1,1), "Barcelona - IBI", "IBI")
IBImap2 <- createMap2(laIBI, laSF, GreenRed, c(-1,-0.2,-0.1,0,0.1,1), "Los Angeles - IBI", "IBI")
IBIplot <- tmap_arrange(IBImap1, IBImap2, ncol=2)

LSTmap1 <- createMap(barcLST, barcelonaSF, GreenRed, c(14,22, 24, 25, 26,30, 36),"Barcelona - LST", "LST (Celsius)")
LSTmap2 <- createMap2(laLST, laSF, GreenRed, c(14,22,30,33,36,40,48), "Los Angeles - LST", "LST (Celsius)")
LSTmap <- tmap_arrange(LSTmap1, LSTmap2, ncol=2)
##

#Crop LST to uhi in order to work out midpoints for plotting

bUHIcrop <- mask(barcLST, barcUHI, maskvalue = 0, updatevalue = NA)
bNUHIcrop <- mask(barcLST,barcUHI, maskvalue = 1, updatevalue = NA)
lUHIcrop <- mask(laLST, laUHI, maskvalue = 0, updatevalue = NA)
lNUHIcrop <- mask(laLST, laUHI, maskvalue = 1, updatevalue = NA)
bUHScrop <- mask(barcLST, barcUHS, maskvalue = 0, updatevalue = NA)
lUHScrop <- mask(laLST, laUHS, maskvalue = 0, updatevalue = NA)

UHImap1 <- tm_shape(barcLST) +
  tm_raster(style = "cont", title = "LST (Celsius)", midpoint = 25.764, palette = "-RdYlGn",n = 9, contrast = c(0.4, 1)) +
  tm_legend(outside = FALSE, legend.position = c("right", "bottom"))+
  tm_shape(barcelonaSF) + 
  tm_borders(alpha =  0.4, col = "black") + 
  tm_layout(inner.margins = c(0.05, 0.02, 0.13, 0.1),
            title = "Barcelona LST\n UHI & Non-UHI",
            title.size = 1.1,
            title.position = c("center", "top")) +
  tm_compass(type = "4star", size = 2.5, text.size = 0.8,position = c("left", "top")) +
  tm_scale_bar(position = c("left","bottom"))


UHImap2 <- tm_shape(laLST) +
  tm_raster(style = "cont", title = "LST (Celsius)", midpoint = 35.51689, palette = "-RdYlGn",n = 9, contrast = c(0.4, 1)) +
  tm_legend(outside = FALSE, legend.position = c("left", "bottom"))+
  tm_shape(laSF) + 
  tm_borders(alpha =  0.4, col = "black") + 
  tm_layout(inner.margins = c(0.08, 0.1, 0.13, 0.1),
            title = "Los Angeles LST\n UHI & Non-UHI",
            title.size = 1.1,
            title.position = c("center", "top")) +
  tm_compass(type = "4star", size = 2.5, text.size = 0.8,position = c("left", "top")) +
  tm_scale_bar(position = c("right","bottom"))

UHIplot <- tmap_arrange(UHImap1, UHImap2, ncol = 2)

UTFVIbins <- c(-1,0,0.005,0.01,0.015,0.02,1)
UTFVImap1 <- createMap(barcUTFVI, barcelonaSF, GreenRed, UTFVIbins, "Barcelona - UTFVI", "UTFVI")
UTFVImap2 <- createMap2(laUTFVI, laSF, GreenRed, UTFVIbins, "Los Angeles - UTFVI", "UTFVI")
UTFVIplot <- tmap_arrange(UTFVImap1, UTFVImap2, ncol = 2)


UHSmap1 <- tm_shape(barcUHS) +
  tm_raster(style = 'cont',title = "Legend" ,palette = c("white", "red"), labels = c("Barcelona","Urban Hot Spots"))+
  tm_legend(outside = F, legend.position = c("right", "bottom"), frame = "black")+
  tm_shape(barcelonaSF) +
  tm_borders(alpha = 0.4, col = "black") +
  tm_layout(inner.margins = c(0.05, 0.02, 0.13, 0.1), 
            title = "Barcelona - Urban Heat Spots (UHS)",
            title.size = 1.1,
            title.position = c("center", "top")) +
  tm_compass(type = "4star", size = 2.5, text.size = 0.8,position = c("left", "top")) +
  tm_scale_bar(position = c("left","bottom"))

UHSmap2 <- tm_shape(laUHS) +
  tm_raster(style = 'cont',title = "Legend" ,palette = c("white", "red"), labels = c("Los Angeles","Urban Hot Spots"))+
  tm_legend(outside = F, legend.position = c("right", "bottom"), frame = "black")+
  tm_shape(laSF) +
  tm_borders(alpha = 0.4, col = "black") +
  tm_layout(inner.margins = c(0.05, 0.2, 0.13, 0.2), 
            title = "Los Angeles - Urban Heat Spots (UHS)",
            title.size = 1.1,
            title.position = c("center", "top")) +
  tm_compass(type = "4star", size = 2.5, text.size = 0.8,position = c("left", "top")) +
  tm_scale_bar(position = c("left","bottom"))

#calculate area of UTFVI
utfviDF1 <- as.data.frame(na.omit(values(barcUTFVI)))
utfviDF2 <- as.data.frame(na.omit(values(laUTFVI)))

bclass1 <- sum(utfviDF1[] < 0)
bclass2 <- sum(utfviDF1[] >= 0 & utfviDF1[] < 0.005)
bclass3 <- sum(utfviDF1[] >= 0.005 & utfviDF1[] < 0.01)
bclass4 <- sum(utfviDF1[] >= 0.01 & utfviDF1[] < 0.015)
bclass5 <- sum(utfviDF1[] >= 0.015 & utfviDF1[] < 0.02)
bclass6 <- sum(utfviDF1[] >= 0.02)

bTot <- bclass1+bclass2+bclass3+bclass4+bclass5+bclass6
b1 <- (bclass1/bTot)*100
b2 <- (bclass2/bTot)*100
b3 <- (bclass3/bTot)*100
b4 <- (bclass4/bTot)*100
b5 <- (bclass5/bTot)*100
b6 <- (bclass6/bTot)*100

lclass1 <- sum(utfviDF2[] < 0)
lclass2 <- sum(utfviDF2[] >= 0 & utfviDF2[] < 0.005)
lclass3 <- sum(utfviDF2[] >= 0.005 & utfviDF2[] < 0.01)
lclass4 <- sum(utfviDF2[] >= 0.01 & utfviDF2[] < 0.015)
lclass5 <- sum(utfviDF2[] >= 0.015 & utfviDF2[] < 0.02)
lclass6 <- sum(utfviDF2[] >= 0.02)

lTot <- lclass1+lclass2+lclass3+lclass4+lclass5+lclass6
l1 <- (lclass1/lTot)*100
l2 <- (lclass2/lTot)*100
l3 <- (lclass3/lTot)*100
l4 <- (lclass4/lTot)*100
l5 <- (lclass5/lTot)*100
l6 <- (lclass6/lTot)*100

####Linear Relationships 
#Stack all data to explore and turn to dataframe
barcStack <- stack(barcLST, barcelonaSAVI, barcelonaIBI)
barcDF <- as.data.frame(na.omit(values(barcStack)))
laStack <- stack(laLST, laSAVI, laIBI)
laDF <- as.data.frame(na.omit(values(laStack)))
#barcelona
bTempVeg <- lm(layer.2~layer.1, data = barcDF)
summary(bTempVeg)
bTempIBI <- lm(layer.3 ~ layer.1, data = barcDF)
summary(bTempIBI)
#LA
lTempVeg <- lm(layer.2~layer.1, data = laDF)
summary(lTempVeg)
lTempIBI <- lm(layer.3 ~ layer.1, data = laDF)
summary(lTempIBI)

#write out all plots
tmap_save(SAVIplot, filename = "Plots/SAVIplot.png")
tmap_save(UHIplot, filename = "Plots/UHIplot.png")
tmap_save(UHSmap1, filename = "Plots/UHSbarcelona.png")
tmap_save(UHSmap2, filename = "Plots/UHSla.png")
tmap_save(IBIplot, filename = "Plots/IBIplot.png")
tmap_save(LSTmap, filename = "Plots/LSTplot.png")
tmap_save(UTFVIplot, filename = "Plots/UTFVIplot.png")
