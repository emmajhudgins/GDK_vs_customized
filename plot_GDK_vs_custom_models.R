## Plotting functions to compare GDK and customized models ##
## Written as part of Hudgins et al. "Comparing generalized to customized models for United States forest pest spread" in prep. J Ecol.

## Written by Emma J. Hudgins
## PhD candidate, McGill University
## Montreal, QC, Canada
## emma.hudgins@mail.mcgill.ca

## Example input is given for the best-fitting HWA customized and GDK models, but additional output can be generated using the other scripts in this folder

rm(list=ls()) 
library(sp)
library(ggplot2)
library(maptools)
library(maps)
library(wesanderson)

spp=49
nm_spp<-rep(0,64)
nm_spp[c(49,51,54)]<-c("HWA", "GM", "BBD")
i=7
cc=1
step_num=1
#setwd("~/Documents/GitHub/GDK_vs_customized/") # set to your own directory

#Read in Data
data<-read.csv('countydatanorm_march.csv', stringsAsFactors = FALSE) # spatial data
data2<-read.csv('spdat_clean_gdk.csv', stringsAsFactors = FALSE) # species data, see Hudgins et al. corrigendum for information on why ALB (spp=3) cannot be accurately fit
prez<-read.csv('prez_clean_gdk.csv') # invasible range by pest species (based on FIA)

if (spp==49)
{
  YEAR=50
  total_YEAR=55
  pfull2<-read.csv('hwa_3yrs.csv')[,1:3]
  prez2<-read.csv('hwa_5yr_fitting_july.csv')[,3:12]
}
if (spp==51)
{
  YEAR=135
  total_YEAR=140
  pfull2<-read.csv('gm_3yrs.csv')[,1:3]
  prez2<-read.csv('gm_5yr_fitting_july.csv')[,3:29]
}
if (spp==54)
{
  YEAR=110
  total_YEAR=115
  pfull2<-read.csv('bbd_3yrs.csv')[,1:3]
  prez2<-read.csv('bbd_5yr_fitting_july.csv')[,3:24]
}

for (m in 1:3)
{
  pfull2[,m]<-c(intersect(prez[which(prez[,spp]!=0),spp], pfull2[which(pfull2[,m]!=0),m]), rep(0,(3372-length(intersect(prez[which(prez[,spp]!=0),spp], pfull2[which(pfull2[,m]!=0),m])))))
}
for (sppp in 1:ncol(prez2))
{prez2[,sppp]<-c(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp]), rep(0,(3372-length(intersect(prez[which(prez[,spp]!=0),spp], prez2[which(prez2[,sppp]!=0),sppp])))))}


m<-SpatialPointsDataFrame(coords=cbind(data$X_coord, data$Y_coord), data=data)
USA_merged<-map('usa', fill=TRUE, plot=F)
IDs <- sapply(strsplit(USA_merged$names, ":"), function(x) x[1])
sp_map_usa <- map2SpatialPolygons(USA_merged, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
transform_usa<-spTransform(sp_map_usa, CRS("+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) # simulation grid is in United States Equidistant Conic projection

Presences<-list()
title<-c("Customized", "GDK (uncor)", "GDK (cor)")
palette<-wes_palette("Rushmore", 5, type="discrete")
mat <- matrix(c(1,2,3,3,4,4),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = mat,heights = c(0.45,0.45,0.1))
par(mar=c(2,1,1,1))


Presences[[1]]<-read.csv(paste("presences_customized",step_num,i,cc,"csv", sep="."))
Presences[[2]]<-read.csv(paste("presences_GDK",nm_spp[spp], "_u__fit_.csv", sep="_"))
if (spp==49)
{
Presences[[3]]<-read.csv(paste("presences_GDK",nm_spp[spp],"s_ic_c_fit_.csv", sep="_"))
}
if (spp!=49)
{
  Presences[[3]]<-read.csv(paste(spp_nm, "presences_GDK_s_ic_fit_.csv", sep="_"))
}


for (i in 1:3)
{
  Pfull<-Presences[[i]]
  plot(transform_usa, lwd=0.2, main=title[i])
  points(cbind(data$X_coord[prez[,spp]], data$Y_coord[prez[,spp]]), pch=15, cex=0.5, col=alpha(palette[2], 0.5))
  true_locations<-cbind(data$X_coord[pfull2[,3]], data$Y_coord[pfull2[,3]])
  points(true_locations, pch=15, cex=0.5, col=alpha(palette[4],0.5))
  pred_locations<-cbind(data$X_coord[Pfull[,(total_YEAR/5)]], data$Y_coord[Pfull[,(total_YEAR/5)]])
  points(pred_locations, pch=15, cex=0.5, col=alpha(palette[5], 0.35))
  plot(transform_usa,add=TRUE, lwd=2)
}
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("topleft",legend = c("Known Host Distribution","Known 2005 Spread", "Predicted 2005 Spread"), 
       col=c(palette[c(2,4,5)]), horiz = TRUE,xpd=TRUE,cex=1, pch=15, bty="n", x.intersp=0.5, pt.cex=2,text.width=c(0,0.25,0.25))
