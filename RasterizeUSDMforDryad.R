## This is an R script to rasterize the USDM to a gridded support.
##
## 1. Libraries
#################################
library(sf)
library(raster) ## install version 3.5.29 from source.
library(ggplot2)
library(RColorBrewer)
library(maps)
library(fields)
#################################

## 2. Code to load one USDM file, set things up, make adjacency matrix
#################################
usdm=st_read("USDM/USDM_20140729_M/USDM_20140729.shp") ## Chosen due to severity of drought
head(usdm)

usa <- maps::map("state", fill=TRUE)
usa2<-map2SpatialPolygons(usa,
                          IDs=sapply(strsplit(usa$names, ":"), "[", 1L),
                          proj4string=CRS("+proj=longlat +datum=WGS84"))
orng=brewer.pal(6, "Oranges")

## Rasterize this to just a regular rectangular grid
rst_template <- raster(ncols = 120, nrows = 50,
                       crs = projection(usdm),
                       ext = extent(c(-125,-65,25,50)))
raster_data <- rasterize(usdm, rst_template)
plot(raster_data)

## Quiltplot with state lines
lon=seq(-124.75,-65.25,by=0.5)
lat=seq(49.75,25.25,by=-0.5)
quilt.plot(expand.grid(lon, lat),c(raster_data[]), nx=120, ny=50, xlab="Longitude", ylab="Latitude", main="USDM", col=orng)
map("state",lwd=1,add=T)

## Need to intersect raster with outer shape of US
us.grid=rasterize(usa2,rst_template)
us.grid=(us.grid>0)-2 ## Make a 0-map for valid locations

name1=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
        "AA","BB","CC","DD","EE","FF","GG","HH","II","JJ","KK","LL","MM","NN","OO","PP","QQ","RR","SS","TT","UU","VV","WW","XX")
name2=seq(1,120,by=1)

## We need to redact S2, S105, EE97, KK88, and RR85 because these locations are over oceans
## but the shapefile usa2 doesn't redact them.
us.grid[which(name1=="S"),which(name2==2)] <- NA
us.grid[which(name1=="S"),which(name2==105)] <- NA
us.grid[which(name1=="EE"),which(name2==97)] <- NA
us.grid[which(name1=="KK"),which(name2==88)] <- NA
us.grid[which(name1=="RR"),which(name2==85)] <- NA

plot(us.grid) ## Should be -1 if in the contiguous US

## Can restrict plot using mask function
usdmr <- mask(crop(raster_data, us.grid), us.grid)
plot(usdmr, col=orng)
map("state",lwd=1,add=T)

## Differentiate between NA (outside US) and level 0.
dm=usdmr
for (i in 1:50){
  print(i)
  for (j in 1:120) {
    if(is.na(us.grid[i,j])=="FALSE") {
      if(is.na(usdmr[i,j])=="TRUE") {dm[i,j]=0} else {}
    }
  }
}

## NA means outside US.  0 means no response.  D0 through D4 are coded 1 to 5.
plot(dm, col=orng)

## Need response vector (0,1,2, etc.) for valid cells
dm[is.na(us.grid)==FALSE]

## Adjancancy Matrix
wh=which(is.na(us.grid[])=="FALSE") ## wh count relevant cells on the full extent
coord=expand.grid(lon, lat)
coord.use=coord[wh,]

plot(coord.use, cex=0.2, pch=16)

name3a=expand.grid(name2,name1)
name3=data.frame(name3a[,2],name3a[,1])

name4=character()
for (i in 1:6000){
  name4[i]=paste(name3[i,1],name3[i,2],sep="")
}
head(name4)
labels=name4[is.na(us.grid[])=="FALSE"]

dz=as.matrix(dist(coord.use))

A = matrix(0,dim(dz)[1],dim(dz)[1])
rownames(A)=wh
colnames(A)=wh

for (i in 1:dim(dz)[1]){
  print(i)
  tmp=which((dz[i,] > 0.01) & (dz[i,]<0.71))
  A[i,tmp]=1
}

rownames(A)=labels
colnames(A)=labels

write.csv(A, "A.csv")
############################


## 4. Read in all years of USDM shapefiles, rasterize, save .csv
#############################
yearseq = seq(2000,2022,by=1)
filenames=character()
setwd("C:/Users/erhardrj/Desktop/Everything/Research/Drought/Full Dataset/USDM")
for (i in 1:length(yearseq)){
  dirs = list.dirs(paste(yearseq[i],"_USDM_M", sep = ""))
  for (j in 2:length(dirs)){
    newfilenames <- Sys.glob(paste(dirs[j],"/*.shp",sep=""))
    filenames=append(filenames, newfilenames)
  }
}

filenames

dm=list()
drought=numeric()
for (q in 1:length(filenames)){
  usdm=st_read(filenames[q], quiet=TRUE)
  cat(q)
  ##   usdm=dplyr::select(usdm, "DM")
  raster_data <- rasterize(usdm, rst_template)
  usdmr <- mask(crop(raster_data, us.grid), us.grid)

  dm[[q]]=usdmr
  dm[[q]][(is.na(us.grid[])=="FALSE") & (is.na(usdmr[])=="TRUE")]=0
  for (i in 1:50){
    print(q+i/100)
    for (j in 1:120) {
      if(is.na(us.grid[i,j])=="FALSE") {
        if(is.na(usdmr[i,j])=="TRUE") {dm[[q]][i,j]=0} else {}
      }
    }
  }
  tmp=dm[[q]][]
  tmp=tmp[is.na(tmp)=="FALSE"]
  tmp2=tmp
  for (m in 1:length(data.frame(usdm)[,1])){
    tmp2[tmp2==data.frame(usdm)[m,1]]=paste("D",data.frame(usdm)[m,2],sep="")
  }
  drought=append(drought,tmp2)
}
drought

time=character()
for (i in 1:length(filenames)){
  cat(i, sep="\n")
  tmp=dm[[i]][]
  tmp=tmp[is.na(tmp)=="FALSE"]
  hold=rep(substring(filenames[i],18,25),length(tmp))
  time=append(time,hold)
}
time

name5=rep(labels,length(filenames))
coord1=rep(coord.use[,1],length(filenames))
coord2=rep(coord.use[,2],length(filenames))

data=data.frame(time,name5,coord1,coord2,drought[1:length(time)])
colnames(data)=c("time","grid","lon","lat","drought")
table(data$drought)
save.image(file="USDMRasterized.rda")

write.csv(data, file="usdmfinal.csv")
##########################################


