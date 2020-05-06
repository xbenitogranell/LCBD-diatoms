#load libraries used
library(raster)
library(tidyr)


#Geodiversity predictors
## code from Antonelli et al 2018 Nature Geosciences
mountain <- read.table("data/DB2017.txt", header=T, sep="\t") 
mountain$HERMratio<-mountain$Avg_HERM02/mountain$Avg_HERM46

#Check and exclude data Avg_USP >0 and RELIEF_MN != -9999
summary(mountain$Avg_USP)
mountain_a<-subset(mountain, mountain$Avg_USP >0)
str(mountain_a)
summary(mountain_a$Countpnts_t)

mountain_b<-subset(mountain_a, mountain_a$RELIEF_MN != -9999)
str(mountain_b)

#Assign new name for dataset
mountain0<-mountain_b
str(mountain0)

#Make subset where number of geological datapoints t >= 3
mountain3<-subset(mountain0, mountain0$Countpnts_t>2)

#Ensure that long-term erosion rate is not zero
mountain3b<-subset(mountain3, mountain3$Avg_HERMER>0)

#Assign final dataset for analysis
mountain_t3<-mountain3b

# subset Andes
mountain_andes<-subset(mountain_t3, mountain_t3$REGION_2017 == 'ANDES')
write.csv(mountain_andes, "data/mountain_andes.csv")

#Dataset for multi-predictor models
mountain_andes_subset <- read.csv("data/mountain_andes.csv", row.names = 1) %>%
  transform(Longitude=Centroid_X, Latitude=Centroid_Y) %>%
  dplyr::select(Longitude, Latitude, grid_bio1, grid_bio7, grid_bio12, grid_bio15,
         RELIEF_MN, Avg_HERMER, Avg_USP, soil_var, Rugg_avg, Rugg_range)



#Correlations among predictor variables
round(cor(mountain_andes_subset, method="spearman"), 2)   # correlation matrix for appendix


#panel correlation plots to assess data distribution
panel.hist <- function(x, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(usr[1:2], 0, 1.5) )     
  h <- hist(x, plot = FALSE)     
  breaks <- h$breaks; nB <- length(breaks)     
  y <- h$counts; y <- y/max(y)     
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(0, 1, 0, 1))     
  r <- abs(cor(x, y, use = "complete"))   
  txt <- format(c(r, 0.123456789), digits = digits)[1]     
  txt <- paste0(prefix, txt)     
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
  text(0.5, 0.5, txt, cex = cex.cor * r) }


# Plot
pairs(mountain_andes_subset, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 



#### HYDE DB ####

##download HYDE 3.2 database from ftp server
times <- as.character(seq(from=1700, to=2000, by=10)) # Define all needed time steps

#population density
ZIPfilenames <- as.list(paste(rep(c("lower","baseline", "upper"), time=10),
                              "/zip/",times, "AD_pop.zip", sep="")) # create the list of filenames

#landuse
ZIPfilenames <- as.list(paste(rep(c("lower","baseline", "upper"), time=10),
                              "/zip/",times, "AD_lu.zip", sep="")) # create the list of filenames

# anthromes
ZIPfilenames <- as.list(paste("anthromes/","zip/",times, "AD_anthromes.zip", sep="")) # create the list of filenames

names(ZIPfilenames)<- ZIPfilenames

#for all information on the HYDE 3.2 data used, see  ftp://ftp.pbl.nl/hyde/hyde3.2/readme_release_HYDE3.2.1.txt
hydeRasters <-
  plyr::llply(ZIPfilenames, function(x){    #iterate through filename list
    temp<- tempfile() #create temporary file
    download.file(url=as.character(paste("ftp://ftp.pbl.nl/hyde/hyde3.2/",x, sep="")), destfile=temp)    #download file from ftp using my filename list
    ## [[i]] list xxxx_lu 
    raster(unzip(temp)[[1]], crs=CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0"))  #unzip and directly import raster to avoid stocking too much data
    #raster(unzip(temp), crs=CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0"))  #unzip and directly import raster to avoid stocking too much data
    unlink(temp)
  }, .progress="text")

hydeStack <- raster::stack(hydeRasters)
save(hydeRasters, file="hyde/pop.count.Holocene.hyde3.31.rasters.Rdata")
save(hydeRasters, file="hyde/landuse.Holocene.hyde3.31.rasters.Rdata")
save(hydeRasters, file="hyde/anthromes.Holocene.hyde3.31.rasters.Rdata")

nlayers(hydeStack)
plot(hydeStack)


writeRaster(hydeStack, filename="hyde/humanpopAD17002000.Holocene.hyde3.21.rasters.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(hydeStack, filename="hyde/landuseAD17002000.Holocene.hyde3.21.rasters.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


coorSites <- lonLatSites
coorSites$region <- NULL
coorSites$ID <- rownames(coorSites)
colnames(coorSites) <- c("Latitude", "Longitude", "sites")
coorSites <- coorSites[,c(2,1,3)] #reorder columns
#write.csv(coorSites, "coorSitesNAndes.csv")
coordinates(coorSites) <- c('Longitude', 'Latitude')

#coordinates(coorSites) <- ~Longitude+Latitude

#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")


load("hyde/pop.count.Holocene.hyde3.31.rasters.Rdata")
load("hyde/landuse.Holocene.hyde3.31.rasters.Rdata") #this is cropland area
load("hyde/anthromes.Holocene.hyde3.31.rasters.Rdata")
hydeStack <- stack(hydeRasters)

sitesHyde <- data.frame(coorSites@data[["sites"]],raster::extract(hydeStack, coorSites))

  colnames(sitesHyde) <- gsub(".zip", "", colnames(sitesHyde))
  colnames(sitesHyde) <- gsub("_lu", "", colnames(sitesHyde))
  
  colnames(sitesHyde) <- gsub("AD_pop", "", colnames(sitesHyde) )
  colnames(sitesHyde) <- gsub("anthromes.", "", colnames(sitesHyde) )
  
  
sitesHydeL <- sitesHyde %>% gather(key="name", "pop", 2:32)
sitesHydeL <- sitesHyde %>% gather(key="name", "anthromes", 2:32)
sitesHydeL <- sitesHyde %>% gather(key="name", "cropland", 2:32)

sitesHydeL <- data.frame(sitesHydeL, plyr::ldply(strsplit(sitesHydeL$name,".",fixed=TRUE)))

sitesHydeL[,4] <- as.factor(sitesHydeL[,4])
sitesHydeL[,5] <- as.numeric(sitesHydeL[,5])
colnames(sitesHydeL)[c(1,4,5)] <- c("site", "scenario", "time")

sitesHydeL[,3] <- as.factor(sitesHydeL[,3]) #anthromes
sitesHydeL$V1 <- NULL #anthromes
colnames(sitesHydeL) <- c("site", "time", "anthromes")

sitesHydeL[,4] <- as.factor(sitesHydeL[,4]) #cropland
sitesHydeL[,5] <- as.factor(sitesHydeL[,5])
sitesHydeL$name <- NULL
colnames(sitesHydeL)[c(1,3,4)] <- c("site", "scenario", "time")

save(sitesHydeL, file="results/sitesHydesvaluesLong.Rdata")
save(sitesHydeL, file="results/sitesAnthromesLong.Rdata")
save(sitesHydeL, file="results/sitesCroplandLong.Rdata")

#summarise
load("results/sitesHydesHumanPopLong.Rdata")
load("results/sitesAnthromesLong.Rdata")
load("results/sitesCroplandLong.Rdata")

      humanPop <- sitesHydeL %>% group_by(site, scenario) %>% 
        summarise(pop.ave=mean(pop,na.rm=T)) %>% 
        spread(scenario, pop.ave)
      
      anthromes <- sitesHydeL %>% group_by(site, time) %>% 
        summarise(anthromes_sum=sum(anthromes,na.rm=T)) %>% 
        spread(site, pop.ave)
      
      cropland <- sitesHydeL %>% group_by(site, time) %>% 
        summarise(crop_avg=mean(cropland,na.rm=T)) %>% 
        spread(site, crop_avg)
      
      
      predictors$humanDensity  <- humanPop$pop.ave
      predictors$humanPresence <- "human_presence"
      predictors$humanPresence[predictors$humanDensity==0] <- "human_absence"
      
      save(predictors, file="results/predictors.Rdata")
      
## spatial joining merging lake sites and mountain andes dataset
## https://stat.ethz.ch/pipermail/r-help/2008-September/173773.html 

source("code/functions.R")
      
#demo
track <- data.frame(xt=runif(10,0,360), yt=rnorm(10,-90, 90))
classif <- data.frame(xc=runif(10,0,360), yc=rnorm(10,-90, 90),
v1=letters[1:20], v2=1:20)

dist.merge(track, classif, 'xt', 'yt', 'xc', 'yc')

#own data
mountain_andes_subset2 <- mountain_andes_subset %>%
  transform(LongitudeD=Longitude, LatitudeD=Latitude)

sitesAndes <- dist.merge(coorSites, mountain_andes_subset2, 
                   'Longitude', 'Latitude', 'LongitudeD', 'LatitudeD')

sitesAndes <- sitesAndes[,-ncol(sitesAndes)] #remove last column

predictors$sites <- rownames(predictors)

all_predictors <- plyr::join(sitesAndes,predictors, by="sites")
all_predictors <- all_predictors[,-c(4,5,6,17,18)]

write.csv(all_predictors, file="results/lakes_predictors.csv")

## map predictors

human<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=humanDensity),size=4, shape=21)+
  scale_fill_viridis(trans="sqrt")+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

croplands<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=croplandhistoric),size=4, shape=21)+
  scale_fill_viridis(trans="sqrt")+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

erosion<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=Avg_HERMER),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

relief<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=RELIEF_MN),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

rugged<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=Rugg_avg),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

soil<- 
  southamerica + 
  geom_point(data=all_predictors, aes(x=Longitude, y=Latitude, fill=soil_var),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position=c(0.85,0.2))

## plot histogram data distributions
hist_humanDensity<-
  ggplot(all_predictors, aes(x=humanDensity, fill=..x..))+
  geom_histogram()+
  theme_classic()+
  #scale_fill_distiller(palette = "YlGnBu")+
  scale_fill_viridis()+
  theme(legend.position="none")+
  scale_x_continuous(trans="log1p")+
  theme_classic()

hist_croplands<-
  ggplot(all_predictors, aes(x=croplandhistoric, fill=..x..))+
  geom_histogram()+
  theme_classic()+
  #scale_fill_distiller(palette = "YlGnBu")+
  scale_fill_viridis()+
  theme(legend.position="none")+
  scale_x_continuous(trans="log1p")+
  theme_classic()

hist_erosion<-
  ggplot(all_predictors, aes(x=Avg_HERMER, fill=..x..))+
  geom_histogram()+
  theme_classic()+
  scale_fill_viridis()+
  theme(legend.position="none")+
  theme_classic()

hist_RELIEF<-
  ggplot(all_predictors, aes(x=RELIEF_MN, fill=..x..))+
  geom_histogram()+
  theme_classic()+
  scale_fill_viridis()+
  theme(legend.position="none")+
  theme_classic()

hist_rugged<-
  ggplot(all_predictors, aes(x=Rugg_avg, fill=..x..))+
  geom_histogram()+
  theme_classic()+
  scale_fill_viridis()+
  theme(legend.position="none")+
  theme_classic()
