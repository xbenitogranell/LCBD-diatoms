#clear workspace
rm(list=ls(all=TRUE))
dev.off()
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

#load functions used
library(adespatial)
library(tidyverse)
library(raster)
library(ggplot2)
library(viridis)
library(vegan)
library(cluster)


## HYDE DB
#summarise
load("results/sitesHydesHumanPopLong.Rdata")
humanPop <- sitesHydeL %>% group_by(site, time) %>% 
  summarise(pop_avg=mean(pop,na.rm=T)) %>% 
  spread(time, pop_avg)

load("results/sitesCroplandLong.Rdata")
cropland <- sitesHydeL %>% group_by(site, scenario) %>% 
  summarise(crop_avg=mean(cropland,na.rm=T)) %>% 
  spread(scenario, crop_avg)

load("results/sitesAnthromesLong.Rdata")


predictors <- env_data_lakes

predictors$humanDensity  <- humanPop$baseline
predictors$croplandhistoric <- cropland$baseline

save(predictors, file="results/predictors.Rdata")

#######
## LCBD and SCBD
lake_predictors <- read.csv("results/lakes_predictors.csv", row.names=1) %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes"))

chile_coastal <- paste(c("Chl-Crr_Aculeo","Chl-Crr_Batuco","Chl-Crr_CingsdNm", "Chl-Crr_LgnAlbfr"), collapse = '|')
diatoms <- read.csv("results/diatoms.csv") %>%
  filter(!str_detect(X, chile_coastal))

rownames(diatoms) <- diatoms[,1] #remove first column
diatoms[,1] <- NULL
diatoms <- diatoms[,colSums(diatoms)>0]

#remove diatom names from the dataframe
#lake_predictors2 <- lake_predictors[, -which(names(lake_predictors) %in% names(diatoms))]
#lake_predictors3 <- lake_predictors2[,-124:ncol(lake_predictors2)]

?beta.div
diatomsLCBD <- beta.div(diatoms, method = "hellinger", nperm = 999, adj = TRUE, sqrt.D=FALSE)

## plot a map of the LCBD indices
LCBD <- data.frame(diatomsLCBD$LCBD)
LCBD$Longitude <- lake_predictors$Longitude
LCBD$Latitude <- lake_predictors$Latitude
colnames(LCBD)[1] <- c("LCBD")
LCBD$sign <- data.frame(diatomsLCBD$p.LCBD)
colnames(LCBD[,4]) <- "sign"
LCBD$sign.ad <- data.frame(diatomsLCBD$p.adj)

#richness calculation
richness<- apply(diatoms>0,1,sum)

#LCBD - richness relationship
LCBD$richness <- richness
plot(LCBD$richness, LCBD$LCBD)
cor.test(LCBD$richness, LCBD$LCBD, method = "spearman")

LCBD_df <- data.frame(LCBD)
LCBDrich <- ggplot(data=LCBD_df, aes(x=richness, y=LCBD))+
  geom_point()+
  geom_smooth(method=lm)+
  annotate("text", x = 35, y = 0.0065, label = "Spearman rho = -0.58")+
  annotate("text", x = 35, y = 0.0064, label = "italic(p) < 0.01", parse=TRUE)+
  xlab("Species richness")+
  theme_bw()

ggsave("figures/LCBD_rich.png", LCBDrich, height = 8, width = 10)


##
LCBD_map <- southamerica + 
  geom_point(data=LCBD, aes(x=Longitude, y=Latitude, colour=sign.ad<0.05,size=LCBD$LBCD))+
  theme(text=element_text(size=15),
        #legend.direction="vertical",
        legend.position = "top")+
  theme_classic() +
  scale_size(range = c(1,4), breaks = c(0.004,0.005, 0.0065), "LCBD")
LCBD_map

ggsave("figures/LCBD_map.png", LCBD_map, height = 8, width = 10)

## Decompose beta diversity components
?beta.div.comp
A <- beta.div.comp(decostand(diatoms, method="hellinger"), coef = "S", quant = TRUE)

Replmatrix <- as.matrix(A$repl)
repl <- Replmatrix[row(Replmatrix) == col(Replmatrix) + 1]

Richmatrix <- as.matrix(A$rich)
rich <- Richmatrix[row(Richmatrix) == col(Richmatrix) + 1]


LCBDrepl <- LCBD.comp(Replmatrix, sqrt.D = FALSE)
LCBDrich <- LCBD.comp(Richmatrix, sqrt.D = FALSE)
  plot(LCBD$richness, LCBDrich$LCBD)
  cor.test(LCBD$richness, LCBDrich$LCBD, method="spearman")

div_predictors <- cbind(LCBD$LBCD, LCBDrepl$LCBD, LCBDrich$LCBD, lake_predictors)
div_predictors <- div_predictors[,-which(names(div_predictors) %in% c("Row.names", "sites"))] #remove extraneous cols

colnames(div_predictors)[1:3] <- c("LCBD", "LCBDrepl", "LCBDrich")
write.csv(div_predictors, "results/div_predictors2.csv")


##
# richLCBD <- list()
# 
# richLCBDts <- function(i, cores,...) {
#   core <- cores[[i]]
#   core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
#   core[is.na(core)] <- 0
#   core <- core[, colSums(core) > 0] #select only present species
#   richLCBD$richness[[i]] <- apply(core>0,1,sum)
#   #LCBDts$richness2[[i]] <- specnumber(core)
#   #raremax <- min(rowSums(core))
#   #Srare <- rarefy(core, raremax)
#   core_prop <- tran(core, method = "percentage")
#   LCBD <- beta.div(core_prop, method = "hellinger", nperm = 999, adj = TRUE)
#   richLCBD$LCBD[[i]] <- LCBD[["LCBD"]][-1]
#   richLCBD$SCBD[[i]] <- LCBD[["SCBD"]]
#   
#   
# }

### LCBD per region (as in)
## split diatom training set by regions and reassemble



    # diatom_regions <- merge(diatoms, lake_regions, by="row.names")
    # rownames(diatom_regions) <- diatom_regions$Row.names
    # diatom_regions$Row.names <- NULL
    
    # diatomRegions[["Argentina-MarChiquita"]]<-NULL
    # diatomRegions[["Colombia-Lowlands"]]<-NULL
    # diatomRegions[["Argentina-Andean"]] <- NULL
    # diatomRegions[["Bolivia-Amazon Lowlands"]] <- NULL
    # diatomRegions[["Bolivia-Beni"]]<- NULL
    # diatomRegions[["Cochabamba"]]<- NULL
    # diatomRegions[["Colombia-lowlands"]]<- NULL
    # diatomRegions[["Ecuador-Amazonia"]]<- NULL
    # diatomRegions[["Ecuador-Coastal"]]<- NULL
    # diatomRegions[["Lauca Basin"]]<- NULL
    # diatomRegions[["Pantanal"]]<- NULL
    # diatomRegions[["Lauca Basin"]]<- NULL
    # diatomRegions[["Peru-Amazon Lowlands"]] <- NULL
    # diatomRegions[["Sehuencas-WesternAndes"]] <- NULL
    # diatomRegions[["Sorata-WesternAndes"]] <- NULL
    # diatomRegions[["Tierra del Fuego"]] <- NULL
    # diatomRegions[["Titicaca-lake"]] <- NULL
    # diatomRegions[["Titicaca-watershed"]] <- NULL


diatomRegions <- split(diatoms, lake_predictors$region)

# function to calculate LCBD for each Andean region
LCBDts <- function(i, cores, ...) {
  core <- cores[[i]]
  #core <- core[ , -which(names(core) %in% c("X"))] # drop year & depths vars
  core[is.na(core)] <- 0
  core <- core[, colSums(core) > 0] #select only present species
  #core <- analogue::tran(core, method="hellinger") #Hellinger transform relative abundance data
  diss <- analogue::distance(core, method = "bray")
  LCBD <- beta.div(core, method = "hellinger", nperm = 999, adj = TRUE)
  LCBD <- LCBD[["LCBD"]][-1]

  #SCBD <- LCBD[["SCBD"]]
  
  A <- beta.div.comp(core, coef = "S", quant = TRUE)
  repl <- as.matrix(A$repl)
  LCBDrepl <- LCBD.comp(repl, sqrt.D = TRUE)
  LCBDrepl <- LCBDrepl[["LCBD"]][-1]

  repl <- repl[row(repl) == col(repl) + 1]
  
  rich <- as.matrix(A$rich)
  LCBDrich <- LCBD.comp(rich, sqrt.D = TRUE)
  LCBDrich <- LCBDrich[["LCBD"]][-1]

  rich <- rich[row(rich) == col(rich) + 1]
  
  # LCBDrich <- LCBDrich[row(LCBDrich) == col(LCBDrich) + 1]
  
  cbind.data.frame(LCBD,LCBDrepl,LCBDrich,repl,rich) #combine extracted columns and remove first row to match with scd
  
}

## apply  function to each core
regionLCBD <- lapply(seq_along(diatomRegions), LCBDts, cores=diatomRegions)
names(regionLCBD) <- names(diatomRegions)

#extract dataframes from list
data_plt <- plyr::ldply(regionLCBD, data.frame) %>%
  mutate(region=factor(.id))

cols <- c("LCBDrepl"="#f04546","LCBDrich"="#3591d1")

## plot the data
regionLCBD_plt <- ggplot(data_plt, aes(x=region, y=repl))+
  geom_boxplot()+
  theme_classic()+
  coord_flip()
regionLCBD_plt
