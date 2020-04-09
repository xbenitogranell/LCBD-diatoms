
library(vegan)
library(cluster)
library(tidyverse)

## Read in data
lake_predictors <- read.csv("results/lakes_predictors.csv", row.names=1) %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region,"Cusco","Peruvian Andes-Cusco")) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Chile-Andean", "Chile-Andes"))

rownames(lake_predictors) <- lake_predictors$Row.names

chile_coastal <- paste(c("Chl-Crr_Aculeo","Chl-Crr_Batuco","Chl-Crr_CingsdNm", "Chl-Crr_LgnAlbfr"), collapse = '|')
diatoms <- read.csv("results/diatoms.csv") %>%
  filter(!str_detect(X, chile_coastal))

rownames(diatoms) <- diatoms[,1] #remove first column
diatoms[,1] <- NULL

diatoms <- diatoms[,colSums(diatoms)>0]


#This is for prepare data to map diatom community cluster composition
##this is to group species by ecological groups
###Read data diatom training set and environmental data
environmental_data_lakes <- read.csv("data/environmental_data_lakes.csv") %>%
  mutate(lake_depth_ratio=Lake_area/Depth_avg) %>%
  mutate(lake_catch_ratio=Lake_area/Wshd_area) %>%
  mutate(catch_vol_ratio=Wshd_area/Vol_total)

rownames(environmental_data_lakes) <- environmental_data_lakes$code

training <- read.csv("data/diatomsTrainingSetDEF.csv", row.names = 1) #with updated diatom taxonomy and selected spp (>3% of RA and present in >2 samples) plus Miriam's Llaviucu slides
richness<- apply(training>0,1,sum)

#Regions
lake_regions <- read.csv("data/regions.csv", row.names = 1)

##Merge training set and regions datasets
modern_lakes <- merge(training, lake_regions, by="row.names")

#transform dataframe to tidy format
df_thin <- modern_lakes %>%
  gather(key = taxa, value = count, -Row.names, -region)#don't gather region

#import dataframe wiht old and new names to group
changes_training <- read.csv("data/old_new_nms_trainingset.csv", stringsAsFactors = FALSE)

#spread
new <- df_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes_training$old, to = changes_training$new_1)) %>%
  group_by(region, Row.names, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)

levels(new$region)

# comment spp when analyzing ecological groups data
new <- new %>%
  filter(Fragilaria.crotonensis < 40) %>%
  filter(Cyclostephanos.tholiformis < 40) %>%
  filter(Cyclostephanos.andinus < 40) %>%
  as.data.frame()

training <- new[, -which(names(new) %in% c("Row.names", "region"))]

#Merge diatom training set and environmental data of lakes
row.names(training) <- new[, which(names(new) %in% c("Row.names"))]
env_surf <- merge(training,environmental_data_lakes, by="row.names")

#For extracting spp from trainingset
training2 <- env_surf[,2:235]

rowSums(training2)
#write.csv(training2, "results/diatoms2.csv")

#For extracting environmental variables from diatom training set
env_data_lakes <- env_surf[,236:ncol(env_surf)]

row.names(env_data_lakes) <- env_data_lakes$code

# For merging environmental dataset with lake regions
env_data_lakes <- merge(env_data_lakes,lake_regions, by="row.names")
row.names(env_data_lakes) <- env_data_lakes$code
env_data_lakes <- env_data_lakes[,-1]

#Save dataframe with Titicaca region
env_data_lakes <- merge(lake_regions, env_data_lakes, by="row.names")
write.csv(env_data_lakes, "results/lake_predictors2.csv")
  
new <- df_thin %>% 
  mutate(taxa = plyr::mapvalues(taxa, from = changes_training$old, to = changes_training$new_2)) %>% #ecological grouping
  group_by(region, Row.names, taxa) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  ungroup() %>%
  group_by(Row.names, region) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(plank=sum(count[taxa=="freshwater_planktic" | taxa=="tycoplanktonic"])) %>%
  mutate(benthic=sum(count[taxa=="epiphytics"| taxa== "saline" | taxa=="benthic"])) %>%
  mutate(P_B=plank/benthic) %>%
  mutate(P_B2=(plank-benthic)/(plank+benthic)) %>% #[-1(benthic dominated) to 1(planktic dominated)]
  ungroup() 

#make it wide
lake_diatom_ratios <- new %>%
  dplyr::select(region, Row.names, taxa, P_B2, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) 

lake_diatom_ratios[is.na(lake_diatom_ratios)] <- 0

new <- lake_diatom_ratios %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region,"Cusco","Peruvian Andes-Cusco")) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes-Wet Puna"))%>%
  mutate(region=str_replace(region, "Chile-Andean", "Chile-Andes")) %>%
  as.data.frame()

cluster_data <- new[,colnames(new[,3:9])]
region <- factor(new$region, levels = c("Colombia-Andes", "Ecuador-Andean",
                                        "Ecuador-Interandean", "JuninPlain", "Peruvian Andes-Cusco", "Peruvian Andes-Wet Puna", "Titicaca-lake",
                                        "Desaguadero-SAltiplano", "Sud Lipez", "Chile-Andes"))


top_spp <- colnames(cluster_data)
top_spp <- match(top_spp, names(cluster_data))
region_n <- as.numeric(region)

#### Building aggregate composition data for each region ####
## Averaging over each cluster, and setting each row to sum to one,then taking  column means
## code from Pedersen et al 2017
library(plyr)
cluster_composition = ddply(cluster_data, .(region),
                            function(x){
                              x = decostand(x,method ="total") 
                              out_data = data.frame(species = c(names(x)[top_spp]),
                                                    proportion = rep(0,times=length(top_spp)))
                              for(i in 1:length(top_spp)){
                                out_data[i,2] = mean(x[,top_spp[i]])
                              }
                              return(out_data)
                            })

cluster_composition <- cluster_composition %>% filter(!species=="P_B2") %>%
  mutate(species=str_replace(species, "benthic", "Benthic")) %>%
  mutate(species=str_replace(species, "ephiphytic", "Ephiphytic")) %>%
  mutate(species=str_replace(species, "freshwater_planktic", "Freshwater planktic")) %>%
  mutate(species=str_replace(species, "oligosaline_planktic", "Oligosaline planktic")) %>%
  mutate(species=str_replace(species, "saline", "Saline")) %>%
  mutate(species=str_replace(species, "tycoplanktonic", "Tychoplanktonic"))
  
cluster_composition$species_f <- factor(cluster_composition$species,
                      levels = c("Benthic","Ephiphytic", "Saline", "OligoSaline planktic",
                                 "Tychoplanktonic", "Freshwater planktic"))

cluster_order <- plyr::ddply(cluster_composition,.(region_n), function(x){
  proportion <- sum(x$proportion)
  return(data.frame(proportion=proportion))
})


# cluster_composition$region <- levels(region)
# cluster_order <- cluster_order$region_n[order(cluster_order$proportion,
#                                                 decreasing = T)]
# cluster_composition$region <- factor(cluster_composition$region_n,
#                                         levels = cluster_order)

library(viridis)
region_palette <- viridis(length(levels(cluster_composition$region))) #levels of sigClust

# Plot
comp_plot <- ggplot(aes(x=region, y=proportion),
                    data=cluster_composition)+
  geom_bar(stat="identity", aes(fill=species_f,order=region))+
  labs(fill="Groups")+
  scale_fill_viridis(discrete = TRUE, option = "D")+
  # annotate(x=factor(1:10),y = rep(-0.1,times=10),
  #          colour=region_palette,
  #          geom="point",size=10,shape=18)+
  #coord_cartesian(ylim=c(-0.2,1))+
  scale_y_continuous("Average % of community",breaks = c(0,0.25,0.5,0.75,1))+
  scale_x_discrete("Region")+
  theme_bw()+
  theme(text= element_text(size=12),panel.grid.major.x = element_blank(),
        axis.text.x=element_text(size=10, angle = 45, hjust = 1),
        axis.ticks.x=element_blank(),legend.direction="horizontal",
        legend.position = c(0.5,0.75),
        legend.text=element_text(size=8),
        legend.title = element_blank())
comp_plot

#plot modern lake database with lake regions
library(maps)
library(rwordlmap)
world <- map_data("world")

interest <- c("Colombia", "Ecuador", "Peru", "Bolivia", "Chile", "Argentina")
countries <- world %>% filter(str_detect(region, interest))

lonlatLakes <- read.csv("results/lake_predictors2.csv", row.names=1) %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region, "Cusco","Peruvian Andes-Cusco")) %>%
  mutate(region=str_replace(region, "Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Chile-Andean", "Chile-Andes"))%>%
  mutate(region_f=factor(region)) %>%
  dplyr::select(long, lat, region_f)

#region_palette <- viridis(length(levels(lonlatLakes$region_f))) #levels of sigClust

#Plot elevation base map (raster library)
library(raster)
DEM <- raster("dem/dem2.bil")
ext<-extent(-82,-35,-40,15)
altmod<-crop(DEM,ext)
plot(altmod)

#convert the raster to points for plotting
map.p <- rasterToPoints(altmod)

#Make the points a dataframe for ggplot
df <- data.frame(map.p)

#Make appropriate column headings
colnames(df) <- c("long", "lat", "Elevation")

## DEM South America
plot_sa <- ggplot(data=df, aes(y=lat, x=long)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="", color="black")+
  geom_raster(aes(fill=Elevation)) +
  scale_fill_distiller(palette = "RdYlBu") +
  #scale_colour_gradient(high = "red") +
  theme(legend.position = "right")+
  #coord_equal() +
  coord_map("albers", parameters = c(-100, -100),  ylim=c(-40,15), xlim=c(-82,-40)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()

# South America basemap
southamerica <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="lightgrey", color="black")+
  theme(legend.position = "right")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(-40,15), xlim=c(-82,-40)) +
  #coord_map("albers", parameters = c(-100, -100),  ylim=c(-10,15), xlim=c(-82,-60)) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()

lonlatLakes$region_f <- factor(lonlatLakes$region, levels = c("Colombia-Andes", "Ecuador-Andean",
                                        "Ecuador-Interandean", "JuninPlain", "Peruvian Andes-Cusco", "Peruvian Andes-Wet Puna", "Titicaca-lake",
                                        "Desaguadero-SAltiplano", "Sud Lipez", "Chile-Andes"))

cluster_map_plt <- southamerica + 
  geom_point(data=lonlatLakes, aes(x=long, y=lat, col=region_f),
                                             shape=18, size=4)+
  #this comes from 4-code-timeseries
  geom_point(data=lake_cores2, aes(x=Longitude, y=Latitude, shape=lake_f), size=3) +
  theme(text= element_text(size=12),
        legend.direction="vertical",
        legend.position = c(0.75,0.5),
        legend.box.just = c("top"), 
        legend.title = element_blank(),
        legend.background =element_rect(fill=alpha('white', 0.5)))
cluster_map_plt


library(cowplot)
plt <- plot_grid(cluster_map_plt,comp_plot, align="hv", axis="t", ncol = 2,
                 rel_widths = c(1,1),labels = "auto")
plt
ggsave("figures/study_region_groups_NEW.eps", plt, height = 8, width = 10)

#plot modern lake database 
#read abiotic variables (Antonelli et al 2018 Nature Geosciences)
mountain_andes <- read.csv("data/mountain_andes.csv", row.names = 1)

mountain_andes_coord <- mountain_andes %>% 
  dplyr::select(Centroid_X, Centroid_Y)

lakes_map_plt <- southamerica + geom_point(data=lonlatLakes, aes(x=long, y=lat, col=region_f),
                                           shape=18, size=4)+
  scale_color_manual(values = region_palette) + 
  geom_point(data=mountain_andes_coord, aes(x=Centroid_X, y=Centroid_Y))+
  theme(text= element_text(size=15),
        legend.direction="vertical",
        legend.position = "right")
lakes_map_plt

ggsave("figures/plot_lakes_regions.png", lakes_map_plt, height = 8, width = 10)



