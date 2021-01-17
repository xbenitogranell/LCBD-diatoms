library(vegan)
library(cluster)
library(tidyverse)

## Read in data
# lakes_predictors from
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

# select our outlier samples with high planktic spp abundance
new <- new %>%
  filter(Fragilaria.crotonensis < 40) %>%
  filter(Cyclostephanos.tholiformis < 40) %>%
  filter(Cyclostephanos.andinus < 40) %>%
  as.data.frame()

training <- new[, -which(names(new) %in% c("Row.names", "region"))]

cluster_data <- new[,colnames(new[,3:ncol(new)])]

new <- new %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region,"Cusco","Peruvian Andes-Cusco")) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes-Wet Puna")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes-Wet Puna"))%>%
  mutate(region=str_replace(region, "Chile-Andean", "Chile-Andes")) %>%
  as.data.frame()

cluster_data <- new[,colnames(new[,3:ncol(new)])]
region <- factor(new$region, levels = c("Colombia-Andes", "Ecuador-Andean",
                                        "Ecuador-Interandean", "JuninPlain", "Peruvian Andes-Cusco", "Peruvian Andes-Wet Puna", "Titicaca-lake",
                                        "Desaguadero-SAltiplano", "Sud Lipez", "Chile-Andes"))


top_spp <- colnames(cluster_data)
top_spp <- match(top_spp, names(cluster_data))
region_n <- as.numeric(region)

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

# here select most abundant spp per region
most_abundant <- cluster_composition %>% 
  group_by(region, species) %>%
  filter(proportion>0.07) %>% 
  ungroup() 

most_abundant$region <- factor(most_abundant$region, levels = c("Colombia-Andes", "Ecuador-Andean",
                                        "Ecuador-Interandean", "JuninPlain", "Peruvian Andes-Cusco", "Peruvian Andes-Wet Puna", "Titicaca-lake",
                                        "Desaguadero-SAltiplano", "Sud Lipez", "Chile-Andes"))

plt <- ggplot(data=most_abundant, aes(x=region, y=proportion)) +
  geom_bar(stat="identity", aes(order=region,fill=species)) +
  scale_y_continuous("Relative abundance (%)",breaks = c(0,0.25,0.5,0.75,1))+
  scale_x_discrete("Region")+
  theme_bw()+
  theme(text= element_text(size=12),panel.grid.major.x = element_blank(),
        axis.text.x=element_text(size=10, angle = 45, hjust = 1),
        axis.ticks.x=element_blank(),legend.direction="horizontal",
        legend.position = "top",
        legend.text=element_text(size=8),
        legend.title = element_blank())
plt  

ggsave("figures/Fig_abundantspp.png", plt, height = 8, width = 12)


