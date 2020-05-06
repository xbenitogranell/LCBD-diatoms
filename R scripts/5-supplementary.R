library(tidyverse)
library(analogue)
library(ggplot2)

sitesdb<- read.csv("data/biogeographySites.csv", row.names=1) %>%
  select(SiteName, code, Lat.DD.S, Long.DD.W, Year)

#env_data_lakes come from R script 1-code-community composition.R
sites <- merge(sitesdb, env_data_lakes, by="code") %>%
  select(SiteName, code, long, lat, region.x, Year) %>%
  rename(Lake=SiteName) %>%
  rename(Longitude=long) %>%
  rename(Latitude=lat) %>%
  rename(Region=region.x) %>%
  filter(!str_detect(Region, "Chile-Coastal")) %>%
  mutate(Region=str_replace(Region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(Region=str_replace(Region, "Peru-Andean eastern", "Peruvian Andes")) %>%
  mutate(Region=str_replace(Region, "Peru-Andean Puna", "Peruvian Andes"))

names(sites)
unique(sites$Region)
write.csv(sites, "results/sitesdb.csv")

kable(sites) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


#### Stratigraphic plots
mergedCores <- read.csv("mergedCores_counts4.csv")[,-1] 
diatoms_save <- mergedCores #save dataframe

#Gather
spp_thin <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake)#don't gather depths, ages and lake variables


#import dataframe wiht old and new names to group
changes <- read.csv("old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names

#spread--> wide format
spp_wide <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)

#no spread --> long format
spp_long <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count))


#filter cores
select <- c("llaviucu")
select <- c("pinan")
select <- c("yahuarcocha")
select <- c("fondococha")
select <- c("titicaca")
select <- c("umayo")


core_lake <- spp_long %>%
  filter(str_detect(lake, select)) %>% #select lake
  filter(!upper_age==0.0) %>%
  group_by(taxa) %>%
  filter(count > 0) %>% #remove species with 0 counts
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>% #calculate RA
  ungroup()

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- core_lake %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 3) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- core_lake %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide
core_counts_wide <- core_counts_common %>%
  dplyr::select(depth, lake, upper_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent)

length(core_common_taxa)


## plot using analogue Stratiplot (need data input in wide format --> spreaded)
png("Stratiplot_lakename.png", width = 11, height = 8, res = 300, units = "in")

Stratiplot(
  core_counts_wide %>% select(-depth, -lake, -upper_age),
  core_counts_wide$upper_age,
  ylab = "Cal yr BP", 
  xlab = "Relative abundance (%)",
  # adds padding to the top of the plot
  # to fix cut-off taxa names
  topPad = 10, 
  # make the plot type a "bar" plot
  type = "h", 
  #sort = "wa",
  # make the bar colour black
  col = "black")

dev.off()  