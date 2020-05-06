library(tidyverse)

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
