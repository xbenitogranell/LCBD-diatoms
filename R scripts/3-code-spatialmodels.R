#load libraries for functions used
source("R scripts/functions.R")

library(mgcv)
library(visreg)
library(usdm)
library(dismo)
library(gbm)
library(tidyverse)
library(ggRandomForests)
library(randomForestSRC)

## Read in data
lake_predictors <- read.csv("results/lakes_predictors.csv", row.names=1) %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes"))

chile_coastal <- paste(c("Chl-Crr_Aculeo","Chl-Crr_Batuco","Chl-Crr_CingsdNm", "Chl-Crr_LgnAlbfr"), collapse = '|')
diatoms <- read.csv("results/diatoms.csv") %>%
  filter(!str_detect(X, chile_coastal))

#Subset variables of interest
names(lake_predictors)
data_df <- lake_predictors %>% dplyr::select(Longitude, Latitude, 
                                          Rugg_avg, Avg_HERMER, RELIEF_MN, soil_var, humanDensity, croplandhistoric,
                                          Water.T, Cond, pH, TP, Cl, Na,SO4, Elevation, area_waterbody,X..aquatic.habitat,
                                          MAT, P.season, MAP, T.season)


# Plot
pairs(data_df, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

#transform variables to meet assumptions of homogenity of variances
data_transf <- transform(data_df, Water.T=log10(Water.T+0.25), Cond=log10(Cond+0.25), TP=log10(TP+0.25), Cl=log10(Cl+0.25), Na=log10(Na+0.25),
                         SO4=log10(SO4+0.25),MAP=log10(MAP+0.25), MAT=log10(MAT+0.25), P.season=log10(P.season+0.25), T.season=log10(T.season+0.25),
                         Rugg_avg=log10(Rugg_avg+0.25), Avg_HERMER=log10(Avg_HERMER+0.25), RELIEF_MN=log10(RELIEF_MN+0.25), 
                         humanDensity=log10(humanDensity+0.25), croplandhistoric=log10(croplandhistoric+0.25),
                         elevation=sqrt(Elevation), area_waterbody=log10(area_waterbody+0.25), aquatic.habitat=log10(X..aquatic.habitat+0.25))

pairs(data_transf, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

library(usdm)
vifstep(data_transf, th=10) # threshold set to VIF below 15

#transform data to improve linearity
data_transf_red <- data_transf %>% dplyr::select(Longitude, soil_var, humanDensity, croplandhistoric, 
                                                 Water.T, Cond, TP, Cl, Na, pH, area_waterbody, MAP,P.season,T.season, 
                                                 aquatic.habitat)

vifstep(data_transf_red, th=10)



### GAM Models
LCBD_data <- read.csv("results/div_predictors.csv") %>%
  filter(!str_detect(region, "Chile-Coastal")) %>%
  dplyr::select(LCBD, LCBDrepl, LCBDrich) %>%
  cbind(data_transf_red) %>%
  mutate(region=lake_predictors$region) %>%
  mutate(Rugg_avg=lake_predictors$Rugg_avg) %>%
  mutate(Latitude=lake_predictors$Latitude) %>%
  mutate(region=str_replace(region,"Cochabamba","Desaguadero-SAltiplano")) %>%
  mutate(region=str_replace(region, "Peru-Andean eastern", "Peruvian Andes")) %>%
  mutate(region=str_replace(region, "Peru-Andean Puna", "Peruvian Andes")) %>%
  mutate(region_f=as.factor(region))
  
names(LCBD_data)

# relationships between LCBDrepl and LBCDrich
LCBD_data$region_f <- factor(LCBD_data$region, levels = c("Colombia-Andes", "Ecuador-Andean",
                                        "Ecuador-Interandean", "JuninPlain", "Cusco", "Peruvian Andes", "Titicaca-lake",
                                        "Desaguadero-SAltiplano", "Sud Lipez", "Chile-Andean"))


LCBD_data_red <- LCBD_data %>% filter(LCBDrepl > 0.005)

LCBDcorr <- ggplot(data = LCBD_data) +  
  geom_point(aes(x=LCBDrepl, y = LCBDrich, color = factor(region_f))) +
  geom_smooth(aes(x=LCBDrepl, y = LCBDrich))+
  labs(color = "Lake regions")+
  theme_bw()

ggsave("figures/LCBDrelationship.png", plot=LCBDcorr, height=8, width=10,units="in",
       dpi = 400)

cor.test(LCBD_data$LCBDrepl, LCBD_data$LCBDrich, method = "spearman")

# relationships betweeb LCBD and environmental predictors
LCBD_env <- LCBD_data %>% select(LCBD, LCBDrepl, LCBDrich, humanDensity, pH, Cond, MAP, Rugg_avg, T.season) %>%
  mutate(Rugg_avg_log=log10(Rugg_avg+0.25))%>%
  mutate(Rugg_avg_sqrt=sqrt(Rugg_avg))

par(mfrow=c(2,2))
plot(LCBD_env$MAP, LCBD_env$LCBD, xlab = "log10 MAP", ylab="LCBD")
plot(LCBD_env$T.season, LCBD_env$LCBD, xlab = "log10 Tseason", ylab="LCBD")

cor <- cor.test(LCBD_env$Rugg_avg, LCBD_env$LCBDrepl, method = "spearman")
plot(LCBD_env$Rugg_avg, LCBD_env$LCBDrepl, xlab = "Ruggedness", ylab="LCBDrepl")
title(paste("Spearman rho=",round(cor$estimate, digits=2), "p=", round(cor$p.value, digits=3)))

cor <- cor.test(LCBD_env$Rugg_avg, LCBD_env$LCBDrich, method = "spearman")
plot(LCBD_env$Rugg_avg, LCBD_env$LCBDrich, xlab = "Ruggedness", ylab="LCBDrich")
title(paste("Spearman rho=",round(cor$estimate, digits=2), "p=", round(cor$p.value, digits=3)))

cor <- cor.test(LCBD_env$T.season, LCBD_env$LCBDrepl, method = "spearman")
plot(LCBD_env$T.season, LCBD_env$LCBDrepl, xlab = "log10 Tseason", ylab="LCBDrepl")
title(paste("Spearman rho=",round(cor$estimate, digits=2), "p=", round(cor$p.value, digits=5)))

cor <- cor.test(LCBD_env$T.season, LCBD_env$LCBDrich, method = "spearman")
plot(LCBD_env$T.season, LCBD_env$LCBDrich, xlab = "log10 Tseason", ylab="LCBDrich")
title(paste("Spearman rho=",round(cor$estimate, digits=2), "p=", round(cor$p.value, digits=3)))


cor.test(LCBD_env$LCBDrepl, LCBD_env$T.season, method = "spearman")
cor.test(LCBD_env$LCBDrich, LCBD_env$T.season, method = "spearman")




# Linear model + spatial smooths + random effects for region
set.seed(10) #set a seed so this is repeatable

mod1 <- gam(LCBD ~ s(Latitude, Longitude) + s(region_f, bs="re") +
                scale(MAP) + scale(T.season) + scale(P.season) +
                scale(soil_var) + scale(area_waterbody) + scale(Rugg_avg)+    
                scale(humanDensity) + scale(croplandhistoric) + scale(Water.T) + scale(pH) + scale(Cond)
                + scale(TP) + scale(Cl) + scale(Na), method="REML", data=LCBD_data, select=TRUE, na.action = 'na.omit', family=gaussian)
summary(mod1)
gam.check(mod1)

visreg(mod1,'humanDensity', gg=T)+
  ylab("n partial residuals")+
  theme_classic()

# Linear model + spatial smooths + random effects for region
mod2 <- gam(LCBDrepl ~ s(Latitude, Longitude) + s(region_f, bs="re") +
              scale(MAP) + scale(T.season) + scale(P.season) +
              scale(Rugg_avg) + scale(soil_var) + scale(area_waterbody) +    
              scale(humanDensity) + scale(croplandhistoric) + scale(Water.T) + scale(pH) + scale(Cond) + scale(TP) + 
              scale(Cl) + scale(Na),method="REML", data = LCBD_data, select=TRUE, na.action = 'na.omit')

summary(mod2)
visreg(mod2,'area_waterbody', gg=T)+
  ylab("n partial residuals")+
  theme_classic()

# Linear model + spatial smooths + random effects for region
mod3 <- gam(LCBDrich ~ s(Latitude, Longitude) + s(region_f, bs="re") +
              scale(MAP) + scale(T.season) + scale(P.season) +
              scale(Rugg_avg) + scale(soil_var) + scale(area_waterbody)+  
              scale(humanDensity) + scale(croplandhistoric) + scale(Water.T) + scale(pH) + scale(Cond) + scale(TP) + 
              scale(Cl) + scale(Na),method="REML", data = LCBD_data, select=TRUE, na.action = 'na.omit')

summary(mod3)

#Check spatial autocorrelation of model residuals
# my.coord <- coordinates(coord)
# 
# #no duplicate rows
# my.coord.nd <- my.coord[!duplicated(my.coord), ]
# 
# duplicated(my.coord)
# sel <- c(113)
# 
# 
# model1.res2 <- as.matrix(model1.res)
# model1.res2 <- model1.res[-sel,]
# 
# coord.nb <- tri2nb(my.coord.nd, row.names = NULL)
# my.coord.listw <- nb2listw (coord.nb, glist=NULL, style="W",
#                             zero.policy=FALSE)
# 
# 
# moran.test(as.vector(model1.res2), my.coord.listw)


#### Summary of models
## Pierre's idea!!

modPred_res <-as.data.frame(rbind(
  round(summary(mod1)$p.table, digits=4)[-1,],
  round(summary(mod2)$p.table, digits=4)[-1,],
  round(summary(mod3)$p.table, digits=4)[-1,]))
  

# modPred_res$predictor<- rep(c("grid_bio1", "grid_bio7", "grid_bio12", "grid_bio15", "Ruggedeness", "Erosion", "Relief", 
#                               "Historical Human Density", "Historic cropland", "Water T", "pH", "Conductivity", "TP"))
#modPred_res$predictor <- rownames(modPred_res)[1:15] 
modPred_res$predictor <- rep(c("MAP", "T.season", "P.season", "Ruggedness", "Soil type", "Lake area", "Historical human density",
                              "Historical cropland area","Water T", "pH", "Conductivity", "TP", "Chloride", "Sodium"))
modPred_res$expl.var <- rep(c("LCBD", "LCBDrepl", "LCBDrich"), each=14)
colnames(modPred_res)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")

modPred_res$coefficient<- as.numeric(as.character(modPred_res$coefficient))
modPred_res$SE<- as.numeric(as.character(modPred_res$SE))
modPred_res$sig <- "ns"
modPred_res$sig[modPred_res$P_value < 0.1] <- "P<0.1"
modPred_res$sig[modPred_res$P_value < 0.05] <- "P<0.05"
modPred_res$sig[modPred_res$P_value < 0.01] <- "P<0.01"

str(modPred_res)
modPred_res$predictor <- factor(modPred_res$predictor, 
                      levels = c("Water T", "pH", "Conductivity", "TP", "Chloride", "Sodium",
                                 "MAP", "T.season", "P.season",
                                 "Ruggedness", "Soil type", "Lake area", 
                                 "Historical human density","Historical cropland area"))


#write.csv(modPred_res, file="results/modPred_results_gamLCBD.csv")

pltLCBD <- ggplot(modPred_res, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0, linetype=2)+
  facet_wrap(~expl.var,nrow=1, scales='free_y')+
  #scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  scale_color_brewer(palette='RdBu')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=1), axis.title.x=element_blank())
pltLCBD

ggsave("figures/gamLCBD_modplot_DEF.png", plot=pltLCBD, height=8, width=10,units="in",
       dpi = 400)




