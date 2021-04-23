#load libraries for functions used
library(tidyverse)
library(ggplot2)
library(adespatial)
library(mgcv)
source("R scripts/functions.R")


#Spatio-temporal LCBD models
coordLakes <- read.csv("results/lakes_predictors.csv", row.names=1) %>%
  dplyr::select(Row.names, Longitude, Latitude)
colnames(coordLakes) <- c("sites", "Longitude", "Latitude")

#add latlong coordinates for missing lakes
coordLakes <- add_row(coordLakes, sites="EpNGEO-J_Piñan1", Longitude = -78.444, Latitude = 0.506)
coordLakes <- add_row(coordLakes, sites="GslngBA_Umayo", Longitude = -70.152, Latitude = -15.720)
coordLakes <- add_row(coordLakes, sites="Titicac_Ttcc-l-13b", Longitude = -69.521, Latitude = -15.655)

#filter core lakes
lake_cores <- paste(c("EpNGEO-J_Piñan1", "GslngBA_Umayo", "EpNGEO-J_Llaviucu", "EpNGEO-F_Yahurcch", "EpNGEO-J_Fondocch1", "Titicac_Ttcc-l-13b"), collapse = '|')
lake_cores2 <- coordLakes %>% filter(str_detect(sites, lake_cores))

#rename observations
lake_cores2[,1] <- c("Yahuarcocha", "Fondococha", "Llaviucu", "Piñan", "Umayo", "Titicaca")
colnames(lake_cores2) <- c("lake", "Longitude", "Latitude")
lake_cores2$lake_f <- factor(lake_cores2$lake, levels = c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu", "Umayo", "Titicaca"))


#read diatom core datasets
mergedCores <- read.csv("data/mergedCores_counts4.csv")[-1] #this is a dataframe with absolute counts containing all the spp

agedepth <- mergedCores[, names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- mergedCores[, !names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0

diatoms_save <- cbind(agedepth, diat)

changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names

#this is to transform to tidy format, calculate % and subset more common species
new <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, taxa, lake, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  filter(upper_age<2000) %>%
  ungroup() %>%
  group_by(depth, lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- new %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 3) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- new %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide
core_counts_wide <- core_counts_common %>%
  dplyr::select(depth, lake, upper_age, lower_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) 

## split cores by lakes and reassemble
coresList <- split(core_counts_wide, core_counts_wide$lake)


# function to calculate time series of LCBD indices and Spearman rho between indices
LCBDts <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  core[is.na(core)] <- 0
  
  core <- core[, colSums(core) > 0] #select only present species
  speciesrichness <- apply(core>0,1,sum)[-1]
  
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

  upper_age <- coresList[[i]]$upper_age
  lower_age <- coresList[[i]]$lower_age
  elapsedTime <- abs(upper_age - lower_age)
  age <- cbind(upper_age, lower_age, elapsedTime)[-1,]
  elapsed <- age[,3] #elapsed time nc==3
  diss <- diss[row(diss) == col(diss) + 1]
  roc <- diss/elapsed 
  #roc <- roc*100
  cor<-cor.test(LCBDrepl, LCBDrich, method = "spearman")$estimate
  cor.p<-cor.test(LCBDrepl, LCBDrich, method = "spearman")$p.value
  cor2<-cor.test(LCBDrich, speciesrichness, method = "spearman")$estimate
  cor2.p<-cor.test(LCBDrich, speciesrichness, method = "spearman")$p.value
  
  cbind.data.frame(age,diss,roc,LCBD,LCBDrepl,LCBDrich,repl,rich,speciesrichness,cor,cor.p,cor2,cor2.p) #combine extracted columns and remove first row to match with scd
}

## apply function to each core
coresLCB <- lapply(seq_along(coresList), LCBDts, cores=coresList)
names(coresLCB) <- names(coresList)

par(mfrow=c(3,3))
#par(mar=c(2,2,2,2))
for (i in 1:length(coresLCB)) {
  plot(coresLCB[[i]][,c("speciesrichness","LCBD")], main=names(coresLCB[i]),
       xlab="species richness", ylab="LCBD")
}

## correlate LCBD and richness for each core
richnessLCBD <- function(i, cores, ...){
  core <- cores[[i]]
  core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  core[is.na(core)] <- 0
  core <- core[, colSums(core) > 0] #select only present species
  LCBD <- beta.div(core, method = "hellinger", nperm = 999, adj = TRUE)
  LCBD <- LCBD[["LCBD"]]
  A <- beta.div.comp(core, coef = "S", quant = TRUE)
  repl <- as.matrix(A$repl)
  
  LCBDrepl <- LCBD.comp(repl, sqrt.D = TRUE)
  LCBDrepl <- LCBDrepl[["LCBD"]]
  
  upper_age <- coresList[[i]]$upper_age
  
  richness <- apply(core>0,1,sum)
  cor<-cor.test(LCBD, richness, method = "spearman")$estimate
  cor.p<-cor.test(LCBD, richness, method = "spearman")$p.value
  cbind.data.frame(upper_age,richness,LCBDrepl,LCBD,cor, cor.p)
}

richnesscore <-lapply(seq_along(coresList), richnessLCBD, cores=coresList)
names(richnesscore) <- names(coresList)

# LCBD vs richness relationship
par(mfrow=c(3,3))
par(mar=c(2,2,2,2))
for (i in 1:length(richnesscore)) {
  plot(richnesscore[[i]][,c(1,3)], main=names(richnesscore[i]))
}
  
#extract dataframes from list
LCBD_richness_data <- plyr::ldply(richnesscore, data.frame) %>%
  filter(.id %in% c("fondococha", "llaviucu", "pinan", "titicaca", "umayo", "yahuarcocha")) %>%
  mutate(lake=factor(.id)) %>%
  #left_join(lake_cores2, by="lake") %>% #here join latlong coordinates
  mutate(lake=str_replace(lake,"fondococha","Fondococha")) %>%
  mutate(lake=str_replace(lake,"yahuarcocha","Yahuarcocha")) %>%
  mutate(lake=str_replace(lake,"pinan","Piñan")) %>%
  mutate(lake=str_replace(lake,"llaviucu","Llaviucu")) %>%
  mutate(lake=str_replace(lake, "umayo", "Umayo")) %>%
  mutate(lake=str_replace(lake,"titicaca","Titicaca")) %>%
  mutate(lake_f=as.factor(lake))

##plot
LCBDcorr <- ggplot(data = LCBD_richness_data) +  
  geom_point(aes(x=richness, y = LCBD, col=upper_age)) +
  facet_wrap(factor(lake,levels=c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu",
                                  "Umayo", "Titicaca")) ~ ., ncol = 2, scales="free") +
  labs(colour="Cal yr BP")+
  xlab("Species richness")+
  theme_bw()
  
ggsave("figures/LCBD_rich_timeseries.png", LCBDcorr, height = 8, width = 10)
  
#extract dataframes from list
LCBD_ts_data <- plyr::ldply(coresLCB, data.frame) %>%
  filter(.id %in% c("fondococha", "llaviucu", "pinan", "titicaca", "umayo", "yahuarcocha")) %>%
  mutate(lake=factor(.id)) %>%
  #left_join(lake_cores2, by="lake") %>% #here join latlong coordinates
  mutate(lake=str_replace(lake,"fondococha","Fondococha")) %>%
  mutate(lake=str_replace(lake,"yahuarcocha","Yahuarcocha")) %>%
  mutate(lake=str_replace(lake,"pinan","Piñan")) %>%
  mutate(lake=str_replace(lake,"llaviucu","Llaviucu")) %>%
  mutate(lake=str_replace(lake, "umayo", "Umayo")) %>%
  mutate(lake=str_replace(lake,"titicaca","Titicaca")) %>%
  mutate(lake_f=as.factor(lake))
  
cols <- c("LCBDrepl"="#f04546","LCBDrich"="#3591d1")

### read South America hydroclimate records
pumacocha <- read.csv("data/SouthAmerica_hydroclimate_records.csv", row.names = 1) %>%
  filter(core %in% c("pumacocha")) %>% filter(!d13C==-9999.00)

plt_pumac <- ggplot(pumacocha) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(0,2000)) +
  theme(axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,1.7,0,0.6),"cm")) +
  theme(legend.position = "none") +
  ggtitle("Pumacocha") +
  xlab("Cal years BP") +
  ylab(expression(paste (delta^18, "O \u2030")))
plt_pumac

## plot the data
#create a scaling factor to apply to the second y axis
scaleFactor <- max(LCBD_ts_data$LCBDrepl) / max(LCBD_ts_data$LCBDrich)

beta_plt <- ggplot(LCBD_ts_data) +
  geom_line(aes(x=upper_age, y=LCBDrepl, colour="LCBDrepl")) +
  geom_line(aes(x=upper_age, y=LCBDrich*scaleFactor, colour="LCBDrich")) +
  scale_y_continuous(name="LCBDrepl", sec.axis=sec_axis(~./scaleFactor, name="LCBDrich"))+
  facet_wrap(factor(lake,levels=c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu",
                                     "Umayo", "Titicaca")) ~ ., ncol = 1, scales="free_y") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  scale_colour_manual(name="",values=cols) +
  ylab("") + xlab("Cal years BP") +
  theme(legend.position = "top")
beta_plt

# create a composite time series plot for LCBD and paleoclimatic record (Figure 3)
plt <- plot_grid(beta_plt, plt_pumac, rel_heights = c(4,1), ncol=1)
plt

ggsave("figures/Fig3_LCBDtimeseries_paleoclimate.png", plt, height = 8, width = 10)


##plot LCBDrepl LCBD rich
LCBDrepl_rich <- ggplot(data = LCBD_ts_data, aes(x=LCBDrepl,y=LCBDrich, col=upper_age)) +  
  geom_point() +
  facet_wrap(factor(lake,levels=c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu",
                                  "Umayo", "Titicaca")) ~ ., ncol = 2, scales="free") +
  labs(color = "Cal yr BP")+
  theme_bw()
LCBDrepl_rich

ggsave("figures/LCBDrepl_rich.png", LCBDrepl_rich, height = 8, width = 10)

## model LCBD across time and space
library(mgcv)
library(gratia)
library(cowplot)

set.seed(10) #set a seed so this is repeatable

# model GI HGAM: a global smooth plus group-level smoothers with different smootheness for each taxa
model_gam_GI <- gam(LCBD ~ s(upper_age) +
                    s(upper_age, by=lake_f, bs="fs"),
                    weights = elapsedTime / mean(elapsedTime),family=Gamma(link ="log"),
                    data=LCBD_ts_data, method="REML")
summary(model_gam_GI)
gam.check(model_gam_GI)
draw(model_gam_GI)

# #model GS HGAM: a global smooth plus group-level smoothers having the same wigliness 
model_gam_GS <- gam(LCBD ~ s(upper_age) +
                     s(upper_age, lake_f, bs="fs"),
                     weights = elapsedTime / mean(elapsedTime),family=Gamma(link = "log"),
                     data=LCBD_ts_data, method="REML")

summary(model_gam_GS)
gam.check(model_gam_GS)
draw(model_gam_GS)

#qqplot, using gratia's qq_plot function, with simulated confidence intervals
pltGS <- qq_plot(model_gam_GS, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y="Deviance residuals")
pltGI <- qq_plot(model_gam_GI, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y=NULL)

qq_plt <- plot_grid(pltGS, pltGI, ncol = 2, align = "hv", axis = "lrtb",labels="auto")
ggsave("figures/qqplots_GS_GI_new.png", qq_plt, height = 6, width = 8)

#Compare different model fits using AIC
AIC_table <- AIC(model_gam_GI, model_gam_GS)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))
AIC_table

#Create synthetic data to predict over a range of ages
lake_plot_data <- with(LCBD_ts_data, as_tibble(expand.grid(upper_age = seq(min(LCBD_ts_data$upper_age), 
                                                                           max(LCBD_ts_data$upper_age)),
                                                           lake_f = levels(LCBD_ts_data$lake_f))))

lake_modGS_fit <- predict(model_gam_GS,
                         newdata = lake_plot_data,
                         se.fit = TRUE,
                         exclude = "s(lake_f,upper_age)")


lake_modGI_fit <- predict(model_gam_GI, 
                         newdata = lake_plot_data,
                         se.fit = TRUE,
                         exclude = "s(lake_f,upper_age)")

#lake_plot_data$modS_fit <- as.numeric(lake_modS_fit$fit)
#lake_plot_data$modI_fit <- as.numeric(lake_modI_fit$fit)

lake_plot_data$modGS_fit <- as.numeric(lake_modGS_fit$fit)
lake_plot_data$modGI_fit <- as.numeric(lake_modGI_fit$fit)

# 
lake_plot_data <- gather(lake_plot_data, key=model, value=fit, modGS_fit, modGI_fit)
lake_plot_data <- mutate(lake_plot_data, se=c(as.numeric(lake_modGS_fit$se.fit),
                                              as.numeric(lake_modGI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

    
#Plot the model output for each lake, with means plus standard deviations for each model.
plot_model_labels <- paste("Model", c("GS", "GI"))
plot_model_labels <- factor(plot_model_labels, levels = plot_model_labels)

#Plot Figure 4
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

lake_plot <- ggplot(lake_plot_data)+
  facet_wrap(factor(lake_f,levels=c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu",
                                  "Umayo", "Titicaca")) ~ ., nrow = 2, scales="free_y") +
  geom_ribbon(aes(x=upper_age,
                  ymin = lower,
                  ymax = upper,
                  fill=model),
              alpha=0.2)+
  geom_point(data=LCBD_ts_data, aes(x = upper_age, y = LCBD), size=0.06) +
  geom_line(aes(x = upper_age, y = fit, color = model))+
  labs(y = "Beta replacement", x = "Age (cal yr BP)")+
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = plot_model_labels)
  
lake_plot

ggsave("figures/Fig4_HGAM_betarepl_GS_GI.png", lake_plot, height = 8, width = 10)



#This function calculates the deviance of out-of-sample data,
#conditional on their mean predicted value from the model
#code from Pedersen et al. 2019
get_deviance <- function(model, y_pred, y_obs, weights = NULL){
  stopifnot(length(y_obs)==length(y_pred))
  #We don't use the weights term in this paper, but it can be useful if
  #how well the model matters more for some sample points than others
  if(is.null(weights)) weights = rep(1, times= length(y_obs))
  #this uses the deviance residual function from the model family to
  #calculate deviances for individual points
  dev_residuals = model$family$dev.resids(y_obs, y_pred, weights)
  return(sum(dev_residuals))
}

# we need to compare how well this model fits with a null model. here we'll use
# an intercept-only model
lake_mod0 <- gam(repl~s(lake_f, bs="re"),
                     data=LCBD_ts_data,
                     family=Gamma(link ="log"),
                     weights = elapsedTime / mean(elapsedTime),
                     method="REML")


lake_predictive_summary <- LCBD_ts_data %>%
  mutate(
    #get out-of-sample predicted fits
    mod0 = as.numeric(predict(lake_mod0,.,type="response")),
    modGS = as.numeric(predict(model_gam_GS,.,type="response")),
    modGI = as.numeric(predict(model_gam_GI,.,type="response")))%>%
    group_by(lake_f)%>%
    summarise(`Intercept only` = format(get_deviance(lake_mod0, 
                                                   mod0, 
                                                   repl), 
                                      scientific = FALSE, 
                                      digits=2),
            `Model GS` = format(get_deviance(model_gam_GS, 
                                             modGS, 
                                             repl), 
                                scientific = FALSE, 
                                digits=2),
            `Model GI` = format(get_deviance(model_gam_GI, 
                                             modGI, 
                                             repl), 
                                scientific = FALSE, 
                                digits=2))%>%
    rename(Lake = lake_f)


write.csv(lake_predictive_summary, "results/HGAM_repl_GS_GI.csv")
write.csv(lake_predictive_summary, "results/HGAM_rich_GS_GI.csv")
write.csv(lake_predictive_summary, "results/HGAM_LCBDrepl_GS_GI.csv")
write.csv(lake_predictive_summary, "results/HGAM_LCBDrich_GS_GI.csv")




