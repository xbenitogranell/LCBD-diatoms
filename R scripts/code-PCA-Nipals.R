

library(ade4) #to perfom PCA NIPALS analysis
library(mice) #to compute missing values (NA) for multivariate ordination

#read data
limno <- env_data_lakes[,3:21]
limno <- limno[, -which(names(limno) %in% c("Turb", "Chl", "Alkalinity", "Secchi", "Si",
                                            "NO2", "PO4", "TN", "NO3", "TP"))]

#read lake predictors (from code-LCBD script)
lake_predictors <- read.csv("results/lakes_predictors.csv", row.names=1)

#Check NAs and remove rows with all NAs
row.all.na <- apply(limno, 1, function(x) all(is.na(x)))
sum(row.all.na)
limno <- limno[ !row.all.na, ]


limno <- transform(limno, Water.T=log10(Water.T+0.25), Cond=log10(Cond+0.25), Ca=log10(Ca+0.25), 
                   Mg=log10(Mg+0.25), K=log10(K+0.25), Na=log10(Na+0.25), Cl=log10(Cl+0.25), 
                   SO4=log10(SO4+0.25))


limno_subset <- limno[!rownames(limno) %in% c("Brdb-SA_SlcrTrqm", "Brdb-SA_LgBS.RJU1",
                                    "Chl-Crr_LagunCrn"),]

PCA.data <- limno_subset
#PCA.data <- as.data.frame(scale(limno, center = TRUE, scale=TRUE))

#Compute the PCA
#PCA NIPALS with ade4
PCA.nipals <- nipals(PCA.data, nf = 4, rec = FALSE, niter = 1000, tol = 1e-09)


#Save summary matrix results  
PCA.summary <- NULL
PCA.summary <- matrix(data = NA, nrow = 8, ncol = 3, byrow = FALSE, dimnames = NULL)
colnames(PCA.summary) <- c("Parameter", "Value", "Explained variance")


PCA.summary[1, 1] <- c("Number of extracted factors")
PCA.summary[1, 2] <- PCA.nipals$nf

PCA.summary[2, 1] <- c("Number of variables")
PCA.summary[2, 2] <- nrow(PCA.nipals$co)

PCA.summary[3, 1] <- c("Eigenvalue 1")
PCA.summary[3, 2] <- PCA.nipals$eig[1]
PCA.summary[3, 3] <- (PCA.nipals$eig[1] / nrow(PCA.nipals$co)) * 100

PCA.summary[4, 1] <- c("Eigenvalue 2")
PCA.summary[4, 2] <- PCA.nipals$eig[2]
PCA.summary[4, 3] <- (PCA.nipals$eig[2] / nrow(PCA.nipals$co)) * 100

PCA.summary[6, 1] <- c("Number of iterations axis 1")
PCA.summary[6, 2] <- PCA.nipals$nb[1]

PCA.summary[7, 1] <- c("Number of iterations axis 2")
PCA.summary[7, 2] <- PCA.nipals$nb[2]

write.csv(PCA.summary, file = "PCA summary table.csv", sep = "\t", col.names = TRUE, row.names=FALSE, na="")

#4 components
#Save the column coordinates (Component matrix) #ncol= n variables 
Component.matrix <- NULL
Component.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 4, byrow = FALSE, dimnames = NULL)
Component.matrix[,1] <- rownames(PCA.nipals$co)
Component.matrix[,2] <- (((PCA.nipals$co[,1] ^ 2) * PCA.nipals$eig[1]) / PCA.nipals$eig[1]) + (((PCA.nipals$co[,2] ^ 2) * PCA.nipals$eig[2]) / PCA.nipals$eig[2])
Component.matrix <- cbind(Component.matrix[,1], PCA.nipals$co, Component.matrix[,2])
colnames(Component.matrix) <- c("Variable", "Component 1", "Component 2", "Component 3", "Component 4", "Power Importance")
#unlink("PCA component scores.txt", recursive = FALSE)
write.csv(Component.matrix, file = "PCA component scores.csv", sep = "\t", col.names = TRUE, row.names=FALSE, na="")

#4 components    
#Save the row coordinates (Factor Scores)
colnames(PCA.nipals$li) <- c("Component 1", "Component 2", "Component 3", "Component 4")
Factor.scores <- data.frame(cbind(PCA.data, PCA.nipals$li))
#unlink("PCA factor scores.txt", recursive = FALSE)
write.csv(Factor.scores, file = "PCA factor scores.csv", sep = "\t", col.names = TRUE, row.names=FALSE, na="")


#4 components  
#Save the column normed scores (Component scores coefficient matrix)
Component.coefficient.matrix <- NULL
Component.coefficient.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 1, byrow = FALSE, dimnames = NULL)
Component.coefficient.matrix[,1] <- rownames(PCA.nipals$c1)
Component.coefficient.matrix <- cbind(Component.coefficient.matrix, PCA.nipals$c1)
colnames(Component.coefficient.matrix) <- c("Variable", "Component 1", "Component 2", "Component 3", "Component 4")
#unlink("PCA component scores coefficient matrix.txt", recursive = FALSE)
write.csv(Component.coefficient.matrix, file = "PCA component scores coefficient matrix.csv", sep = "\t", col.names = TRUE, row.names=FALSE, na="")


#Plotting PCA

par(mfrow=c(1,2))
par(mar=c(3,3,2,2)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text. 

#Factor scores (samples)

#Create data frame with factor scores, regions and site groupings
PCA.scores <- data.frame(component1=Factor.scores$Component.1, component2=Factor.scores$Component.2)


plot(PCA.scores$component1, PCA.scores$component2, type = "n", xlab = "PCA1", ylab = "PCA2")
abline(h=0, col="grey")
abline(v=0, col="grey")

#
points(PCA.scores, pch=20)
text(PCA.scores$component1, PCA.scores$component2, labels = rownames(Factor.scores), pos = 1, cex = 0.5, offset = 0.2)


#plotting variables
  comp1 <- as.numeric(Component.coefficient.matrix[,2])
  comp2 <- as.numeric(Component.coefficient.matrix[,3])

#Labels
  labels <- as.character(Component.coefficient.matrix[,1])

  plot(comp1, comp2, pch=16, col="black", xlab = "PCA1", ylab = "PCA2")
  abline(h=0, col="grey")
  abline(v=0, col="grey")

  text(comp1, comp2, labels = labels, pos = 1, cex = 0.5, offset = 0.2)

## Merge mutiple datasets
PCA.factors <- Factor.scores[,which(names(Factor.scores) %in% c("Component.1", "Component.2", "Component.3", "Component.4"))]   
rownames(lake_predictors) <- lake_predictors$sites
  
lake_predictors2 <- subset(lake_predictors, lake_predictors$sites %in% rownames(PCA.factors))

t <- merge(lake_predictors2, PCA.factors, by=0, all=TRUE)
rownames(t) <- t$Row.names
write.csv(t, "results/lakes_predictors.csv")

#subset diatom trainingset
t2 <- subset(training, rownames(training) %in% rownames(t))
rowSums(t2)
write.csv(t2, "results/diatoms.csv")




  
