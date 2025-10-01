### GUERLINGUETUS ATLANTIC FOREST
## Loading packages
library(ecodist)
library(ade4)
library(fossil)
library(raster)
library(RStoolbox)
library(vegan)


# -------------------------------------- #


######### INPUT FILES #########
# 1- Individual coordinates
# 2- PCA scores
# 3- Environmental variables

######### MANTEL TEST #########
setwd(" ") # Set working directory

scores <- read.csv("PCA_scores.csv")
rownames(scores) <- scores[[1]]
scores <- scores[ , -1]
scores <- scores[,1:3] # select the PCs that explain 70% of the variance
coords <- read.csv ("Guerlinguetus_coords.csv")
rownames(coords) <- coords[[1]]
coords <- coords[ , -1]

# Preparing the Mahalanobis and greographical distance matrices
dist_scores <- ecodist::distance(scores, method = "mahalanobis") 
dist_coords <- earth.dist(coords, dist=TRUE) 

# Running Mantel
mantel_test <- mantel.rtest(dist_scores, dist_coords, nrepet = 10000)
print(mantel_test)


# -------------------------------------- #


######### dbRDA #########
### Geographical distance
dist_coords <- as.dist(dist_coords)
dist_coords_vec <- as.vector(dist_coords)
hist(dist_coords_vec)
dist_coords2<-as.matrix(dist(dist_coords))

pcnm_coords<-pcnm(dist_coords2, dist.ret = FALSE)

### Mahalanobis distance
dist_scores2 <- as.matrix(dist(dist_scores))

### Environmental variables
setwd(" ") # set working directory with the environmental variables for your area

predictors <- stack(sapply(list.files(pattern='\\.tif$', recursive = F), raster))
names(predictors)

#Perform the PCA:
pca_env <- rasterPCA(predictors, nComp = 1, spca = TRUE)
summary(pca_env$model) # PC1 50% of the total variance
pca_env$model$loadings

PC1_mor<-extract(pca_env$map$PC1, coords)

write.csv(PC1_mor, "PC1_env_mor.csv") # save the PCA env distance

######### Redudancy Analysis #########
dbrda4<-capscale(dist_scores2 ~ PC1_mor$PC1, dist = "man") 
summary(dbrda4)
sig4<- anova(dbrda4)
sig4

dbrda5<-capscale(dist_scores2 ~ pcnm_coords$vectors, dist = "man") 
dbrda5
summary(dbrda5)
sig5<- anova(dbrda5)
sig5

dbrda6<-capscale(dist_scores2 ~ PC1_mor$PC1 + Condition(pcnm_coords$vectors), dist = "man") 
summary(dbrda6)
sig6<- anova(dbrda6)
sig6

# -------------------------------------- #
