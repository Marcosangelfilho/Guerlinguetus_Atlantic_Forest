#### MODIFIED FROM PRADO ET AL., 2022 (https://doi.org/10.1093/biolinnean/blab132) ########

### GUERLINGUETUS ATLANTIC FOREST
## Loading packages
library(raster)
library(vegan)


# -------------------------------------- #


######### INPUT FILES #########
# 1- Environmental layers
# 2- Mask for the study area


######### MASK #########
setwd(" ") # Set working directory with the mask for the study area 

# Load the mask for the study area
mask_AF<- shapefile("Masking_guerling.shp")
crs(mask_AF) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(mask_AF)


######### ENVIRONMENTAL LAYERS #########
## Cut bioclimatic layers by mask 

# Set working directory 
setwd(" ")
dir.create("AF_cut")

# Load bioclimatic layers
# Set working directory with the layers from the present
setwd(" ")
layers <- list.files(, pattern = '.tif')

outpath_AF <- "..../AF_cut"

outfiles_AF <- paste0(outpath_AF, layers)

for(i in 1:length(layers)) {
  e <- extent(mask_AF)
  r <-raster(layers[i])
  rc <- crop(r, e)
  rc <- extend(rc,mask_AF)
  rc<- mask(rc, mask_AF)
  rc <- writeRaster(rc, outfiles_AF[i], format = "ascii", overwrite=TRUE)
}

### Re-scale variable when needed after the clip
# Set working directory where the cut variables are
setwd ("./AF_cut/")

layers <- list.files(, pattern='asc', full.names=TRUE )

r1<- raster(layers[12]) #choose variable to be re-scaled
crs(r1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
r2<- raster(layers[13]) #choose one variable in the correct scale
crs(r2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

r1resampled <- projectRaster(r1,r2,method = 'ngb')

writeRaster(r1resampled, "layer_resc.asc")

######### SELECTING VARIABLES FOR ENMS #########
######### PCA #########

# Set working directory where the cut variables are
setwd ("~/AF_cut/")

# Calculating the Correlation among layers 
# Load all layers
layers_AF = stack(sapply(list.files(pattern='asc$', recursive = F), raster))

# PCA of the layers
raw_values_AF = values(layers_AF)
raw_values_AF = na.omit(raw_values_AF) #remove the NAs
raw_values_AF = as.data.frame(raw_values_AF)
head(raw_values_AF)
summary(raw_values_AF) 

# Standardize data 
values_trans_AF = decostand(raw_values_AF, method="standardize")
summary(values_trans_AF)

# PCA
pca_AF = prcomp (values_trans_AF)

# % PC
std_list_AF = pca_AF$sdev^2
std_list_pct_AF = std_list_AF / sum (std_list_AF) * 100
std_list_pct_sum_AF = round(std_list_pct_AF, 2)
write.csv(std_list_pct_sum_AF, file = "PC_variance_AF.csv")

contr_AF= as.data.frame(pca_AF$rotation)
write.csv(contr_AF, file = "variables_PCA_AF.csv")

# Look the number of PCs that explains 90% of the variation

n.pc_AF = length(std_list_pct_sum_AF)
sum_pca_AF = matrix(NA, n.pc_AF, 2)
for (i in 1:n.pc_AF) {sum_pca_AF[i,] = c(i, sum(std_list_pct_sum_AF[1:i]))}
colnames(sum_pca_AF) = c("PCs","% Variance")
sum_pca_AF
pc_AF = 5 ## CHANGE HERE DEPENDING OF THE RESULT
write.csv(sum_pca_AF, file = "PC_variance_sum_AF.csv", row.names = F)


## Within the PCs chosen above, choose variables which values are above the fixed value of 
# 0.32 (Based on Dormann et al., 2012 and Wang et al. 2016)


tab_AF = abs(contr_AF) 
lista_AF = list()
for (i in 1:pc_AF)
{
  linhas_AF=tab_AF[tab_AF[, i] > 0.32, ]
  lista_AF[[i]]=row.names(linhas_AF)
}

linhas.resultado_AF=unlist(lista_AF)
linhas.resultado_AF=unique(linhas.resultado_AF)
write.table(linhas.resultado_AF, file = "Variables_10%_PCA_AF.csv", sep = "\n", row.names = F, col.names = F)


################################# CORRELATION ############################################
################################# Atlantic Forest ########################################

# Correlation analysis to select only the variables with less than 0.7 of correlation 

# Load the variables selected above (linhas.resultado_AF)

setwd ("~/AF_cut/")
## Correlations among layers:
# Load all layers
predictors_AF = stack(".asc",".asc", ...) # Choose here according to the results (linhas.resultado_AF)


# Calculate correlation:
correlation_AF = layerStats(predictors_AF, "pearson", na.rm = T)
correlation_AF = correlation_AF$`pearson correlation coefficient`
correlation_AF = round(correlation_AF,2)

write.csv(correlation_AF, file = "correlation_AF.csv")

# -------------------------------------- #

