### GUERLINGUETUS ATLANTIC FOREST
## Loading packages
library("vegan")
library("ade4")
library("factoextra")
library("dplyr")

# -------------------------------------- #

######### INPUT FILES #########
# 1- Dataset with skull measurements, imputed data, and Geometric Mean

######### PRINCIPAL COMPONENT ANALYSIS #########
setwd(" ") # Set working directory

data <- read.csv("Guerlinguetus_AF_imputed_dataset.csv")

# Selecting the numeric columns
data_pca = data[,4:16]
nvariable = length(data_pca)
nind = nrow(data_pca)
sapply(data_pca, class)
head(data_pca)
summary(data_pca)

# Normalization and standardization of data 
data_pca = log10(data_pca)
data_pca = decostand (data_pca, method="standardize")
summary(data_pca)

## PCA
output = data_pca
pca_size = dudi.pca(output, center = TRUE, scannf = FALSE, nf = 13)

# Contribution of each PC
PC1 = paste0("PC1 (",round(pca_size$eig[1]/sum(pca_size$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_size$eig[2]/sum(pca_size$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_size$eig[3]/sum(pca_size$eig)*100,2),"%)")
PC4 = paste0("PC4 (",round(pca_size$eig[4]/sum(pca_size$eig)*100,2),"%)")
PC5 = paste0("PC5 (",round(pca_size$eig[5]/sum(pca_size$eig)*100,2),"%)")

# Saving the results
pc_names <- paste0("PC", 1:nvariable) # PCs names
var_names = colnames(data_pca) # Variables names
sample_names = data$Voucher # Sample names

# PC contribution
perc_pca = get_eigenvalue(pca_size)
rownames(perc_pca) = pc_names
colnames(perc_pca) = c("Eigenvalues", "Variance (%)", "Cumulative Variance (%)")
head(perc_pca)
write.csv(perc_pca, file = "PCA_eigenvalues.csv")

# Loadings
var = get_pca_var(pca_size)
var_loadings = var$coord
colnames(var_loadings) = pc_names
write.csv(var_loadings, file = "PCA_variable_loadings.csv")

# Individual Scores
ind_scores = pca_size$li
rownames(ind_scores) = sample_names
colnames(ind_scores) = pc_names
write.csv(ind_scores, file = "PCA_scores.csv")

# -------------------------------------- #

######### GEOMETRIC MEAN DIFFERENCES #########

# Testing Normality
data %>%
  group_by(Pooled.sample) %>%
  filter(n() >= 3) %>%
  summarise(
    W = shapiro.test(GM)$statistic,
    p_value = shapiro.test(GM)$p.value
  )

# Levene's test
levene_test <- car::leveneTest(GM ~ Pooled.sample, data = data)
print(levene_test)

# ANOVA
anova_result <- aov(GM ~ Pooled.sample, data = data)
summary(anova_result)

# POST HOC TUKEY THSD
tukey_result <- TukeyHSD(anova_result, "Pooled.sample")
print(tukey_result)

# -------------------------------------- #
