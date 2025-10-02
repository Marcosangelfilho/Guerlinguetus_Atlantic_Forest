### GUERLINGUETUS ATLANTIC FOREST
## Loading packages
library(MVN)
library(purrr)
library(pander)

# -------------------------------------- #

######### INPUT FILES #########
# 1- Dataset with skull measurements and imputed data for each pooled sample

######### CHECKING FOR NORMALITY #########
setwd(" ") # Set working directory

Litoral <- read.csv("Litoral_MW.csv")
Parana <- read.csv("Parana_MW.csv")
SerraMar <- read.csv("SerraMar_MW.csv")
SerraParanapiacaba <- read.csv("SerraParanapiacaba_MW.csv")
SouthRioDoce <- read.csv("SouthRioDoce_MW.csv")
Tiete <- read.csv("Tiete_MW.csv")


normality = mvn(data = SerraMar[,2:15], subset = "Sex", mvn_test = "mardia",
             univariate_test = "SW", 
             multivariate_outlier_method = "adj",
             show_new_data = TRUE)

## Univariate Normality Result (Shapiro-Wilk)
summary(normality, select = "univariate")

## Multivariate Normality Result (Mardia's Kurtosis)
summary(normality, select = "mvn")


######### MANN-WHITNEY #########
# Inserting "NA" on blank strings
Litoral$Sex[Litoral$Sex == ""] <- NA
Parana$Sex[Parana$Sex == ""] <- NA
SerraMar$Sex[SerraMar$Sex == ""] <- NA
SerraParanapiacaba$Sex[SerraParanapiacaba$Sex == ""] <- NA
SouthRioDoce$Sex[SouthRioDoce$Sex == ""] <- NA
Tiete$Sex[Tiete$Sex == ""] <- NA

vars <- c("SL","ZW","IO","PC","BW","NL","NW","DT","UTB","MM","MBL","MBH","LT") #  Listing variables

# wilcox.test for every variable by sex in each pooled sample
dfs <- list(
  Litoral = Litoral,
  Parana = Parana,
  SerraMar = SerraMar,
  SerraParanapiacaba = SerraParanapiacaba,
  SouthRioDoce = SouthRioDoce,
  Tiete = Tiete
)

walk2(dfs, names(dfs), ~{
  df <- .x
  dataset_name <- .y
  cat("\n\n## Dataset:", dataset_name, "\n")
  
res <- map(vars, ~ wilcox.test(df[[.x]] ~ df$Sex))
  names(res) <- vars
  
# Print results
  walk2(res, vars, ~{
    cat("\n###", .y, "\n")
    pander(.x)
  })
})


######### MANOVA #########
manova_results <- manova (cbind(SL,ZW,IO,PC,BW,NL,NW,DT,UTB,MM,MBL,MBH,LT) ~ Sex, data = SerraMar)
summary(manova_results, test = "Pillai")

# -------------------------------------- #
