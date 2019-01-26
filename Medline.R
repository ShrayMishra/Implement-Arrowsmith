# Implementation of “Quantitative model for linking two disparate sets of article in MEDLINE” research 
# from Bioinformatics journal using logistic regression

# By Shray Mishra 


#Installing packages 

#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("reshape")

#Importing libraries for assessment

library(readxl)
library(ggplot2)

# Loading Arrowsmith dataset in R
# col.names parameter is used inside read.csv which will create vectors. 

arwsmt <- read.csv("Arrowsmith.csv", header=TRUE, skip = 4,
                   col.names=c("Arrowsmith_search", "A_litsize", "C_litsize", "B_term", "target",
                               "nA", "nC", "nof_MeSH", "nofsemantic_cat",
                               "cohesion_score", "n_in_MEDLINE", "Firstyr_MEDLINE", "pAC",
                               "onmedium_stoplist.", "onlong_stoplist.", "X1", "X2", "X3", "X4",
                               "X5", "X6", "X7", "I1", "I2", "I3", "I4", "I5", "I6", "Y"))

# Feature creation for modelling

# X1= 1 if (nA > 1 or A_litsize < 1000) and (nC > 1 or C_litsize < 1000), otherwise it is 0 
arwsmt$X1 <- ifelse(((arwsmt$nA > 1)|(arwsmt$A_litsize < 1000)) &
                      ((arwsmt$nC > 1)|(arwsmt$C_litsize < 1000)), 1, 0)

# X2 = 1 if nof_MeSH > 0 and < 99999 and  0.5 if nof_MeSH = 99999 and otherwise it is 0 
arwsmt$X2 <- ifelse((arwsmt$nof_MeSH > 0 & arwsmt$nof_MeSH < 99999), 1,
                    ifelse(arwsmt$nof_MeSH == 99999, 0.5, 0))

# X3 = 1 if nof semantic categories > 0, 0 otherwise
arwsmt$X3 <- ifelse(arwsmt$nofsemantic_cat > 0, 1, 0)

# X4 = cohesion score if cohesion score < 0.3, 0.3 otherwise
arwsmt$X4 <- sapply(arwsmt$cohesion_score, function(x) min(c(x, 0.3)))

# X5 = -|log10(n in MEDLINE) – 3|
arwsmt$X5 <- -abs(log10(arwsmt$n_in_MEDLINE)-3)

# X6 = max(min(1st year in MEDLINE,2005),1950)
arwsmt$X6 <- sapply(arwsmt$Firstyr_MEDLINE, function(x) max(c(min(c(x, 2005)), 1950)))

# X7 = min(8,-log10(pAC+0.000000001))
arwsmt$X7 <- sapply(arwsmt$pAC, function(x) min(c(8, -log10(x+0.000000001))))

# I1 = 1 if Arrowsmith search = ‘retinal detachment vs aortic aneurysm’, 0 otherwise
arwsmt$I1 <- ifelse(arwsmt$Arrowsmith_search == "retinal detachment vs aortic aneurysm", 1, 0)

# I2 = 1 if Arrowsmith search = ‘NO and mitochondria vs PSD’
arwsmt$I2 <- ifelse(arwsmt$Arrowsmith_search == "NO and mitochondria vs PSD", 1, 0)

# I3 = 1 if Arrowsmith search = ‘mGluR5 vs lewy bodies’
arwsmt$I3 <- ifelse(arwsmt$Arrowsmith_search == "mGluR5 vs lewy bodies", 1, 0)

# I4 = 1 if Arrowsmith search = ‘magnesium vs migraine’
arwsmt$I4 <- ifelse(arwsmt$Arrowsmith_search == "magnesium vs migraine", 1, 0)

# I5 = 1 if Arrowsmith search = ‘Calpain vs PSD’
arwsmt$I5 <- ifelse(arwsmt$Arrowsmith_search == "Calpain vs PSD", 1, 0)

# I6 = 1 if Arrowsmith search = ‘APP vs reelin’
arwsmt$I6 <- ifelse(arwsmt$Arrowsmith_search == "APP vs reelin", 1, 0)

# Y = 1 if target = 0 or 2, 0 otherwise
arwsmt$Y <- ifelse(arwsmt$target %in% c(0, 2), 1, 0)

#Calculating summary for arrowsmith dataset
summary(arwsmt)

# Finding Missing Values
# "na" is a missing value indicator and "which" used along with will give the position of missing values 
which(is.na(arwsmt) == TRUE) 

# Missing integer are 0

# Finding out outliers in dataset
# Plotting the variables to find outliers

scatter.smooth(arwsmt$nof_MeSH)
scatter.smooth(arwsmt$cohesion_score)
scatter.smooth(arwsmt$nA)
scatter.smooth(arwsmt$nC)
scatter.smooth(arwsmt$Firstyr_MEDLINE)


which(arwsmt$nof_MeSH == 99999 |
        arwsmt$cohesion_score > 0.6 |
        arwsmt$nA > 1000 |
        arwsmt$nC > 1000 |
        arwsmt$Firstyr_MEDLINE == 9999) 

#Plotting different variables of the database

library(ggplot2)
library(gridExtra)

# hist_fucs function for graphical presentations of histogram
hist_funcs <- function(variable){
  ggplot(data = arwsmt) +
    geom_histogram(mapping = aes_string(x = variable),
                   bins = 30) +
    ggtitle(paste0("Histogram of ", variable))
}

# Plots before feature selection
grid.arrange(hist_funcs("A_litsize"))
grid.arrange(hist_funcs("nA"))
grid.arrange(hist_funcs("C_litsize"))
grid.arrange(hist_funcs("nC"))
grid.arrange(hist_funcs("nof_MeSH"))
grid.arrange(hist_funcs("cohesion_score"))
grid.arrange(hist_funcs("n_in_MEDLINE"))
grid.arrange(hist_funcs("Firstyr_MEDLINE"))
grid.arrange(hist_funcs("pAC"))
grid.arrange(ggplot(data = arwsmt) +
               geom_bar(mapping = aes(x = Arrowsmith_search)) +
               ggtitle("BarPlot of Arrowsmith_search") + scale_x_discrete(labels = abbreviate))


# Plots after feature selection
grid.arrange(hist_funcs("X1"))
grid.arrange(hist_funcs("X2"))
grid.arrange(hist_funcs("X3"))
grid.arrange(hist_funcs("X4"))
grid.arrange(hist_funcs("X5"))
grid.arrange(hist_funcs("X6"))
grid.arrange(hist_funcs("X7"))
grid.arrange(hist_funcs("I1"))
grid.arrange(hist_funcs("I2"))
grid.arrange(hist_funcs("I3"))
grid.arrange(hist_funcs("I4"))
grid.arrange(hist_funcs("I5"))
grid.arrange(hist_funcs("I6"))

# sp function function for graphical presentations of scatterplot

sp <- function(variable1, variable2){
  ggplot(data = arwsmt) + 
    geom_point(mapping = aes_string(x = variable1,
                                    y = variable2)) +
    ggtitle(paste0(variable2, " vs ", variable1))
}

# Plotting X1 vs nA, nC, A_litsize, C_litsize
# melt is used to convert the object into molten dataframe

library(reshape)
arwsmt1 <- melt(subset(arwsmt, select = c("X1", "nA", "nC", "A_litsize", "C_litsize")), id.vars="X1")
ggplot(data = arwsmt1) +
  geom_point(mapping = aes(x = X1,
                           y = value)) +
  facet_grid(variable~.) +
  ggtitle("nA, nC, A_litsize, C_litsize vs X1")

# Scatterplot for variables
grid.arrange(sp("X2", "nof_MeSH"),
             sp("X3", "nofsemantic_cat"),
             sp("X4", "cohesion_score"),
             sp("X5", "n_in_MEDLINE"),
             nrow=2, ncol=2)

grid.arrange(sp("X6", "Firstyr_MEDLINE"),
             sp("X7", "pAC"),
             nrow=2, ncol=2)

# Fitting logistic regression model in the dataset

log_reg <- glm(Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + I1 + I2 + I3 + I4 + I5 + I6,
               data = arwsmt,
               family=binomial)
# After looking at the summary statistcs, we can ignore I6 in the model.

summary(log_reg)
