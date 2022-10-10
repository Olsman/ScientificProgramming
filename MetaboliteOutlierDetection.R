# Code 04-10-2022 - METABOLITES
# Scientific Programming
# Rosan Olsman


########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

# Packages
# Required CRAN packages:
CRANpackages <- c("tidyverse", "readxl", "devtools", "ggplot2", "dplyr", "MASS")

# Required Bioconductor packages:
BiocPackages <- c("vioplot", "plotly", "pcaMethods", "limma")

# Install (if not yet installed) and load the required packages: 
for (pkg in CRANpackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (pkg in BiocPackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}

# Import data
MetaData <- data.frame(read_excel("./Dataset/MetaTumourData.xlsx"), row.names = "Mouse.ID")
MetabolitesData <- data.frame(read_excel("./Dataset/StoolMetabolites.xlsx"), check.names = FALSE)
MicrobiomeData <- data.frame(read_excel("./Dataset/OTUTable.xlsx"), row.names = "ID", check.names = FALSE)

# Sanity check - identify common and different samples for each dataset
all(colnames(MicrobiomeData) %in% rownames(MetaData)) # FALSE: one sample with wrong label - 243 is missing
setdiff(rownames(MetaData), colnames(MicrobiomeData))
all(colnames(MetabolitesData[,14:43]) %in% rownames(MetaData))
setdiff(rownames(MetaData), colnames(MetabolitesData[,14:43]))

# most likely that 233 should have been 243; change 233 to 243
setdiff(colnames(MicrobiomeData), rownames(MetaData))
names(MicrobiomeData)[names(MicrobiomeData) == '233'] <- '243'
all(colnames(MicrobiomeData) %in% rownames(MetaData))

# check for NAs in dataset - no missing values in microbiome and metabolite data
which(colSums(is.na(MetaData))>0)
which(colSums(is.na(MetabolitesData))>0)
which(colSums(is.na(MicrobiomeData))>0)


########## Basic data exploration - metabolite ########## 
# Data visualization - weight
Weight <- ggplot(MetaData, aes(x = Category, y = Weight, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Body Weight (g)",
       title = "Violin Plot: Body Weight") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Weight
# two possible outliers: 151 (WT) & 13 (WT)

# Data visualization - fat mass
FatMass <- ggplot(MetaData, aes(x = Category, y = Fat.mass, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Fat Mass (g)",
       title = "Violin Plot: Fat Mass") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
FatMass
# no clear outliers

# Data visualization - lean mass
LeanMass <- ggplot(MetaData, aes(x = Category, y = Lean.mass, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Lean Mass (g)",
       title = "Violin Plot: Lean Mass") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
LeanMass
# two possible outliers: 108 (HF) & 115 (LF)

# Data visualization - sex
Gender <- ggplot(data = MetaData, aes(x = Category, fill = Sex)) +
  geom_bar(position = position_dodge()) +
  theme_classic() +
  labs(title = "Gender per Condition", x = "Condition", y = "Count") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Gender
# boxplot shows clear class imbalance for sex

# Data visualization - tumor volume
TumorVolume <- ggplot(MetaData, aes(x = Category, y = Tumor.volume, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Tumor Volume",
       title = "Violin Plot: Tumor Volume") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
TumorVolume

# Data visualization - number of tumors
Tumors <- ggplot(MetaData, aes(x = Category, y = Tumors, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Tumours",
       title = "Violin Plot: Tumor count") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Tumors


########## Metabolite outlier removal ########## 
# PCA and outlier detection 
# scale data 
scaledMetabolites <- (MetabolitesData[,14:43] - rowMeans(MetabolitesData[,14:43]))/(apply(MetabolitesData[,14:43],1,sd))

# sanity check - scaling
mean(rowMeans(scaledMetabolites))            # close to 0
mean(sd(as.matrix(scaledMetabolites)))       # close to 1

# PCA
pcaResults <- pca(t(scaledMetabolites))
summary(pcaResults)

# index - select metadata from samples included in metabolomics data 
indexMets <- as.list(rownames(pcaResults@scores))
sub_df <- MetaData %>%
  filter(rownames(MetaData) %in% indexMets)

# sanity check
all(rownames(pcaResults@scores) %in% rownames(sub_df))

# data for PCA plots
plotData <- merge(data.frame(pcaResults@scores), sub_df, by="row.names", all.x=TRUE)
rownames(plotData) <- plotData$Row.names
plotData$Row.names <- NULL

# PCA plot
ggplot(plotData, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
       title = "Scaled metabolomics data - no preprocessing") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# no clear outliers visible


########## Isolation forest for outlier removal (without further pre-processing metabolomics data) ########## 
install.packages("IsolationForest", repos="http://R-Forge.R-project.org")
# yes - packages that need compilation
library("IsolationForest")

# obtain dataframe with scaled metabolites and ID 
metabolites <- as.data.frame(t(scaledMetabolites))
colnames(metabolites) <- MetabolitesData$COMP_ID

# building the isolation tree for outlier detection
IsolationTree <- IsolationTrees(metabolites, rFactor = 0) #rFactor = 0; meaning fully deterministic
AnomalyScore <- AnomalyScore(metabolites, IsolationTree)

# data for PCA plot including the anomaly score 
plotDataIT <- cbind(data.frame(plotData, as.data.frame(AnomalyScore[["outF"]])))
colnames(plotDataIT)[which(names(plotDataIT) == rev(names(plotDataIT))[1])] <- "AnomalyScore"

# identify an anomaly score ≥ 0.8 as an outlier
plotDataIT$Outlier <- as.factor(ifelse(plotDataIT$AnomalyScore >=0.8, "Outlier", "Normal"))
sum(plotDataIT$Outlier == "Outlier")    # 1 outlier detected

# pca using anomaly score
ggplot(plotDataIT, aes(x = PC1, y = PC2, color = AnomalyScore, shape = Category)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
       title = "Scaled metabolites data coloured by anomaly score") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# pca when setting anomaly score to > 0.8
ggplot(plotDataIT, aes(x = PC1, y = PC2, color = Outlier)) + 
  geom_point() +
  scale_color_manual(values = c("Normal" = "#00BFC4", "Outlier" = "#F8766D")) +
  labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
       title = "PCA: Metabolites coloured by Anonamly Score Treshold") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# only one outlier visible

# select normal data
inclusion <- plotDataIT %>%
  filter(Outlier != "Outlier")

index <- row.names(inclusion)
MetabolitesDataMet <- as.data.frame(t(MetabolitesData[,14:43]))
metabolitesFiltered <- MetabolitesDataMet %>%
  filter(rownames(MetabolitesDataMet) %in% index)
all(rownames(inclusion) %in% rownames(metabolitesFiltered))

## remove (?)?(?(?)?(?(?)?(?(?)?)))
#metabolitesFiltered <- as.data.frame(scale(metabolitesFiltered), scale = TRUE)

# PCA - now filtering the outliers
scaledMetabolitesFiltered <- (t(metabolitesFiltered) - rowMeans(t(metabolitesFiltered)))/(apply(t(metabolitesFiltered),1,sd))
mean(rowMeans(scaledMetabolitesFiltered)) # close to 0
mean(sd(scaledMetabolitesFiltered))       # close to 1

pcaResultsIT <- pca(t(scaledMetabolitesFiltered), nPcs = 8)
summary(pcaResultsIT)

# index 
indexMetsIT <- as.list(rownames(pcaResultsIT@scores))
sub_dfIT <- MetaData %>%
  filter(rownames(MetaData) %in% indexMetsIT)

# sanity check
all(rownames(pcaResultsIT@scores) %in% rownames(sub_dfIT))

# data for PCA plots
plotDataITexluded <- merge(data.frame(pcaResultsIT@scores), sub_dfIT, by="row.names", all.x=TRUE)
rownames(plotDataITexluded) <- plotDataITexluded$Row.names
plotDataITexluded$Row.names <- NULL

# PCA plot
ggplot(plotDataITexluded, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsIT@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsIT@R2[2] * 100, 1), "% of the variance"),
       title = "PCA plot: after outlier removal by Isolation Tree") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


########## Isolation forest for outlier removal (with pre-processing metabolomics data) ########## 
# normalize
source("PQNfunction.R")
# metabolites - run PQN function
pqnMetabolites <- pqn(as.matrix(MetabolitesData[,14:43]))

# transform
transMetabolites <- log2(pqnMetabolites)

# scale
# scale data 
scaledMetabolitesPP <- (transMetabolites - rowMeans(transMetabolites))/(apply(transMetabolites,1,sd))
mean(rowMeans(scaledMetabolitesPP)) # close to 0
mean(sd(scaledMetabolitesPP))       # close to 1

# PCA
pcaResultsPP <- pca(t(scaledMetabolitesPP))
summary(pcaResultsPP)

# index - select metadata from samples included in metabolomics data 
indexMetsPP <- as.list(rownames(pcaResultsPP@scores))
sub_dfPP <- MetaData %>%
  filter(rownames(MetaData) %in% indexMetsPP)

# sanity check
all(rownames(pcaResultsPP@scores) %in% rownames(sub_dfPP))

# data for PCA plots
plotDataPP <- merge(data.frame(pcaResultsPP@scores), sub_dfPP, by="row.names", all.x=TRUE)
rownames(plotDataPP) <- plotDataPP$Row.names
plotDataPP$Row.names <- NULL

# PCA plot
ggplot(plotDataPP, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA plot: scaled metabolomics data") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# two possible outliers


########## Outlier removal - isolatin tree with pre-processed metabolomics data ##########
# what happens if i do it on the data
metabolitesPP <- as.data.frame(t(scaledMetabolitesPP))
colnames(metabolitesPP) <- MetabolitesData$COMP_ID
IsolationTreePP <- IsolationTrees(metabolitesPP, rFactor = 0) #rFactor = 0; meaning fully deterministic
AnomalyScorePP <- AnomalyScore(metabolitesPP, IsolationTreePP)
plotDataITPP <- cbind(data.frame(plotDataPP, as.data.frame(AnomalyScorePP[["outF"]])))
colnames(plotDataITPP)[which(names(plotDataITPP) == rev(names(plotDataITPP))[1])] <- "AnomalyScore"
plotDataITPP$Outlier <- as.factor(ifelse(plotDataITPP$AnomalyScore >=0.8, "Outlier", "Normal"))
sum(plotDataITPP$Outlier == "Outlier")

# pca using anomaly score
ggplot(plotDataITPP, aes(x = PC1, y = PC2, color = AnomalyScore, shape = Category)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA: Metabolites coloured by Anomaly Score") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# pca when setting anomaly score to > 0.8
ggplot(plotDataITPP, aes(x = PC1, y = PC2, color = Outlier)) + 
  geom_point() +
  scale_color_manual(values = c("Normal" = "#00BFC4", "Outlier" = "#F8766D")) +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA: Metabolites coloured by Anonamly Score Treshold") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# only one outlier visible: sample 49

# select normal data
inclusionPP <- plotDataITPP %>%
  filter(Outlier != "Outlier")

indexPP <- row.names(inclusionPP)
metabolitesFilteredPP <- MetabolitesDataMet %>%
  filter(rownames(MetabolitesDataMet) %in% indexPP)
all(rownames(inclusion) %in% rownames(metabolitesFiltered))
# remove sample 4 and 7 because they appear to be far away from the rest
#metabolitesFilteredPP <- metabolitesFilteredPP[!(rownames(metabolitesFilteredPP) == "4" | rownames(metabolitesFilteredPP) == "7" ),]

#metabolitesFilteredPP <- as.data.frame(scale(metabolitesFiltered), scale = TRUE)

## HERE SOMETHING GOES WRONG I THINK BECAUSE I GET THE SAME PCA

# PCA - now filtering the outliers


# norm??
# trans??
# metabolites - run PQN function
pqnMetabolitesPP <- pqn(as.matrix(metabolitesFilteredPP))

# transform
transMetabolites <- log2(pqnMetabolitesPP)


#
scaledMetabolitesFiltered <- (t(transMetabolites) - rowMeans(t(transMetabolites)))/(apply(t(transMetabolites),1,sd))
mean(rowMeans(scaledMetabolitesFiltered)) # close to 0
mean(sd(scaledMetabolitesFiltered))       # close to 1

pcaResultsIT <- pca(t(scaledMetabolitesFiltered))
summary(pcaResultsIT)

# index 
indexMetsIT <- as.list(rownames(pcaResultsIT@scores))
sub_dfIT <- MetaData %>%
  filter(rownames(MetaData) %in% indexMetsIT)

# sanity check
all(rownames(pcaResultsIT@scores) %in% rownames(sub_dfIT))

# data for PCA plots
plotDataITexluded <- merge(data.frame(pcaResultsIT@scores), sub_dfIT, by="row.names", all.x=TRUE)
rownames(plotDataITexluded) <- plotDataITexluded$Row.names
plotDataITexluded$Row.names <- NULL

# PCA plot
ggplot(plotDataITexluded, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsIT@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsIT@R2[2] * 100, 1), "% of the variance"),
       title = "PCA plot: after outlier removal by Isolation Tree") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# IT IS CLEAR: ONLY SAMPLE 49 IS AN OUTLIER - PCA before normalization and transformation appears to be better

########## DATA ANALYSIS ##########
normalizedMetabolitest <- pqn(t(metabolitesFilteredPP))
transformationMetabolites <- log2(normalizedMetabolitest)
scaledMetabolitesAnalysis <- (t(transformationMetabolites) - rowMeans(t(transformationMetabolites)))/(apply(t(transformationMetabolites),1,sd))

