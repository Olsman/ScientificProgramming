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

# Sanity check - identify common and different samples for each dataset
all(colnames(MetabolitesData[,14:43]) %in% rownames(MetaData))
setdiff(rownames(MetaData), colnames(MetabolitesData[,14:43]))

# check for NAs in dataset - no missing values in microbiome and metabolite data
which(colSums(is.na(MetaData))>0)
which(colSums(is.na(MetabolitesData))>0)


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
MetaData$ID <- rownames(MetaData)
sub_df <- MetaData %>%
  filter(MetaData$ID %in% indexMets)

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
       title = "Metabolomics data - before outlier removal") +
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
       title = "Metabolomics data - Anomaly Score") +
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
       title = "Metabolomics data - Outliers") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# only one outlier visible

# select data that are not outliers
inclusion <- plotDataIT %>%
  filter(Outlier != "Outlier")

# obtain sample ID and samples from included values
index <- row.names(inclusion)
MetabolitesDataMet <- MetabolitesData[,14:43]
rownames(MetabolitesDataMet) <- MetabolitesData$COMP_ID
MetabolitesDataMet <- as.data.frame(t(MetabolitesDataMet))

# new data frame with only included samples
metabolitesFiltered <- MetabolitesDataMet %>%
  filter(rownames(MetabolitesDataMet) %in% index)

# sanity check
all(rownames(inclusion) %in% rownames(metabolitesFiltered))

# PCA - now filtering the outliers
# scale data without oultiers
scaledMetabolitesFiltered <- (t(metabolitesFiltered) - rowMeans(t(metabolitesFiltered)))/(apply(t(metabolitesFiltered),1,sd))

# check if standard scaling worked
mean(rowMeans(scaledMetabolitesFiltered)) # close to 0
mean(sd(scaledMetabolitesFiltered))       # close to 1

# obtain PCA scores for data without outliers
pcaResultsIT <- pca(t(scaledMetabolitesFiltered), nPcs = 8)
summary(pcaResultsIT)

# obtain metadata based on sample ID from the filtered data
indexMetsIT <- as.list(rownames(pcaResultsIT@scores))
sub_dfIT <- MetaData %>%
  filter(MetaData$ID %in% indexMetsIT)

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
       title = "MMetabolomics data - after outlier removal") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# PCA plot - tumors
plotDataITexluded$Tumor <- as.factor(ifelse(plotDataITexluded$Tumors > 0, "Tumor", "No Tumor"))

# PCA
ggplot(plotDataITexluded, aes(x = PC1, y = PC2, color = Tumor)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsIT@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsIT@R2[2] * 100, 1), "% of the variance"),
       title = "Metabolite PCA plot after outlier removal - Tumor/No Tumor") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# save data with removed outliers 
write.csv(t(metabolitesFiltered), file = "FilteredMetabolite.csv")