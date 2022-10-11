# Microbiome outlier detection
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

# Sanity check - identify removed samples
all(colnames(MicrobiomeData) %in% rownames(MetaData)) # FALSE: one sample with wrong label - 243 is missing
setdiff(rownames(MetaData), colnames(MicrobiomeData))
all(colnames(MetabolitesData[,14:43]) %in% rownames(MetaData))
setdiff(rownames(MetaData), colnames(MetabolitesData[,14:43]))

# change 233 to 243
setdiff(colnames(MicrobiomeData), rownames(MetaData))
names(MicrobiomeData)[names(MicrobiomeData) == '233'] <- '243'
all(colnames(MicrobiomeData) %in% rownames(MetaData))

# check for NAs in dataset - no missing values in microbiome and metabolite data
which(colSums(is.na(MetaData))>0)
which(colSums(is.na(MetabolitesData))>0)
which(colSums(is.na(MicrobiomeData))>0)

# set in right format - set scientific notation (e.g., E-04 to numeric)
NumMicrobiome <- as.data.frame(sapply(MicrobiomeData, as.numeric))
row.names(NumMicrobiome) <- rownames(MicrobiomeData)
NumMicrobiome <- as.matrix(NumMicrobiome)

# download vegan here, otherwise you will get an error - conflicting packages
# if you get an error; remove packages, then download packages above before installing vegan
install.packages("vegan")
library(vegan)

# obtain a matrix containing microbiome data set as numeric
m_com = as.matrix(NumMicrobiome)
# set.seed(123)
nmds = metaMDS(t(m_com), distance = "bray")
nmds     # stress: 0.09, which is ok

# obtain the data scores for NMDS plot
data.scores = as.data.frame(scores(nmds)$sites)

# sanity check
all(rownames(data.scores) %in% rownames(t(MicrobiomeData)))
setdiff(rownames(t(MicrobiomeData)), rownames(MetaData))

# obtain index values for samples included in microbiome data analysis
indexMB <- as.data.frame(rownames(data.scores))
sub_dfMB <- MetaData %>%
  filter(rownames(MetaData) %in% indexMB[,1]) 

# add the right category for each sample
data.scores$Category <- sub_dfMB$Category
head(data.scores)

# plot NMDS plot
ggplot(data.scores, aes(x = NMDS1, y = NMDS2, color = Category)) + 
  geom_point() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", title = "Microbiome NMDS plot before pre-processing") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")


# obtain dataframe with numerical microbiome data
microbiome <- as.data.frame(t(NumMicrobiome))

########## Isolation forest for outlier removal (without further pre-processing metabolomics data) ########## 
install.packages("IsolationForest", repos="http://R-Forge.R-project.org")
# yes - packages that need compilation
library("IsolationForest")

# train isolation forest
IsolationTreeMicrobiome <- IsolationTrees(microbiome, rFactor = 0)   #rFactor = 0; fully deterministic
AnomalyScoreMicrobiome <- AnomalyScore(microbiome, IsolationTreeMicrobiome)

# plotting data for NMDS plot including the anomaly score
plotDataMicrobiome <- cbind(data.frame(data.scores, as.data.frame(AnomalyScoreMicrobiome[["outF"]])))
colnames(plotDataMicrobiome)[which(names(plotDataMicrobiome) == rev(names(plotDataMicrobiome))[1])] <- "AnomalyScore"

# new column assigning each sample with an anomaly score ≥ as an outlier
plotDataMicrobiome$Outlier <- as.factor(ifelse(plotDataMicrobiome$AnomalyScore >=0.8, "Outlier", "Normal"))
sum(plotDataMicrobiome$Outlier == "Outlier")

# NMDS plot; samples coloured by the anomaly score
ggplot(plotDataMicrobiome, aes(x = NMDS1, y = NMDS2, color = AnomalyScore, shape = Category)) + 
  geom_point() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", title = "Microbiome NMDS score coloured by anomaly score") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")

# NMDS plot; samples coloured if they are considered to be an outlier
ggplot(plotDataMicrobiome, aes(x = NMDS1, y = NMDS2, color = Outlier)) + 
  geom_point() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", title = "Microbiome NMDS score coloured by outlier") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")


########## REMOVE OUTLIER FROM ISOLATION FOREST - MICROBIOME DATA ##########
# select samples that are not an outlier
inclusionMicrobiome <- plotDataMicrobiome %>%
  filter(Outlier != "Outlier")

# obtain index values and samples that are not an outlier
indexMicrobiome <- row.names(inclusionMicrobiome)
microbiomeFiltered <- as.data.frame(t(MicrobiomeData)) %>%
  filter(rownames(t(MicrobiomeData)) %in% indexMicrobiome)

# sanity check
all(rownames(inclusionMicrobiome) %in% rownames(microbiomeFiltered))

# data to numeric and add the samples
NumMicrobiomeFiltered <- as.data.frame(sapply(microbiomeFiltered, as.numeric))
row.names(NumMicrobiomeFiltered) <- rownames(inclusionMicrobiome)

# NMDS plot without outliers
m_comFiltered = as.matrix(NumMicrobiomeFiltered)
# set.seed(123)
nmdsFiltered = metaMDS(m_comFiltered, distance = "bray")
nmdsFiltered
# plot(nmdsFiltered)

# obtain the data scores for NMDS plot without the outliers
data.scoresFiltered = as.data.frame(scores(nmdsFiltered)$sites)

# sanity check
all(rownames(data.scoresFiltered) %in% rownames(t(MicrobiomeData)))
setdiff(rownames(t(MicrobiomeData)), rownames(MetaData))

# obtain index values for samples included in microbiome data analysis
indexMBFiltered <- as.data.frame(rownames(data.scoresFiltered))
sub_dfMBFiltered <- MetaData %>%
  filter(rownames(MetaData) %in% indexMBFiltered[,1]) 

# add the right category for each sample
data.scoresFiltered$Category <- sub_dfMBFiltered$Category
head(data.scoresFiltered)

# NMDS plot without the outliers
ggplot(data.scoresFiltered, aes(x = NMDS1, y = NMDS2, color = Category)) + 
  geom_point() +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", title = "Microbiome NMDS score without outliers") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")


# save data without the outliers
write.csv(NumMicrobiomeFiltered, file = "FilteredMicrobiome.csv", row.names = TRUE)