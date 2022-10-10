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

# set in right format
NumMicrobiome <- as.data.frame(sapply(MicrobiomeData, as.numeric))
row.names(NumMicrobiome) <- rownames(MicrobiomeData)
NumMicrobiome <- as.matrix(NumMicrobiome)

# transMicrobiome <- asinh(NumMicrobiome)

m_com = as.matrix(NumMicrobiome)
set.seed(123)
nmds = metaMDS(t(m_com), distance = "bray")
nmds
plot(nmds)

# download vegan here, otherwise you will get an error
install.packages("vegan")
library(vegan)

data.scores = as.data.frame(scores(nmds)$sites)
all(rownames(data.scores) %in% rownames(t(MicrobiomeData)))
setdiff(rownames(t(MicrobiomeData)), rownames(MetaData))

indexMB <- as.data.frame(rownames(data.scores))
sub_dfMB <- MetaData %>%
  filter(rownames(MetaData) %in% indexMB[,1]) 

#
data.scores$Category = sub_dfMB$Category
head(data.scores)

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
  labs(x = "NMDS1", colour = "Category", y = "NMDS2") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")


# isolation forest on NMDS score?
microbiome <- as.data.frame(t(NumMicrobiome))
IsolationTreeMicrobiome <- IsolationTrees(microbiome, rFactor = 0)
AnomalyScoreMicrobiome <- AnomalyScore(microbiome, IsolationTreeMicrobiome)
plotDataMicrobiome <- cbind(data.frame(data.scores, as.data.frame(AnomalyScoreMicrobiome[["outF"]])))
colnames(plotDataMicrobiome)[which(names(plotDataMicrobiome) == rev(names(plotDataMicrobiome))[1])] <- "AnomalyScore"
plotDataMicrobiome$Outlier <- as.factor(ifelse(plotDataMicrobiome$AnomalyScore >=0.8, "Outlier", "Normal"))
sum(plotDataMicrobiome$Outlier == "Outlier")

# NMDA PLOT USING THE ANOMALY SCORE 
# pca using anomaly score
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
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", main = "NMDA PLOT: OUTLIER DETECTION") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")

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
  labs(x = "NMDS1", colour = "Category", y = "NMDS2", title = "NMDA PLOT: OUTLIER DETECTION - ANOMALY") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")


########## REMOVE OUTLIER FROM ISOLATION FOREST - MICROBIOME DATA ##########
# select normal data
inclusionMicrobiome <- plotDataMicrobiome %>%
  filter(Outlier != "Outlier")

indexMicrobiome <- row.names(inclusionMicrobiome)
microbiomeFiltered <- as.data.frame(t(MicrobiomeData)) %>%
  filter(rownames(t(MicrobiomeData)) %in% indexMicrobiome)
all(rownames(inclusionMicrobiome) %in% rownames(microbiomeFiltered))

# data to numeric
NumMicrobiomeFiltered <- as.data.frame(sapply(microbiomeFiltered, as.numeric))
row.names(NumMicrobiomeFiltered) <- rownames(inclusionMicrobiome)


# NMDA plot without outliers
m_comFiltered = as.matrix(NumMicrobiomeFiltered)
set.seed(123)
nmdsFiltered = metaMDS(m_comFiltered, distance = "bray")
nmdsFiltered
plot(nmdsFiltered)

data.scoresFiltered = as.data.frame(scores(nmdsFiltered)$sites)
all(rownames(data.scoresFiltered) %in% rownames(t(MicrobiomeData)))
setdiff(rownames(t(MicrobiomeData)), rownames(MetaData))

indexMBFiltered <- as.data.frame(rownames(data.scoresFiltered))
sub_dfMBFiltered <- MetaData %>%
  filter(rownames(MetaData) %in% indexMBFiltered[,1]) 

#
data.scoresFiltered$Category = sub_dfMBFiltered$Category
head(data.scoresFiltered)

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
  labs(x = "NMDS1", colour = "Category", y = "NMDS2") +
  annotate("text", x = -0.6, y = 0.5, label = "Stress = 0.09")

