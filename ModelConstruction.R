# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman


########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

# load data
PPMicrobiomeData <- read.csv("FilteredMicrobiome.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
PPMetabolite <- read.csv("FilteredMetabolite.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
MetaData <- data.frame(read_excel("./Dataset/MetaTumourData.xlsx"), row.names = "Mouse.ID")

# further pre-process data - microbiome
transMicrobiome <- as.data.frame(asinh(PPMicrobiomeData))

# further pre-process data - metabolite
source("PQNfunction.R")
pqnMetabolites <- pqn(as.matrix(PPMetabolite))
transMetabolites <- log2(pqnMetabolites)
scaleMetabolites <- as.data.frame((transMetabolites - rowMeans(transMetabolites))/(apply(transMetabolites,1,sd)))


# Packages
# Required CRAN packages:
CRANpackages <- c("tidyverse", "readxl", "devtools", "ggplot2", "dplyr", "caret")

# Required Bioconductor packages:
BiocPackages <- c("vioplot", "plotly")

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


#concatenate
ConcatenatePP <- as.data.frame(colnames(scaleMetabolites))
transMicrobiome <- as.data.frame(t(transMicrobiome))
ConcatenatePP <- transMicrobiome %>%
  filter(rownames(transMicrobiome) %in% ConcatenatePP[,1])
remove <- setdiff(colnames(scaleMetabolites), rownames(ConcatenatePP))
scaleMetabolites <- scaleMetabolites[ , !(names(scaleMetabolites) %in% remove)]

# sanity check 
setdiff(colnames(scaleMetabolites), rownames(ConcatenatePP))

# concatenate data
Concatenate <- rbind(scaleMetabolites, t(ConcatenatePP))

ConcatenateMeta <- MetaData %>%
  filter(rownames(MetaData) %in% colnames(Concatenate))

ConcatenateTogether <- as.matrix(rbind(Concatenate, ConcatenateMeta$Category))


set.seed(1234) # !!!!!! DO NOT UNDERSTAND DOES NOT WOK????
trainIndex <- createDataPartition(ConcatenateMeta$Category,
                                  p = 0.80, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)

Train <- as.data.frame(t(Concatenate[,trainIndex]))
Test  <- as.data.frame(t(Concatenate[,-trainIndex]))

TrainMeta <- MetaData %>%
  filter(rownames(MetaData) %in% rownames(Train))


RF = 'cforest'

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)



set.seed(85)
gbmFit1 <- train(Train, TrainMeta$Category,
                 method = "gbm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE)
# yes
gbmFit1
