# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman


########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

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

# Loading data
df_micro <- read.csv("FilteredMicrobiome.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE)  # Microbiome rows:columns -> taxa:samples
df_meta <- read.csv("FilteredMetabolite.csv", header = TRUE, row.names = 1, sep = ",", check.names = FALSE) # Metabolomics rows:columns -> metabolites:samples
metadata <- data.frame(read_excel("./Dataset/MetaTumourData.xlsx"), row.names = "Mouse.ID")

# further pre-process data - microbiome
df_micro <- as.data.frame(asinh(df_micro))

# further pre-process data - metabolite
source("PQNfunction.R")
df_meta <- pqn(as.matrix(df_meta))
df_meta <- log2(df_meta)
df_meta <- as.data.frame((df_meta - rowMeans(df_meta))/(apply(df_meta,1,sd)))

remove <- setdiff(colnames(df_meta), colnames(df_micro))
df_meta <- df_meta[, !(names(df_meta) %in% remove)]

remove <- setdiff(colnames(df_micro), colnames(df_meta))
df_micro <- df_micro[, !(names(df_micro) %in% remove)]

df_com_f <- rbind(df_micro, df_meta)

install.packages("rpart.plot")
library(rpart.plot)

str(t(df_com_f))
summary(df_com_f)

data <- as.data.frame(t(df_com_f))
metcat <- metadata %>%
  filter(rownames(metadata) %in% rownames(data))
metcat$TumorB <- as.factor(ifelse(metcat$Tumors > 0, "1", "0"))
data$Tumor <- metcat$TumorB
# data$ID <- rownames(metcat)

set.seed(667)
train <- createDataPartition(data[,"Tumor"],
                             p = 0.8,
                             list = FALSE)

data.trn <- data[train,]
data.tst <- data[-train,]

ctrl <- trainControl(method = "cv",
                     number = 10)

fit.cv <- train(Tumor ~ .,
                data = data.trn,
                methods = "rf",
                trControl = ctrl,
                tuneLength = 50) # only 44?


print(fit.cv)
plot(fit.cv)
fit.cv$results

pred <- predict(fit.cv, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred))

VarImp <- varImp(fit.cv)
VarImp <- VarImp$importance

ggplot(data = metcat, aes(x = Category, fill = Sex)) +
  geom_bar(position = position_dodge()) +
  theme_classic() +
  labs(title = "Gender per Condition", x = "Condition", y = "Count") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

VarImp <- arrange(VarImp, Overall)
top10 <- tail(VarImp, n = 20)
###########
























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
