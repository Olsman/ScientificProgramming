# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman


#### --------------- WORK IN PROGRESS --------------- ####
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#



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
# inverse hyperbolic sine transformation
df_micro <- as.data.frame(asinh(df_micro))

# further pre-process data - metabolite
# probabilistic quotient normalization
source("PQNfunction.R")
df_meta <- pqn(as.matrix(df_meta))
# log2 transformation
df_meta <- log2(df_meta)
# scaling
df_meta <- as.data.frame((df_meta - rowMeans(df_meta))/(apply(df_meta,1,sd)))

# removing samples that are not present in both datasets
# metabolites v microbiome
remove <- setdiff(colnames(df_meta), colnames(df_micro)) # only one sample
df_meta <- df_meta[, !(names(df_meta) %in% remove)]
# microbiome v metabolites
remove <- setdiff(colnames(df_micro), colnames(df_meta))  # 12 samples
df_micro <- df_micro[, !(names(df_micro) %in% remove)]

# remove low feature variances
install.packages('matrixStats')
library(matrixStats)
df_meta = df_meta[rowVars(as.matrix(df_meta)) > 0.01,]
df_micro = df_micro[rowVars(as.matrix(df_micro)) > 0.01,]


# concatenate dataframes - 28 samples remain to train the model
df_com_f <- rbind(df_micro, df_meta)

# install.packages("rpart.plot")
# library(rpart.plot)

# store data and metadata for concatenated datatypes
data <- as.data.frame(t(df_com_f))
data$ID <- rownames(data)
metcat <- metadata %>%
  filter(rownames(metadata) %in% rownames(data))
metcat$ID <- rownames(metcat)

# define new column containing binary classification of tumor presence
metcat$TumorB <- as.factor(ifelse(metcat$Tumors > 0, "1", "0"))
data$Tumor  = metcat$TumorB[match(data$ID,metcat$ID)]

# examine the class imbalance after outlier detection - tumor presence
sum(metcat$TumorB == 1)
sum(metcat$TumorB == 0)

# examine the class imbalance after outlier detection 
ggplot(data = metcat, aes(x = Category, fill = Sex)) +
  geom_bar(position = position_dodge()) +
  theme_classic() +
  labs(title = "Gender per Condition - Concatenated Data", x = "Condition", y = "Count") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# 
data$ID <- NULL

# # DIABLO
# BiocManager::install('mixOmics')
# library(mixOmics)
# X <- list(microbiome = data[,1:10], metabolite = data[,11:425])
# Y <- data[,426]
# 
# result.diablo.tcga <- block.plsda(X, Y)
# plotIndiv(result.diablo.tcga)
# plotVar(result.diablo.tcga,
#         var.names = FALSE)
# 
# plotDiablo(result.diablo.tcga, ncomp = 1)
# cimDiablo(result.diablo.tcga, margin=c(8,20))
# 
# plotLoadings(result.diablo.tcga, comp = 2, contrib = "max")


########### TRAIN MODEL ##########

# divide train and test set
set.seed(12)
train <- createDataPartition(data[,"Tumor"],
                             p = 0.8,
                             list = FALSE)
data.trn <- data[train,]
data.trn <- select(data.trn, -ID)
data.tst <- data[-train,]

# cross validation
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3)

# random forest
fit.cv <- train(Tumor ~ .,
                data = data.trn,
                methods = "rf",
                trControl = ctrl,
                tuneLength = 50) 

# examine output
print(fit.cv)
plot(fit.cv)
fit.cv$results

# examine accuracy of the RF model
pred <- predict(fit.cv, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred))

# select most important features
VarImp <- varImp(fit.cv)
VarImp <- arrange(VarImp$importance, Overall)
featuresRF <- tail(VarImp, n = 20)


########### SVM ###########3
# linear kernel - try others 
fit.cvSVM <- train(Tumor ~ .,
                data = data.trn,
                methods = "svmLinearWeights2",
                trControl = ctrl,
                tuneLength = 50) # only 44? maybe only 44 features unique among samples

# examine output
print(fit.cvSVM)
plot(fit.cvSVM)
fit.cvSVM$results

# examine accuracy of the SVM model
pred2 <- predict(fit.cvSVM, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred2))

# select important features
VarImp2 <- varImp(fit.cvSVM)
VarImp2 <- arrange(VarImp2$importance, Overall)
featuresSVM <- tail(VarImp2, n = 20)


### multinom - can it actually be used for 2 classes; look it up pls
fit.cvMNOM <- train(Tumor ~ .,
             data = data.trn,
             method = "multinom",
             trControl = ctrl)

# examine output
print(fit.cvMNOM)
plot(fit.cvMNOM)
fit.cvMNOM$results

# examine accuracy of the multinom model
pred3 <- predict(fit.cvMNOM, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred3))

# select important features
VarImp3 <- varImp(fit.cvMNOM)
VarImp3 <- arrange(VarImp3$importance, Overall)
featuresMNOM <- tail(VarImp3, n = 20)


############ compare accuracy ############ 

resamps <- resamples(list(RF = fit.cv,
                          SVM = fit.cvSVM,
                          MNOM = fit.cvMNOM))

# compare RF and SVM and MNOM
summary(resamps)


# getModelInfo()$mnom$parameters


##### recursive feature elemination 
# i will probably not do this anymore

# set.seed(123)
# rfeCtrl <- rfeControl(functions = rfFuncs,
#                       method = "cv",
#                       verbose = FALSE)
# 
# # proportion of subsets
# set.seed(123)
# subsets <- c(10, 20, 30, 40, 50, 60)
# 
# drop <- c("Tumor")
# data.trn2 <- data.trn[,!(names(data.trn) %in% drop)]
# 
# rfProfile2 <- rfe(x = data.trn2, 
#                  y = data.trn$Tumor, 
#                  sizes = subsets,
#                  rfeControl = rfeCtrl)
# 
# rfProfile2
# 
# best20features <- predictors(rfProfile2)

# now select these features to train the new model? RF?




##################
# with response as a integer (0/1)
# fit_logistic <- train(Tumor ~.,
#                       data = data.trn,
#                       method = "glmnet",
#                       trControl = ctrl,
#                       family = "binomial")
# print(fit_logistic)
# pred4 <- predict(fit_logistic, data.tst)
# confusionMatrix(table(data.tst[,"Tumor"], pred4))
# 
# VarImp4 <- varImp(fit_logistic)
# VarImp4 <- VarImp4$importance
# 
# VarImp4 <- arrange(VarImp4, Overall)

# featureElemination <- as.data.frame(best20features)
# top20logistic$features <- gsub("`","",rownames(top20logistic))
# rownames(top20rf) <- top20rf$`rownames(top10)`
# top20rf$features <- gsub("`","",rownames(top20rf))



# TO DO:
# - PCA concatenated data
# - assess different model
# - optimize model parameters
# - model fusion; assess mutual and complementary aspects of each model
# - assess feature stability; Jaccard index(?)


# diablo (?)


X <- list(metabolites = df_meta, microbiome = df_micro)
Y <- data.frame(tumor = as.factor(ifelse(metadata$Tumors > 1, "1", "0")),
                sample = rownames(metadata))
Y = Y$sample[match(names(df_micro),Y$sample)]

setdiff(rownames(X[["microbiome"]], 

result.diablo.tcga <- block.plsda(X, Y)          
          
          
          
          