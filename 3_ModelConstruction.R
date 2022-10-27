# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman


########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

# Packages
# Required CRAN packages:
CRANpackages <- c("tidyverse", "readxl", "devtools", "ggplot2", "dplyr", "caret"
                  , "matrixStats", "MLeval")

# Required Bioconductor packages:
BiocPackages <- c("vioplot", "plotly", "pcaMethods")

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
df_micro = df_micro[rowVars(as.matrix(df_micro)) > 0.01,]
df_micro <- as.data.frame(asinh(df_micro))

# further pre-process data - metabolite
# probabilistic quotient normalization
df_meta = df_meta[rowVars(as.matrix(df_meta)) > 0.01,]
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

# concatenate dataframes - 28 samples remain to train the model
df_com_f <- rbind(df_micro, df_meta)

# store data and metadata for concatenated datatypes
data <- as.data.frame(t(df_com_f))
data$ID <- rownames(data)
metcat <- metadata %>%
  filter(rownames(metadata) %in% rownames(data))
metcat$ID <- rownames(metcat)

# define new column containing binary classification of tumor presence
metcat$TumorB <- ifelse(metcat$Tumors > 0, "1", "0")
data$Tumor  = metcat$TumorB[match(data$ID,metcat$ID)]

# examine the class imbalance after outlier detection - tumor presence
sum(metcat$TumorB == 1)
sum(metcat$TumorB == 0)

# remove data id column 
data$ID <- NULL


#########-----------DIABLO/FACTORIZATION - DO NOT RUN----------------############
# skip this part as i was unable to really use it #
# BiocManager::install('mixOmics')
# library(mixOmics)
X <- list(microbiome = data[,1:10], metabolite = data[,11:425])
Y <- data[,426]
result.diablo.tcga <- block.plsda(X, Y)
plotIndiv(result.diablo.tcga,
          ind.names = data$category,
          ellipse = T, 
          legend = T)
plotLoadings
plotVar(result.diablo.tcga,
        var.names = FALSE)
plotDiablo(result.diablo.tcga, ncomp = 1)
cimDiablo(result.diablo.tcga, margin=c(8,20))
plotLoadings(result.diablo.tcga, comp = 2, contrib = "max")

# or use factorization 
install.packages("FactoMineR")
r.mfa <- FactoMineR::MFA(
  t(rbind(df_meta,df_micro)), # binding the omics types together
  c(dim(df_meta)[1], dim(df_micro)[1]), # specifying the dimensions of each
  graph=FALSE)

# first, extract the H and W matrices from the MFA run result
mfa.h <- as.data.frame(r.mfa$global.pca$ind$coord)
mfa.w <- r.mfa$quanti.var$coord

# create a dataframe with the H matrix and the CMS label
mfa_df <- as.data.frame(mfa.h)
mfa_df$subtype <- metcat$Category[match(rownames(mfa_df), metcat$ID)]
mfa_df$tumor <- metcat$TumorB[match(rownames(mfa_df), metcat$ID)]

# create the plot
ggplot2::ggplot(mfa_df, ggplot2::aes(x=Dim.1, y=Dim.2, color=tumor)) +
  ggplot2::geom_point() + ggplot2::ggtitle("Scatter plot of MFA")
#########-----------END----------------############


########### TRAIN MODEL ##########
# cross validation
data$Tumor <- make.names(data$Tumor, unique = FALSE, allow_ = TRUE)

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3,
                     classProbs = T,
                     savePredictions = T,
                     summaryFunction = twoClassSummary)

# random forest 
set.seed(12)
fit.cv <- train(Tumor ~ .,
                data = data,
                methods = "rf",
                trControl = ctrl, 
                tuneLength = 50) 
# examine output
MLeval::evalm(fit.cv)

print(fit.cv)
plot(fit.cv)
fit.cv$results

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
          

# 1. Random forest
data$Tumor <- make.names(data$Tumor, unique = FALSE, allow_ = TRUE)

# cross val
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     classProbs = T,
                     savePredictions = T)
rfGrid = expand.grid(mtry = 1:100)
fit.cv <- train(Tumor ~ .,
                data = data,
                methods = "rf",
                trControl = ctrl,
                tuneGrid = rfGrid) 
fit.cv
summary(fit.cv)
x <- MLeval::evalm(fit.cv)


# 2. gradient boosting
gbmGrid <-  expand.grid(interaction.depth = 10,
                        n.trees = 18000,                                          
                        shrinkage = 0.01,                                         
                        n.minobsinnode = 4) 

# Build using a gradient boosted machine
set.seed(5627)
gbm <- train(Tumor ~ .,
             data = data,
             method = "gbm",
             metric = "ROC",
             tuneGrid = gbmGrid,
             verbose = FALSE,
             trControl = ctrl) 
gbm
summary(gbm)
x <- MLeval::evalm(gbm)


# URF
myColRamp <- colorRampPalette(colors = c("#25591f", "#818c3c"))

install.packages("randomForest")
library(randomForest)
install.packages("cluster")
library(cluster)
dat <- data[,-ncol(data)]

rf2 <- randomForest(x = dat, mtry = 25, ntree = 2000, proximity = TRUE)
rf2
prox <- rf2$proximity

# example data
M = prox

# compute the row-wise and column-wise mean matrices
R = M*0 + rowMeans(M)  # or `do.call(cbind, rep(list(rowMeans(tst)), 3))`
C = t(M*0 + colMeans(M))  # or `do.call(rbind, rep(list(colMeans(tst)), 3))`

# substract them and add the grand mean
M_double_centered = M - R - C + mean(M[])         

# do PCA on scaled data
pcaResults <- pca(M_double_centered)
summary(pcaResults)

indexMets <- as.list(rownames(pcaResults@scores))
sub_df <- metadata %>%
  filter(rownames(metadata) %in% indexMets)          

# sanity check
all(rownames(pcaResults@scores) %in% rownames(sub_df))

# data for PCA plots
plotData <- merge(data.frame(pcaResults@scores), sub_df, by="row.names", all.x=TRUE)
rownames(plotData) <- plotData$Row.names
plotData$Row.names <- NULL

# PCA plot - URF
ggplot(plotData, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
       title = "Unsupervised Random Forest PCA - Category") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# get binary column for tumor presence
plotData$TumorB = ifelse(plotData$Tumors > 0, "1", "0")

# PCA plot - URF
ggplot(plotData, aes(x = PC1, y = PC2, color = TumorB, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
       title = "Unsupervised Random Forest PCA - Tumor Presence") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# proximity
proxAnn = prox
rownames(proxAnn) = metadata$Category[match(rownames(proxAnn), rownames(metadata))]
colnames(proxAnn) = metadata$Category[match(colnames(proxAnn), rownames(metadata))]
proxAnn = 1-proxAnn
pheatmap::pheatmap(proxAnn,
                   main = "Proximity Matrix URF - Category")

# or prox based on tumor presence
metadata$tumorB = ifelse(metadata$Tumors > 0, "Yes", "No")

# aa
proxAnnT = prox
rownames(proxAnnT) = metadata$tumorB[match(rownames(proxAnnT), rownames(metadata))]
colnames(proxAnnT) = metadata$tumorB[match(colnames(proxAnnT), rownames(metadata))]
proxAnnT = 1-proxAnnT
pheatmap::pheatmap(proxAnnT,
                   main = "Proximity Matrix URF")
