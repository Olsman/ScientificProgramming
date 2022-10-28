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

par(mar=c(1,1,1,1))
hist(as.matrix(df_meta))
hist(as.matrix(df_micro))

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

par(mar=c(1,1,1,1))
hist(as.matrix(df_meta))
hist(as.matrix(df_micro))

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


########### TRAIN MODEL ##########
# caret has a hard time dealing with 0 and 1, therefore make new names 
data$Tumor <- make.names(data$Tumor, unique = FALSE, allow_ = TRUE)

# cross validation
fitControl1 <- trainControl(
  method = 'LOOCV',                
  number = 1,                     
  savePredictions = T,        
  classProbs = T ,
  seed = as.list(rep(1,425)),                
  summaryFunction=twoClassSummary 
) 

# model - RF
fit.cv1 <- train(Tumor ~ .,
                data = data,
                methods = "rf",
                trControl = fitControl1,
                tuneGrid = data.frame(mtry=3)) 

# get model output
fit.cv1
summary(fit.cv1)

# plot the ROC curve
ROC = MLeval::evalm(fit.cv1)

# feature importance from fitted model
VarImp <- varImp(fit.cv1)
VarImp <- arrange(VarImp$importance, Overall)
featuresRF <- tail(VarImp, n = 20)


# URF
install.packages("randomForest")
library(randomForest)
install.packages("cluster")
library(cluster)
dat <- data[,-ncol(data)]

# unsupervised random forest
rf2 <- randomForest(x = dat, mtry = 75, ntree = 2000, proximity = TRUE)
rf2
prox <- rf2$proximity

# feature importance URF
VarImpURF <- as.data.frame(rf2[["importance"]])
VarImpURF = VarImpURF %>%                                      
  arrange(desc(MeanDecreaseGini))
featuresURF <- head(VarImpURF, n = 20)

# save proximity in new variable
M = prox

# compute the row-wise and column-wise mean matrices
R = M*0 + rowMeans(M)  # or `do.call(cbind, rep(list(rowMeans(tst)), 3))`
C = t(M*0 + colMeans(M))  # or `do.call(rbind, rep(list(colMeans(tst)), 3))`

# substract them and add the grand mean
M_double_centered = M - R - C + mean(M[])         

# do PCA on scaled data
pcaResults <- pca(M_double_centered)
summary(pcaResults)

# obtain metadata
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

# similarity matrix based on condition
proxAnn = prox
rownames(proxAnn) = metadata$Category[match(rownames(proxAnn), rownames(metadata))]
colnames(proxAnn) = metadata$Category[match(colnames(proxAnn), rownames(metadata))]
proxAnn = 1-proxAnn
pheatmap::pheatmap(proxAnn,
                   main = "Similarity URF - Condition")

# get new variable yes if there is a tumor and no if there is not a tumor
metadata$tumorB = ifelse(metadata$Tumors > 0, "Yes", "No")

# similarity matrix based on tumor presence
proxAnnT = prox
rownames(proxAnnT) = metadata$tumorB[match(rownames(proxAnnT), rownames(metadata))]
colnames(proxAnnT) = metadata$tumorB[match(colnames(proxAnnT), rownames(metadata))]
proxAnnT = 1-proxAnnT
pheatmap::pheatmap(proxAnnT,
                   main = "Similarity Matrix URF - Tumor Presence")


# examine features both present in RF and URF
intersect = intersect(rownames(featuresRF), rownames(featuresURF))
# only trypthophan intersect

s#########-----------DIABLO/FACTORIZATION - DO NOT RUN----------------############
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