# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman


########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

# Packages
# Required CRAN packages:
CRANpackages <- c("tidyverse", "readxl", "devtools", "ggplot2", "dplyr", "caret"
                  , "matrixStats", "MLeval", "randomForest", "cluster", "pheatmap")
# Package version: tidyverse 1.3.2; readxl 1.4.1; devtools 2.4.5; ggplot2 3.3.6;
# Package version: dplyr 1.0.10; caret 6.0-93; matrixStats 0.62.0; MLeval 0.3;
# Package version: randomForest 4.7-1.1; cluster 2.1.4; pheatmap 1.0.12

# Required Bioconductor packages:
BiocPackages <- c("vioplot", "plotly", "pcaMethods")
# Package version: vioplot 0.3.7; plotly 4.10.0; pcaMethods 1.88.0

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

# examine distribution of the data
par(mar=c(2,2,2,2))
hist(as.matrix(df_meta), breaks = 50)
hist(as.matrix(df_micro), breaks = 50)

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

# examine distribution
par(mar=c(2,2,2,2))
hist(as.matrix(df_meta), breaks = 50)
hist(as.matrix(df_micro), breaks = 50)
# looks good for metabolomics data, microbiome not so much

# removing samples that are not present in both datasets
# metabolites v microbiome
remove1 <- setdiff(colnames(df_meta), colnames(df_micro)) # only one sample
df_meta <- df_meta[, !(names(df_meta) %in% remove1)]
# microbiome v metabolites
remove2 <- setdiff(colnames(df_micro), colnames(df_meta))  # 12 samples
df_micro <- df_micro[, !(names(df_micro) %in% remove2)]

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
#####---------- Random Forest ----------#####
# get right names otherwise caret cannot work with it
data$Tumor <- make.names(data$Tumor, unique = FALSE, allow_ = TRUE)

# grid search for mtry
tuneGrid <- expand.grid(.mtry = c(1 : 10))

# cross validation
set.seed(12)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 5,
                     search = 'grid',
                     classProbs = TRUE,
                     savePredictions = "final",
                     summaryFunction = twoClassSummary) #in most cases a better summary for two class problems 

# also ook at ntrees and nodesize
ntrees <- c(500, 1000)    
nodesize <- c(1, 5)

# save hyperparameter values in params
params <- expand.grid(ntrees = ntrees,
                      nodesize = nodesize)

# to save output
store_maxnode <- vector("list", nrow(params))

# RF forest with hyperparameter tuning
for(i in 1:nrow(params)){
  nodesize <- params[i,2]
  ntree <- params[i,1]
  set.seed(65)
  rf_model <- train(Tumor~.,
                    data = data,
                    method = "rf",
                    importance=TRUE,
                    metric = "ROC",
                    tuneGrid = tuneGrid,
                    trControl = ctrl,
                    ntree = ntree,
                    nodesize = nodesize)
  store_maxnode[[i]] <- rf_model # save output of each iteration
}

# save output
names(store_maxnode) <- paste("ntrees:", params$ntrees,
                              "nodesize:", params$nodesize)
results_mtry <- resamples(store_maxnode)

# examine RF output
summary(results_mtry)
lapply(store_maxnode, function(x) x$best)
lapply(store_maxnode, function(x) x$results[x$results$ROC == max(x$results$ROC),])
# first model, i.e., ntrees = 500, nodesize = 1, and mtry = 5 gives an ROC of 0.715


#####---------- Elastic net ----------#####
set.seed(12)
model_en <- train(Tumor~.,
            data = data,
            method = "glmnet",
            metric = "ROC",
            tuneGrid = expand.grid(alpha = seq(0,1,length=10),
                                  lambda = seq(0.0001,0.2,length=5)),
            trControl= ctrl)
# examine model output for varying values of both alpha and lambda
model_en
mean(model_en$resample$ROC)

# plot model output
plot(model_en, main = "Elastic Net Regression")
# best model is obtained from alpha = 0.22222 and lambda = 0.10005
# ROC = 0.715


#####---------- Unsupervised RF ----------#####
dat <- data[,-ncol(data)]
mtry <- c(1,5,10)
paramsURF <- expand.grid(ntrees = ntrees,
                         nodesize = nodesize,
                         mtry = mtry)

store_maxnodeURF <- vector("list", nrow(paramsURF))
store_proximities <- vector("list", nrow(paramsURF))

# unsupervised random forest
for(i in 1:nrow(paramsURF)){
  nodesize <- paramsURF[i,2]
  ntree <- paramsURF[i,1]
  mtry <- paramsURF[i,3]
  set.seed(12)
  urf_model <- randomForest(x = dat,
                            mtry = mtry,
                            ntree = ntree,
                            nodesize = nodesize,
                            proximity = TRUE)
  store_maxnodeURF[[i]] <- urf_model
  store_proximities[[i]] <- urf_model$proximity
}

# save proximity in new variable
for(i in 1:nrow(paramsURF)){
  
  # save proximitiy matrix
  M = store_proximities[[i]]
  
  # compute the row-wise and column-wise mean matrices
  R = M*0 + rowMeans(M)
  C = t(M*0 + colMeans(M))  
  
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
  plots = list()
  
  plots[[i]] = ggplot(plotData, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
    geom_point() +
    labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
         y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
         title = paste("Unsupervised Random Forest PCA - Category", rownames(paramsURF[i,]))) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
  print(plots[[i]])
  
  # get binary column for tumor presence
  plotData$TumorB = ifelse(plotData$Tumors > 0, "1", "0")
  
  plots2 <- list()
  
  # PCA plot - URF
  plots2[[i]] <- ggplot(plotData, aes(x = PC1, y = PC2, color = TumorB, shape = Sex)) +
    geom_point() +
    labs(x=paste("PC1: ", round(pcaResults@R2[1] * 100, 1), "% of the variance"),
         y=paste("PC2: ", round(pcaResults@R2[2] * 100, 1), "% of the variance"),
         title = paste("Unsupervised Random Forest PCA - Tumor Presence", rownames(paramsURF[i,]))) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
  print(plots2[[i]])
}

# proximity 4 best: 1000 ntrees, 1 nodesize, and mtry 5.

# similarity matrix based on condition - proximity 4
proxAnn = store_proximities[[4]]
rownames(proxAnn) = metadata$Category[match(rownames(proxAnn), rownames(metadata))]
colnames(proxAnn) = metadata$Category[match(colnames(proxAnn), rownames(metadata))]
proxAnn = 1-proxAnn
pheatmap::pheatmap(proxAnn,
                   main = "Similarity URF - Condition")
# no clustering visible

# get new variable yes if there is a tumor and no if there is not a tumor
metadata$tumorB = ifelse(metadata$Tumors > 0, "Yes", "No")

# similarity matrix based on tumor presence
proxAnnT = store_proximities[[4]]
rownames(proxAnnT) = metadata$tumorB[match(rownames(proxAnnT), rownames(metadata))]
colnames(proxAnnT) = metadata$tumorB[match(colnames(proxAnnT), rownames(metadata))]
proxAnnT = 1-proxAnnT
pheatmap::pheatmap(proxAnnT,
                   main = "Similarity Matrix URF - Tumor Presence")
# no clustering visible


s#########-----------DIABLO/FACTORIZATION - DO NOT RUN----------------############
# skip this part as i was unable to really use it #
BiocManager::install('mixOmics')
library(mixOmics)
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