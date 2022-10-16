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

# concatenate dataframes - 28 samples remain to train the model
df_com_f <- rbind(df_micro, df_meta)

install.packages("rpart.plot")
library(rpart.plot)

# store data and metadata for concatenated datatypes
data <- as.data.frame(t(df_com_f))
data$ID <- rownames(data)
metcat <- metadata %>%
  filter(rownames(metadata) %in% rownames(data))
metcat$ID <- rownames(metcat)

# define new column containing binary classification of tumor presence
metcat$TumorB <- as.factor(ifelse(metcat$Tumors > 0, "1", "0"))
data$Tumor <- DF1$Col = DF2$Col[match(DF1$Name,DF2$Name)]


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

########### TRAIN MODEL ##########
# set seed for reproducibility
set.seed(667)
train <- createDataPartition(data[,"Tumor"],
                             p = 0.8,
                             list = FALSE)

data.trn <- data[train,]
data.tst <- data[-train,]

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 3)

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



VarImp <- arrange(VarImp, Overall)
top10 <- tail(VarImp, n = 20)
###########

fit.cvSVM <- train(Tumor ~ .,
                data = data.trn,
                methods = "svmLinearWeights2",
                trControl = ctrl,
                tuneLength = 50) # only 44? maybe only 44 features unique among samples

print(fit.cvSVM)
plot(fit.cvSVM)
fit.cvSVM$results

pred2 <- predict(fit.cvSVM, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred2))

VarImp2 <- varImp(fit.cvSVM)
VarImp2 <- VarImp2$importance


### multinom
mnom <- train(Tumor ~ .,
             data = data.trn,
             method = "multinom",
             trControl = ctrl)
print(mnom)
plot(mnom)
mnom$results

pred3 <- predict(mnom, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred3))

VarImp3 <- varImp(mnom)
VarImp3 <- VarImp3$importance


############################

resamps <- resamples(list(RF = fit.cv,
                          SVM = fit.cvSVM,
                          MNOM = mnom))

# compare RF and SVM
summary(resamps)


getModelInfo()$mnom$parameters


##### recursive feature elemination
set.seed(355)

rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "cv",
                      verbose = FALSE)

# proportion of subsets
set.seed(355)
subsets <- c(10, 20, 30, 40, 50, 60)

drop <- c("Tumor")
data.trn2 <- data.trn[,!(names(data.trn) %in% drop)]

rfProfile2 <- rfe(x = data.trn2, 
                 y = data.trn$Tumor, 
                 sizes = subsets,
                 rfeControl = rfeCtrl)

rfProfile2

best20features <- predictors(rfProfile2)

# now select these features to train the new model? RF?


# examine the class imbalance after outlier detection 
ggplot(data = metcat, aes(x = TumorB)) +
  geom_bar(position = position_dodge()) +
  theme_classic() +
  labs(title = "Number of mice with or without tumour", x = "Condition", y = "Count") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

sum(metcat$TumorB == 1)
sum(metcat$TumorB == 0)

##################
# with response as a integer (0/1)
fit_logistic <- train(Tumor ~.,
                      data = data.trn,
                      method = "glmnet",
                      trControl = ctrl,
                      family = "binomial")
print(fit_logistic)
pred4 <- predict(fit_logistic, data.tst)
confusionMatrix(table(data.tst[,"Tumor"], pred4))

VarImp4 <- varImp(fit_logistic)
VarImp4 <- VarImp4$importance

VarImp4 <- arrange(VarImp4, Overall)
top20logistic <- tail(VarImp4, n = 20)

featureElemination <- as.data.frame(best20features)
top20logistic$features <- gsub("`","",rownames(top20logistic))
rownames(top20rf) <- top20rf$`rownames(top10)`
top20rf$features <- gsub("`","",rownames(top20rf))




