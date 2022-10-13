# Data pre-rpocessing and model contruction
# Scientific Programming
# Rosan Olsman

# IMPORTANT !! DO NOT PRE_PROCESS TRAIN AND TEST TOGETHER!!!!!!!!!??


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

rfProfile2 <- rfe(x = data, 
                 y = data[,"Tumor"], 
                 sizes = subsets,
                 rfeControl = rfeCtrl)

rfProfile2



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

## DO PCA ON TUMORS

