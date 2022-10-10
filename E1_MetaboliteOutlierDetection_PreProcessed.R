########## Isolation forest for outlier removal (with pre-processing metabolomics data) ########## 
# normalize
source("PQNfunction.R")
# metabolites - run PQN function
pqnMetabolites <- pqn(as.matrix(MetabolitesData[,14:43]))

# transform
transMetabolites <- log2(pqnMetabolites)

# scale
# scale data 
scaledMetabolitesPP <- (transMetabolites - rowMeans(transMetabolites))/(apply(transMetabolites,1,sd))
mean(rowMeans(scaledMetabolitesPP)) # close to 0
mean(sd(scaledMetabolitesPP))       # close to 1

# PCA
pcaResultsPP <- pca(t(scaledMetabolitesPP))
summary(pcaResultsPP)

# index - select metadata from samples included in metabolomics data 
indexMetsPP <- as.list(rownames(pcaResultsPP@scores))
sub_dfPP <- MetaData %>%
  filter(rownames(MetaData) %in% indexMetsPP)

# sanity check
all(rownames(pcaResultsPP@scores) %in% rownames(sub_dfPP))

# data for PCA plots
plotDataPP <- merge(data.frame(pcaResultsPP@scores), sub_dfPP, by="row.names", all.x=TRUE)
rownames(plotDataPP) <- plotDataPP$Row.names
plotDataPP$Row.names <- NULL

# PCA plot
ggplot(plotDataPP, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA plot: scaled metabolomics data") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# two possible outliers


########## Outlier removal - isolatin tree with pre-processed metabolomics data ##########
# what happens if i do it on the data
metabolitesPP <- as.data.frame(t(scaledMetabolitesPP))
colnames(metabolitesPP) <- MetabolitesData$COMP_ID
IsolationTreePP <- IsolationTrees(metabolitesPP, rFactor = 0) #rFactor = 0; meaning fully deterministic
AnomalyScorePP <- AnomalyScore(metabolitesPP, IsolationTreePP)
plotDataITPP <- cbind(data.frame(plotDataPP, as.data.frame(AnomalyScorePP[["outF"]])))
colnames(plotDataITPP)[which(names(plotDataITPP) == rev(names(plotDataITPP))[1])] <- "AnomalyScore"
plotDataITPP$Outlier <- as.factor(ifelse(plotDataITPP$AnomalyScore >=0.8, "Outlier", "Normal"))
sum(plotDataITPP$Outlier == "Outlier")

# pca using anomaly score
ggplot(plotDataITPP, aes(x = PC1, y = PC2, color = AnomalyScore, shape = Category)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(5)) +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA: Metabolites coloured by Anomaly Score") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# pca when setting anomaly score to > 0.8
ggplot(plotDataITPP, aes(x = PC1, y = PC2, color = Outlier)) + 
  geom_point() +
  scale_color_manual(values = c("Normal" = "#00BFC4", "Outlier" = "#F8766D")) +
  labs(x=paste("PC1: ", round(pcaResultsPP@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsPP@R2[2] * 100, 1), "% of the variance"),
       title = "PCA: Metabolites coloured by Anonamly Score Treshold") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
# only one outlier visible: sample 49

# select normal data
inclusionPP <- plotDataITPP %>%
  filter(Outlier != "Outlier")

indexPP <- row.names(inclusionPP)
metabolitesFilteredPP <- MetabolitesDataMet %>%
  filter(rownames(MetabolitesDataMet) %in% indexPP)
all(rownames(inclusion) %in% rownames(metabolitesFiltered))
# remove sample 4 and 7 because they appear to be far away from the rest
#metabolitesFilteredPP <- metabolitesFilteredPP[!(rownames(metabolitesFilteredPP) == "4" | rownames(metabolitesFilteredPP) == "7" ),]

#metabolitesFilteredPP <- as.data.frame(scale(metabolitesFiltered), scale = TRUE)

## HERE SOMETHING GOES WRONG I THINK BECAUSE I GET THE SAME PCA

# PCA - now filtering the outliers


# norm??
# trans??
# metabolites - run PQN function
pqnMetabolitesPP <- pqn(as.matrix(metabolitesFilteredPP))

# transform
transMetabolites <- log2(pqnMetabolitesPP)


#
scaledMetabolitesFiltered <- (t(transMetabolites) - rowMeans(t(transMetabolites)))/(apply(t(transMetabolites),1,sd))
mean(rowMeans(scaledMetabolitesFiltered)) # close to 0
mean(sd(scaledMetabolitesFiltered))       # close to 1

pcaResultsIT <- pca(t(scaledMetabolitesFiltered))
summary(pcaResultsIT)

# index 
indexMetsIT <- as.list(rownames(pcaResultsIT@scores))
sub_dfIT <- MetaData %>%
  filter(rownames(MetaData) %in% indexMetsIT)

# sanity check
all(rownames(pcaResultsIT@scores) %in% rownames(sub_dfIT))

# data for PCA plots
plotDataITexluded <- merge(data.frame(pcaResultsIT@scores), sub_dfIT, by="row.names", all.x=TRUE)
rownames(plotDataITexluded) <- plotDataITexluded$Row.names
plotDataITexluded$Row.names <- NULL

# PCA plot
ggplot(plotDataITexluded, aes(x = PC1, y = PC2, color = Category, shape = Sex)) +
  geom_point() +
  labs(x=paste("PC1: ", round(pcaResultsIT@R2[1] * 100, 1), "% of the variance"),
       y=paste("PC2: ", round(pcaResultsIT@R2[2] * 100, 1), "% of the variance"),
       title = "PCA plot: after outlier removal by Isolation Tree") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# IT IS CLEAR: ONLY SAMPLE 49 IS AN OUTLIER - PCA before normalization and transformation appears to be better

########## DATA ANALYSIS ##########
normalizedMetabolitest <- pqn(t(metabolitesFilteredPP))
transformationMetabolites <- log2(normalizedMetabolitest)
scaledMetabolitesAnalysis <- (t(transformationMetabolites) - rowMeans(t(transformationMetabolites)))/(apply(t(transformationMetabolites),1,sd))

