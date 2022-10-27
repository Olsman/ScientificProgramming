# Scientific Programming
# Rosan Olsman

# Explorative Data Analysis
# 10-10-2022

########## Import and examine data ########## 
# go to your personal working directory - set working directory
DIR <- setwd("/Users/rosanolsmanx/Documents/Maastricht University/Courses/MSB1015 Scientific Programming")

# Packages
# Required CRAN packages:
CRANpackages <- c("tidyverse", "readxl", "ggplot2", "dplyr")
# Package version: tidyverse 1.3.2; readxl 1.4.1; ggplot2 3.3.6; dplyr 1.0.10

# Required Bioconductor packages:
BiocPackages <- c("vioplot", "plotly")
# Package version: vioplot 0.3.7; plotly 4.10.0

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
# Go to you own working directory, i.e., where you stored the data
MetaData <- data.frame(read_excel("./Dataset/MetaTumourData.xlsx"), row.names = "Mouse.ID")
MetabolitesData <- data.frame(read_excel("./Dataset/StoolMetabolites.xlsx"), check.names = FALSE)
MicrobiomeData <- data.frame(read_excel("./Dataset/OTUTable.xlsx"), row.names = "ID", check.names = FALSE)

# Sanity check - identify common and different samples for each dataset
all(colnames(MicrobiomeData) %in% rownames(MetaData)) # FALSE: one sample with wrong label - 243 is missing
setdiff(rownames(MetaData), colnames(MicrobiomeData))
all(colnames(MetabolitesData[,14:43]) %in% rownames(MetaData))
setdiff(rownames(MetaData), colnames(MetabolitesData[,14:43]))

# most likely that 233 should have been 243; change 233 to 243
setdiff(colnames(MicrobiomeData), rownames(MetaData))
names(MicrobiomeData)[names(MicrobiomeData) == '233'] <- '243'
all(colnames(MicrobiomeData) %in% rownames(MetaData))

# check for NAs in dataset - no missing values in microbiome and metabolite data
which(colSums(is.na(MetaData))>0)
which(colSums(is.na(MetabolitesData))>0)
which(colSums(is.na(MicrobiomeData))>0)


########## Basic data exploration - metabolite ########## 
# Data visualization - weight
Weight <- ggplot(MetaData, aes(x = Category, y = Weight, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Condition",
       y = "Body Weight (g)",
       title = "Body Weight per Condition") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Weight
# two possible outliers: 151 (WT) & 13 (WT)

# Data visualization - fat mass
FatMass <- ggplot(MetaData, aes(x = Category, y = Fat.mass, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Condition",
       y = "Fat Mass (g)",
       title = "Fat Mass per Condition") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
FatMass
# no clear outliers

# Data visualization - lean mass
LeanMass <- ggplot(MetaData, aes(x = Category, y = Lean.mass, fill = Category)) +
  geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Lean Mass (g)",
       title = "Violin Plot: Lean Mass") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
LeanMass
# two possible outliers: 108 (HF) & 115 (LF)

# Data visualization - sex
Gender <- ggplot(data = MetaData, aes(x = Category, fill = Sex)) +
  geom_bar(position = position_dodge()) +
  theme_classic() +
  labs(title = "Gender per Condition", x = "Condition", y = "Count") +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Gender
# boxplot shows clear class imbalance for sex

# Data visualization - tumor volume
TumorVolume <- ggplot(MetaData, aes(x = Category, y = Tumor.volume, fill = Category)) +
  #geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Condition",
       y = "Tumor Volume",
       title = "Tumor Volume per Condition") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
TumorVolume

# Data visualization - number of tumors
Tumors <- ggplot(MetaData, aes(x = Category, y = Tumors, fill = Category)) +
  #geom_violin(width=1, trim = FALSE) +
  geom_boxplot(width=0.1) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  labs(x = "Category",
       y = "Count",
       title = "No. Tumors per Condition") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
Tumors
