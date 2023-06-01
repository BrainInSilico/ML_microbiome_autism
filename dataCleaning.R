library(stringr)
library(randomForest)
library(Boruta)
library(mlbench)
library(caret)

countNbVariant <- function(x) {
 return(length(which(x>0)))
}

isVariant <- function(x, threshold=0.05) {
  bool <- F
  if(countNbVariant(x)/length(x) >= threshold){
    bool <- T
  }
  return(bool)
}

#-------------------------------- main -----------------------------

# import data
path <- 'Microbiome_ASD/GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv'
abundance <- read.delim2(path, sep=',')
otuTable <- read.delim2('Microbiome_ASD/ASD_meta_abundance.csv')

# extract taxonomy, OTU names and targets
taxonomy <- data.frame(otu=abundance$OTU, taxonomy=abundance$taxonomy)
abundance <- abundance[,3:ncol(abundance)]
rownames(abundance) <- taxonomy$otu
target <- str_extract(colnames(abundance), '^[A,B]')

# Compute relative abundance
relativeAbundance <- apply(abundance, 2,
                           FUN = function(x){return(x/sum(x))})

# remove OTUs with less than 5% variability
variability <- apply(relativeAbundance, 2, isVariant)
relativeAbundance <- relativeAbundance[,variability]
relativeAbundance <- data.frame(OTU=rownames(relativeAbundance), relativeAbundance)

# select relevant OTUs by Random forest feature selections
data <- relativeAbundance
otu <- data[,1]
data <- data[,-1]
data <- apply(data, 2, as.numeric)
data <- t(data)
colnames(data) <- otu
y <- sapply(rownames(data), FUN = function(x){substr(x, 1, 1)})
data <- data.frame(data, y=as.factor(y))
boruta <- Boruta(y ~ ., data = data, maxRuns = 1000)

importantOTU <- names(boruta$finalDecision)[boruta$finalDecision %in% c('Confirmed', 'Tentative')]
relativeAbundance <- relativeAbundance[relativeAbundance$OTU %in% importantOTU, ]

write.table(relativeAbundance, file = 'Microbiome_ASD/relativeAbundance.csv',
            row.names = F, quote = F, sep=',')