library(caret)
library(e1071)
library(pROC)

# Load Data
data <- read.delim2('Microbiome_ASD/relativeAbundance.csv', sep=',')
otu <- data[,1]
data <- data[,-1]
data <- apply(data, 2, as.numeric)
data <- t(data)
colnames(data) <- otu

# split 10-fold-cross validations
y <- sapply(rownames(data), FUN = function(x){substr(x, 1, 1)})
folds <- createFolds(y)
resSVM <- NULL
resRF <- NULL
# test SVM
for(i in 1:length(folds)) {
  # train model
  train <- data[-folds[[i]],]
  train.y <- y[-folds[[i]]]
  test <- data[folds[[i]],]
  test.y <- y[folds[[i]]]
  df <- data.frame(train, y=train.y)
  df$y <- as.factor(df$y)
  
  ## SVM
  # optimize parameters
  param <- tune.svm(y~., data=df, gamma=2^(-1:1), cost = 2^(2:4), kernel='radial')
  fit <- svm(y~., data=df, kernel = 'radial',
             gamma = param$best.parameters$gamma,
             cost = param$best.parameters$cost, probability=T)
  # compare to test
  prediction <- attributes(predict(fit, test, probability = T))$probabilities
  rocPred <- roc(test.y, prediction[,1])
  resSVM <- c(resSVM, auc(rocPred))
  
  ## Random Forest
  param <- tune.randomForest(y~., data=df, ntree = c(200,500, 700, 1000),
                             mtry = c(0.1*nrow(data),
                                      0.2*nrow(data),
                                      0.3*nrow(data)))
  fit <- randomForest(y~., data=df, keep.forest=T,
                      ntree = param$best.parameters$ntree,
                      mtry=param$best.parameters$mtry)
  prediction <- predict(fit, test, type='prob')
  rocPred <- roc(test.y, prediction[,1])
  resRF <- c(resRF, auc(rocPred))
}


# evaluate
print(paste('average AUC expected for SVM: ', round(mean(resSVM), 3)))
print(paste('average AUC expected for Random Forest: ', round(mean(resRF), 3)))