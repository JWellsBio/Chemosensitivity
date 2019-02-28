
install.packages('Boruta')
library(Boruta)
gdsc_rna_seq <- read.csv('Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE)
gdsc_rna_seq <- gdsc_rna_seq[, -1]
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))

cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)
cisplatin_lines           <- cisplatin$COSMIC_ID #680

gdsc <- data.frame(t(gdsc_rna_seq)) # puts predictors in columns
rownames(gdsc) <- gsub('X', '', rownames(gdsc))
# dim: 962 x 14209 

# split GDSC in half randomly
set.seed(5)
# get random numbers to use for split
random_sample <- sample(x = rownames(gdsc), size = nrow(gdsc)/2)
# create function opposite of %in%
'%ni%' <- Negate('%in%')

# get training and testing sets
gdsc_train         <- gdsc[random_sample, ] #481 x 14209

gdsc_test          <- gdsc[which(rownames(gdsc) %ni% random_sample), ] #481 x 14209

# make sure zero overlap
intersect(rownames(gdsc_train), rownames(gdsc_test))

cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# 338 x 14209
cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]

cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

cisplatin_rna_seq_train$most_sensitive <- as.factor(cisplatin_train$most_sensitive)

set.seed(5)
cisplatin_boruta <- Boruta(most_sensitive ~ ., data = cisplatin_rna_seq_train, doTrace = 3, maxRuns = 5000, getImp = getImpExtraGini)
print(cisplatin_boruta)
cisplatin_i_want <- getSelectedAttributes(cisplatin_boruta, withTentative = T)
getConfirmedFormula(cisplatin_boruta)
plotImpHistory(cisplatin_boruta)

plot(cisplatin_boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(cisplatin_boruta$ImpHistory),function(i)
  cisplatin_boruta$ImpHistory[is.finite(cisplatin_boruta$ImpHistory[,i]),i])
names(lz) <- colnames(cisplatin_boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(cisplatin_boruta$ImpHistory), cex.axis = 0.7)

cisplatin_final_boruta <- TentativeRoughFix(cisplatin_boruta)
print(cisplatin_final_boruta)
plotImpHistory(cisplatin_final_boruta)
getSelectedAttributes(cisplatin_final_boruta, withTentative = F)
cisplatin_boruta_df <- attStats(cisplatin_final_boruta)
print(cisplatin_boruta_df)
cisplatin_i_want <- getSelectedAttributes(cisplatin_final_boruta, withTentative = F)

install.packages('randomForest')
library(randomForest)
cisplatin_rna_seq_train_short <- cisplatin_rna_seq_train[, cisplatin_i_want]
cisplatin_rna_seq_train_short <- apply(cisplatin_rna_seq_train_short, 2, scale)
cisplatin_rna_seq_train_short <- as.data.frame(cisplatin_rna_seq_train_short)
cisplatin_rna_seq_train_short$most_sensitive <- as.factor(cisplatin_train$most_sensitive)

cisplatin_most_rf_model <- randomForest(most_sensitive ~ ., data = cisplatin_rna_seq_train_short, importance = TRUE)
cisplatin_most_rf_model

cisplatin_most_rf_model2 <- randomForest(most_sensitive ~ ., data = cisplatin_rna_seq_train_short, mtry = 5, ntree = 1000, importance = TRUE)
cisplatin_most_rf_model2
varImpPlot(cisplatin_most_rf_model2, lcolor = 'red', color = c('turquoise', 'darkorange'))

install.packages('mlbench')
library(mlbench)
install.packages('caret')
library(caret)

control <- trainControl(method="repeatedcv", number=10, repeats=10)
seed <- 5
metric <- "Accuracy"
set.seed(seed)
mtry <- 5
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(most_sensitive ~ ., data=cisplatin_rna_seq_train_short, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)

control <- trainControl(method="repeatedcv", number=10, repeats=10, search="random")
set.seed(seed)
mtry <- sqrt(ncol(cisplatin_rna_seq_train_short))
rf_random <- train(most_sensitive ~ ., data=cisplatin_rna_seq_train_short, method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)

control <- trainControl(method="repeatedcv", number=10, repeats=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(most_sensitive ~ ., data=cisplatin_rna_seq_train_short, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)

set.seed(seed)
# this one needs to NOT HAVE the most_sensitive factor in the data
bestmtry <- tuneRF(cisplatin_rna_seq_train, as.factor(cisplatin_train$most_sensitive), stepFactor=1.5, improve=1e-5, ntree=500)
print(bestmtry)

rf_new <- train(most_sensitive ~ ., data = cisplatin_rna_seq_train_short, method = 'rf', metric = metric, mtry = 5, tuneGrid = tunegrid, trControl = control)

# TRYING SVM
library(ROCR)
library(e1071)
library(AUC)
model_svm <- svm(most_sensitive ~ ., data = cisplatin_rna_seq_train_short)
pred <- predict(model_svm, newdata = cisplatin_rna_seq_test)
roc_test <- prediction(pred_num, cisplatin_test$most_sensitive)
perf <- performance(roc_test, 'tpr', 'fpr')
plot(perf)
legend(x = 0.7, y = 0.4, legend = c('AUC = 0.91'), lty = 1, col = 'black', bty = 'n')

levels(as.factor(cisplatin_train$most_sensitive))
auc(as.numeric(levels(pred)), cisplatin_rna_seq_train_short$most_sensitive)
pred_num <- as.numeric(pred)
pred_num <- pred_num - 1
perf_auc <- performance(roc_test, 'auc')
perf_auc
model_svm

train <-sim.data(n = 200, ng = 1000, nsg = 10, seed=123)
Lambda.scad <- seq(0.01,0.05, 0.01)
fit.scad<- svmfs(x=t(train$x),y=train$y, fs.method='scad', lambda1.set=Lambda.scad)
test <- sim.data(n = 200, ng = 1000, nsg = 10, seed=124)
predict(fit.scad,newdata=t(test$x), new.data.labels=test$y)





