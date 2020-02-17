## load necessary packages ----
if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) # for model building

if (!require ('ROSE')) install.packages('ROSE')
library(ROSE)

if (!require ('caret')) install.packages('caret')
library(caret)

if (!require(earth)) install.packages('earth')
library(earth)

if (!require('vip')) install.packages('vip')
library(vip)

if (!require('pdp')) install.pacakges('pdp') 
library(pdp)

if (!require('kernlab')) install.packages('kernlab')
library(kernlab)

if (!require('nnet')) install.packages('nnet')
library(nnet)

if (!require('quantreg')) install.packages('quantreg')
library(quantreg)

if(!require('gridExtra')) install.packages('gridExtra')
library(gridExtra)


### functions needed ----
# create function opposite of %in%
'%ni%' <- Negate('%in%')

methotrexate     <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

methotrexate_lines           <- methotrexate$COSMIC_ID #766


methotrexate_rna_seq <- gdsc[rownames(gdsc) %in% methotrexate$COSMIC_ID, ]
methotrexate_rna_seq <- as.data.frame(methotrexate_rna_seq)
methotrexate_rna_seq$res_sens <- methotrexate$res_sens

methotrexate_rose <- read.csv('Processed_Gene_Expression/methotrexate_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)


## glm model
set.seed(5)
methotrexate_fit_elnet <- cv.glmnet(x = as.matrix(methotrexate_rose[, -14210]), y = methotrexate_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/methotrexate_auc.png')
plot(methotrexate_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/methotrexate_model.rds', methotrexate_fit_elnet)
methotrexate_pred <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#methotrexate_preds <- glmnet::auc(methotrexate$res_sens, methotrexate_pred)
#methotrexate_preds <- round(methotrexate_preds, digits = 2) #0.66

methotrexate_overall_acc <- sum(methotrexate$res_sens == methotrexate_pred)/length(methotrexate_pred) #0.916

methotrexate_glm_confusion <- confusionMatrix(factor(methotrexate_pred), factor(methotrexate$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

methotrexate_svm <- train(factor(res_sens) ~ ., data = methotrexate_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



methotrexate_pred <- predict(methotrexate_svm, newdata = methotrexate_rna_seq)


methotrexate_svm_confusion <- confusionMatrix(methotrexate_pred, factor(methotrexate$res_sens))

saveRDS(file = 'GLM_Models/methotrexate_svm_model_1.rds', methotrexate_svm)



## fit glm mars first degree
set.seed(5)
methotrexate_fit_1 <- earth(res_sens ~ ., data = methotrexate_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
methotrexate_fit_1_summary <- summary(methotrexate_fit_1)
methotrexate_fit_1_evimp <- evimp(methotrexate_fit_1) #plot these w ggplot side by side

methotrexate_fit_1_pred <- predict(methotrexate_fit_1, newdata = as.matrix(methotrexate_rna_seq), type = 'class')

methotrexate_mars_confusion <- confusionMatrix(factor(methotrexate_fit_1_pred), factor(methotrexate$res_sens))


saveRDS(file = 'GLM_Models/methotrexate_cv_mars_glm_model_1.rds', methotrexate_fit_1)