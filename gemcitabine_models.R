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

gemcitabine     <- read.csv('Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

gemcitabine_lines           <- gemcitabine$COSMIC_ID #766


gemcitabine_rna_seq <- gdsc[rownames(gdsc) %in% gemcitabine$COSMIC_ID, ]
gemcitabine_rna_seq <- as.data.frame(gemcitabine_rna_seq)
gemcitabine_rna_seq$res_sens <- gemcitabine$res_sens

gemcitabine_rose <- read.csv('Processed_Gene_Expression/gemcitabine_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)


## glm model
# set.seed(5)
# gemcitabine_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rose[, -14210]), y = gemcitabine_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# #save plot
# png(filename = 'Images/gemcitabine_auc.png')
# plot(gemcitabine_fit_elnet)
# dev.off()
# #save model
# saveRDS(file = 'GLM_Models/gemcitabine_model.rds', gemcitabine_fit_elnet)
# gemcitabine_pred <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
# #gemcitabine_preds <- glmnet::auc(gemcitabine$res_sens, gemcitabine_pred)
# #gemcitabine_preds <- round(gemcitabine_preds, digits = 2) #0.66
# 
# gemcitabine_overall_acc <- sum(gemcitabine$res_sens == gemcitabine_pred)/length(gemcitabine_pred) #0.916
# 
# gemcitabine_glm_confusion <- confusionMatrix(factor(gemcitabine_pred), factor(gemcitabine$res_sens))

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

gemcitabine_glm_caret <- train(factor(res_sens) ~ ., data = gemcitabine_rose, 
                   method = 'glmnet',
                   family = 'binomial',
                   trControl = trctrl, 
                   tuneLength = 10)

gemcitabine_glm_caret$bestTune
summary(glm_caret)

gemcitabine_pred <- predict(gemcitabine_glm_caret, newdata = gemcitabine_rna_seq)


gemcitabine_glm_confusion <- confusionMatrix(gemcitabine_pred, factor(gemcitabine$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

gemcitabine_svm <- train(factor(res_sens) ~ ., data = gemcitabine_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



gemcitabine_pred <- predict(gemcitabine_svm, newdata = gemcitabine_rna_seq)


gemcitabine_svm_confusion <- confusionMatrix(gemcitabine_pred, factor(gemcitabine$res_sens))

saveRDS(file = 'GLM_Models/gemcitabine_svm_model_1.rds', gemcitabine_svm)



## fit glm mars first degree
set.seed(5)
gemcitabine_fit_1 <- earth(res_sens ~ ., data = gemcitabine_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
gemcitabine_fit_1_summary <- summary(gemcitabine_fit_1)
gemcitabine_fit_1_evimp <- evimp(gemcitabine_fit_1) #plot these w ggplot side by side

gemcitabine_fit_1_pred <- predict(gemcitabine_fit_1, newdata = as.matrix(gemcitabine_rna_seq), type = 'class')

gemcitabine_mars_confusion <- confusionMatrix(factor(gemcitabine_fit_1_pred), factor(gemcitabine$res_sens))


saveRDS(file = 'GLM_Models/gemcitabine_cv_mars_glm_model_1.rds', gemcitabine_fit_1)