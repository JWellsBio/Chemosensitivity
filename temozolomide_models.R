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

temozolomide     <- read.csv('Processed_Clinical_Data/temozolomide_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

temozolomide_lines           <- temozolomide$COSMIC_ID #766


temozolomide_rna_seq <- gdsc[rownames(gdsc) %in% temozolomide$COSMIC_ID, ]
temozolomide_rna_seq <- as.data.frame(temozolomide_rna_seq)
temozolomide_rna_seq$res_sens <- temozolomide$res_sens

temozolomide_rose <- read.csv('Processed_Gene_Expression/temozolomide_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)


## glm model
set.seed(5)
temozolomide_fit_elnet <- cv.glmnet(x = as.matrix(temozolomide_rose[, -14210]), y = temozolomide_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/temozolomide_auc.png')
plot(temozolomide_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/temozolomide_model.rds', temozolomide_fit_elnet)
temozolomide_pred <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#temozolomide_preds <- glmnet::auc(temozolomide$res_sens, temozolomide_pred)
#temozolomide_preds <- round(temozolomide_preds, digits = 2) #0.66

temozolomide_overall_acc <- sum(temozolomide$res_sens == temozolomide_pred)/length(temozolomide_pred) #0.916

temozolomide_glm_confusion <- confusionMatrix(factor(temozolomide_pred), factor(temozolomide$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

temozolomide_svm <- train(factor(res_sens) ~ ., data = temozolomide_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



temozolomide_pred <- predict(temozolomide_svm, newdata = temozolomide_rna_seq)


temozolomide_svm_confusion <- confusionMatrix(temozolomide_pred, factor(temozolomide$res_sens))

saveRDS(file = 'GLM_Models/temozolomide_svm_model_1.rds', temozolomide_svm)



## fit glm mars first degree
set.seed(5)
temozolomide_fit_1 <- earth(res_sens ~ ., data = temozolomide_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
temozolomide_fit_1_summary <- summary(temozolomide_fit_1)
temozolomide_fit_1_evimp <- evimp(temozolomide_fit_1) #plot these w ggplot side by side

temozolomide_fit_1_pred <- predict(temozolomide_fit_1, newdata = as.matrix(temozolomide_rna_seq), type = 'class')

temozolomide_mars_confusion <- confusionMatrix(factor(temozolomide_fit_1_pred), factor(temozolomide$res_sens))


saveRDS(file = 'GLM_Models/temozolomide_cv_mars_glm_model_1.rds', temozolomide_fit_1)