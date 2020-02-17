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

sn38     <- read.csv('Processed_Clinical_Data/sn38_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

sn38_lines           <- sn38$COSMIC_ID #766


sn38_rna_seq <- gdsc[rownames(gdsc) %in% sn38$COSMIC_ID, ]
sn38_rna_seq <- as.data.frame(sn38_rna_seq)
sn38_rna_seq$res_sens <- sn38$res_sens

sn38_rose <- read.csv('Processed_Gene_Expression/sn38_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)


## glm model
set.seed(5)
sn38_fit_elnet <- cv.glmnet(x = as.matrix(sn38_rose[, -14210]), y = sn38_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/sn38_auc.png')
plot(sn38_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/sn38_model.rds', sn38_fit_elnet)
sn38_pred <- predict(sn38_fit_elnet, newx = as.matrix(sn38_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#sn38_preds <- glmnet::auc(sn38$res_sens, sn38_pred)
#sn38_preds <- round(sn38_preds, digits = 2) #0.66

sn38_overall_acc <- sum(sn38$res_sens == sn38_pred)/length(sn38_pred) #0.916

sn38_glm_confusion <- confusionMatrix(factor(sn38_pred), factor(sn38$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

sn38_svm <- train(factor(res_sens) ~ ., data = sn38_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



sn38_pred <- predict(sn38_svm, newdata = sn38_rna_seq)


sn38_svm_confusion <- confusionMatrix(sn38_pred, factor(sn38$res_sens))

saveRDS(file = 'GLM_Models/sn38_svm_model_1.rds', sn38_svm)



## fit glm mars first degree
set.seed(5)
sn38_fit_1 <- earth(res_sens ~ ., data = sn38_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
sn38_fit_1_summary <- summary(sn38_fit_1)
sn38_fit_1_evimp <- evimp(sn38_fit_1) #plot these w ggplot side by side

sn38_fit_1_pred <- predict(sn38_fit_1, newdata = as.matrix(sn38_rna_seq), type = 'class')

sn38_mars_confusion <- confusionMatrix(factor(sn38_fit_1_pred), factor(sn38$res_sens))


saveRDS(file = 'GLM_Models/sn38_cv_mars_glm_model_1.rds', sn38_fit_1)