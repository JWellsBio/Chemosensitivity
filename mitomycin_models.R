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

mitomycin     <- read.csv('Processed_Clinical_Data/mitomycin_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

mitomycin_lines           <- mitomycin$COSMIC_ID #766


mitomycin_rna_seq <- gdsc[rownames(gdsc) %in% mitomycin$COSMIC_ID, ]
mitomycin_rna_seq <- as.data.frame(mitomycin_rna_seq)
mitomycin_rna_seq$res_sens <- mitomycin$res_sens

mitomycin_rose <- read.csv('Processed_Gene_Expression/mitomycin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)


## glm model
set.seed(5)
mitomycin_fit_elnet <- cv.glmnet(x = as.matrix(mitomycin_rose[, -14210]), y = mitomycin_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/mitomycin_auc.png')
plot(mitomycin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/mitomycin_model.rds', mitomycin_fit_elnet)
mitomycin_pred <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#mitomycin_preds <- glmnet::auc(mitomycin$res_sens, mitomycin_pred)
#mitomycin_preds <- round(mitomycin_preds, digits = 2) #0.66

mitomycin_overall_acc <- sum(mitomycin$res_sens == mitomycin_pred)/length(mitomycin_pred) #0.916

mitomycin_glm_confusion <- confusionMatrix(factor(mitomycin_pred), factor(mitomycin$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

mitomycin_svm <- train(factor(res_sens) ~ ., data = mitomycin_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



mitomycin_pred <- predict(mitomycin_svm, newdata = mitomycin_rna_seq)


mitomycin_svm_confusion <- confusionMatrix(mitomycin_pred, factor(mitomycin$res_sens))

saveRDS(file = 'GLM_Models/mitomycin_svm_model_1.rds', mitomycin_svm)



## fit glm mars first degree
set.seed(5)
mitomycin_fit_1 <- earth(res_sens ~ ., data = mitomycin_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
mitomycin_fit_1_summary <- summary(mitomycin_fit_1)
mitomycin_fit_1_evimp <- evimp(mitomycin_fit_1) #plot these w ggplot side by side

mitomycin_fit_1_pred <- predict(mitomycin_fit_1, newdata = as.matrix(mitomycin_rna_seq), type = 'class')

mitomycin_mars_confusion <- confusionMatrix(factor(mitomycin_fit_1_pred), factor(mitomycin$res_sens))


saveRDS(file = 'GLM_Models/mitomycin_cv_mars_glm_model_1.rds', mitomycin_fit_1)