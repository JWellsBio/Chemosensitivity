## THIS SCRIPT GETS AUC VALUES FOR OVERALL AND BY CANCER TYPE

## load necessary packages ----
if (!require ('ggplot2')) install.packages('ggplot2')
library(ggplot2) # general plotting functions

if (!require ('ROCR')) install.packages('ROCR')
library(ROCR) # for drawing AUC curves

if (!require ('gplots')) install.packages('gplots')
library(gplots) # for heatmap.2 functionality

if (!require ('survival')) install.packages('survival')
library(survival) # for survival curve functions

if (!require ('survminer')) install.packages('survminer')
library(survminer) # for additional survival curve design

if (!require ('ComplexHeatmap')) BiocInstaller::biocLite('ComplexHeatmap')
library(ComplexHeatmap) # for increased heatmap design

if (!require ('formattable')) install.packages('formattable')
library(formattable) # for table formatting and output

if(!require ('htmltools')) install.packages('htmltools')
library(htmltools) # to support formattable functions

if (!require ('webshot')) install.packages('webshot')
library(webshot) # to support formattavle functions

if (!require ('glmnet')) install.packages('glmnet')
library(glmnet)

library(survminer)
library(survival)

### GDSC ------
## load clinical data ----
cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)
etoposide     <- read.csv('Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv', row.names = 1)
gemcitabine   <- read.csv('Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv', row.names = 1)
methotrexate  <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc_rna_seq <- read.csv('Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE)
gdsc_rna_seq <- gdsc_rna_seq[, -1]
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))

### set up data for model building ----------
# get names of GDSC cell lines treated with each drug
cisplatin_lines           <- cisplatin$COSMIC_ID #680
etoposide_lines           <- etoposide$COSMIC_ID #718
methotrexate_lines        <- methotrexate$COSMIC_ID # 679

# set GDSC to usable format
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

# create training/testing sets for each drug
cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# 338 x 14209
cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]
# 342 x 14209
etoposide_rna_seq_train     <- gdsc_train[intersect(etoposide_lines, rownames(gdsc_train)), ]
# 352 x 14209
etoposide_rna_seq_test      <- gdsc_test[intersect(etoposide_lines, rownames(gdsc_test)), ]
# 366 x 14209
methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# 337 x 14209
methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# 342 x 14209

# split clinical data
cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]

methotrexate_train     <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test      <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]

# scale data
cisplatin_rna_seq_train_scaled          <- apply(cisplatin_rna_seq_train, 2, scale)
cisplatin_rna_seq_test_scaled           <- as.data.frame(apply(cisplatin_rna_seq_test, 2, scale))

etoposide_rna_seq_train_scaled          <- apply(etoposide_rna_seq_train, 2, scale)
etoposide_rna_seq_test_scaled           <- as.data.frame(apply(etoposide_rna_seq_test, 2, scale))

methotrexate_rna_seq_train_scaled       <- apply(methotrexate_rna_seq_train, 2, scale)
methotrexate_rna_seq_test_scaled        <- as.data.frame(apply(methotrexate_rna_seq_test, 2, scale))

### load models ----
cisplatin_most_fit_elnet      <- readRDS('GLM_Models/cisplatin_most_model.rds')
cisplatin_least_fit_elnet     <- readRDS('GLM_Models/cisplatin_least_model.rds')

etoposide_most_fit_elnet      <- readRDS('GLM_Models/etoposide_most_model.rds')
etoposide_least_fit_elnet     <- readRDS('GLM_Models/etoposide_least_model.rds')

methotrexate_most_fit_elnet   <- readRDS('GLM_Models/methotrexate_most_model.rds')
methotrexate_least_fit_elnet  <- readRDS('GLM_Models/methotrexate_least_model.rds')

### get accuracy on testing data --------

## CISPLATIN
new_cisplatin_1se <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_test_gdsc_auc_1se <- auc(cisplatin_test$res_sens, new_cisplatin_1se)
cisplatin_test_gdsc_auc_1se <- round(cisplatin_test_gdsc_auc_1se, digits = 2)


## ETOPOSIDE
new_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_most_test_gdsc_auc_min <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_min)
etoposide_most_test_gdsc_auc_min <- round(etoposide_most_test_gdsc_auc_min, digits = 2)


new_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_least_test_gdsc_auc_min <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_min)
etoposide_least_test_gdsc_auc_min <- round(etoposide_least_test_gdsc_auc_min, digits = 2)


## METHOTREXATE
new_methotrexate_most_sensitive_min <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_most_test_gdsc_auc_min <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_min)
methotrexate_most_test_gdsc_auc_min <- round(methotrexate_most_test_gdsc_auc_min, digits = 2)


new_methotrexate_least_sensitive_1se <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_least_test_gdsc_auc_1se <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_1se)
methotrexate_least_test_gdsc_auc_1se <- round(methotrexate_least_test_gdsc_auc_1se, digits = 2)

#put them together
overall_auc <- c(bleomycin_preds, camptothecin_preds, cisplatin_preds, cytarabine_preds, doxorubicin_preds, 
                 etoposide_preds, gemcitabine_preds, methotrexate_preds, mitomycin_preds, 
                 sn38_preds, temozolomide_preds)

# subsetting lines by cancer type
bleomycin_aero_dig_tract_lines              <- bleomycin_test$Cell_line_tissue_type == 'aero_dig_tract'
bleomycin_bone_lines                        <- bleomycin_test$Cell_line_tissue_type == 'bone'
bleomycin_breast_lines                      <- bleomycin_test$Cell_line_tissue_type == 'breast'
bleomycin_digestive_system_lines            <- bleomycin_test$Cell_line_tissue_type == 'digestive_system'
bleomycin_kidney_lines                      <- bleomycin_test$Cell_line_tissue_type == 'kidney'
bleomycin_large_intestine_lines             <- bleomycin_test$Cell_line_tissue_type == 'large_intestine'
bleomycin_lung_lines                        <- bleomycin_test$Cell_line_tissue_type == 'lung'
bleomycin_lung_NSCLC_lines                  <- bleomycin_test$Cell_line_tissue_type == 'lung_NSCLC'
bleomycin_lung_SCLC_lines                   <- bleomycin_test$Cell_line_tissue_type == 'lung_SCLC'
bleomycin_nervous_system_lines              <- bleomycin_test$Cell_line_tissue_type == 'nervous_system'
bleomycin_neuroblastoma_lines               <- bleomycin_test$Cell_line_tissue_type == 'neuroblastoma'
bleomycin_pancreas_lines                    <- bleomycin_test$Cell_line_tissue_type == 'pancreas'
bleomycin_skin_lines                        <- bleomycin_test$Cell_line_tissue_type == 'skin'
bleomycin_soft_tissue_lines                 <- bleomycin_test$Cell_line_tissue_type == 'soft_tissue'
bleomycin_thyroid_lines                     <- bleomycin_test$Cell_line_tissue_type == 'thyroid'
bleomycin_urogenital_system_lines           <- bleomycin_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
bleomycin_test_aero_dig_tract_exp <- bleomycin_rna_seq_test[bleomycin_aero_dig_tract_lines, ]
bleomycin_test_aero_dig_tract <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_aero_dig_tract_auc <- sum(bleomycin_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_bone_exp <- bleomycin_rna_seq_test[bleomycin_bone_lines, ]
bleomycin_test_bone <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_bone_auc <- sum(bleomycin_test_bone$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_breast_exp <- bleomycin_rna_seq_test[bleomycin_breast_lines, ]
bleomycin_test_breast <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_breast_auc <- sum(bleomycin_test_breast$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_digestive_system_exp <- bleomycin_rna_seq_test[bleomycin_digestive_system_lines, ]
bleomycin_test_digestive_system <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_digestive_system_auc <- sum(bleomycin_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_kidney_exp <- bleomycin_rna_seq_test[bleomycin_kidney_lines, ]
bleomycin_test_kidney <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_kidney_auc <- sum(bleomycin_test_kidney$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_large_intestine_exp <- bleomycin_rna_seq_test[bleomycin_large_intestine_lines, ]
bleomycin_test_large_intestine <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_large_intestine_auc <- sum(bleomycin_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_lung_exp <- bleomycin_rna_seq_test[bleomycin_lung_lines, ]
bleomycin_test_lung <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_lung_auc <- sum(bleomycin_test_lung$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_lung_NSCLC_exp <- bleomycin_rna_seq_test[bleomycin_lung_NSCLC_lines, ]
bleomycin_test_lung_NSCLC <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_lung_NSCLC_auc <- sum(bleomycin_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_lung_SCLC_exp <- bleomycin_rna_seq_test[bleomycin_lung_SCLC_lines, ]
bleomycin_test_lung_SCLC <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_lung_SCLC_auc <- sum(bleomycin_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_nervous_system_exp <- bleomycin_rna_seq_test[bleomycin_nervous_system_lines, ]
bleomycin_test_nervous_system <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_nervous_system_auc <- sum(bleomycin_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_neuroblastoma_exp <- bleomycin_rna_seq_test[bleomycin_neuroblastoma_lines, ]
bleomycin_test_neuroblastoma <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_neuroblastoma_auc <- sum(bleomycin_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_pancreas_exp <- bleomycin_rna_seq_test[bleomycin_pancreas_lines, ]
bleomycin_test_pancreas <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_pancreas_auc <- sum(bleomycin_test_pancreas$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_skin_exp <- bleomycin_rna_seq_test[bleomycin_skin_lines, ]
bleomycin_test_skin <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_skin_auc <- sum(bleomycin_test_skin$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_soft_tissue_exp <- bleomycin_rna_seq_test[bleomycin_soft_tissue_lines, ]
bleomycin_test_soft_tissue <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_soft_tissue_auc <- sum(bleomycin_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_thyroid_exp <- bleomycin_rna_seq_test[bleomycin_thyroid_lines, ]
bleomycin_test_thyroid <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_thyroid_auc <- sum(bleomycin_test_thyroid$res_sens == new_ic50)/length(new_ic50)

bleomycin_test_urogenital_system_exp <- bleomycin_rna_seq_test[bleomycin_urogenital_system_lines, ]
bleomycin_test_urogenital_system <- bleomycin_test[which(bleomycin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
bleomycin_urogenital_system_auc <- sum(bleomycin_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

bleomycin_auc <- c(bleomycin_aero_dig_tract_auc, bleomycin_bone_auc, 
                   bleomycin_breast_auc, bleomycin_digestive_system_auc, 
                   bleomycin_kidney_auc, bleomycin_large_intestine_auc, 
                   bleomycin_lung_auc, bleomycin_lung_NSCLC_auc, bleomycin_lung_SCLC_auc, 
                   bleomycin_nervous_system_auc, bleomycin_neuroblastoma_auc, 
                   bleomycin_pancreas_auc, bleomycin_skin_auc, 
                   bleomycin_soft_tissue_auc, bleomycin_thyroid_auc, 
                   bleomycin_urogenital_system_auc)

camptothecin_aero_dig_tract_lines              <- camptothecin_test$Cell_line_tissue_type == 'aero_dig_tract'
camptothecin_bone_lines                        <- camptothecin_test$Cell_line_tissue_type == 'bone'
camptothecin_breast_lines                      <- camptothecin_test$Cell_line_tissue_type == 'breast'
camptothecin_digestive_system_lines            <- camptothecin_test$Cell_line_tissue_type == 'digestive_system'
camptothecin_kidney_lines                      <- camptothecin_test$Cell_line_tissue_type == 'kidney'
camptothecin_large_intestine_lines             <- camptothecin_test$Cell_line_tissue_type == 'large_intestine'
camptothecin_lung_lines                        <- camptothecin_test$Cell_line_tissue_type == 'lung'
camptothecin_lung_NSCLC_lines                  <- camptothecin_test$Cell_line_tissue_type == 'lung_NSCLC'
camptothecin_lung_SCLC_lines                   <- camptothecin_test$Cell_line_tissue_type == 'lung_SCLC'
camptothecin_nervous_system_lines              <- camptothecin_test$Cell_line_tissue_type == 'nervous_system'
camptothecin_neuroblastoma_lines               <- camptothecin_test$Cell_line_tissue_type == 'neuroblastoma'
camptothecin_pancreas_lines                    <- camptothecin_test$Cell_line_tissue_type == 'pancreas'
camptothecin_skin_lines                        <- camptothecin_test$Cell_line_tissue_type == 'skin'
camptothecin_soft_tissue_lines                 <- camptothecin_test$Cell_line_tissue_type == 'soft_tissue'
camptothecin_thyroid_lines                     <- camptothecin_test$Cell_line_tissue_type == 'thyroid'
camptothecin_urogenital_system_lines           <- camptothecin_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
camptothecin_test_aero_dig_tract_exp <- camptothecin_rna_seq_test[camptothecin_aero_dig_tract_lines, ]
camptothecin_test_aero_dig_tract <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_aero_dig_tract_auc <- sum(camptothecin_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_bone_exp <- camptothecin_rna_seq_test[camptothecin_bone_lines, ]
camptothecin_test_bone <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_bone_auc <- sum(camptothecin_test_bone$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_breast_exp <- camptothecin_rna_seq_test[camptothecin_breast_lines, ]
camptothecin_test_breast <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_breast_auc <- sum(camptothecin_test_breast$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_digestive_system_exp <- camptothecin_rna_seq_test[camptothecin_digestive_system_lines, ]
camptothecin_test_digestive_system <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_digestive_system_auc <- sum(camptothecin_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_kidney_exp <- camptothecin_rna_seq_test[camptothecin_kidney_lines, ]
camptothecin_test_kidney <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_kidney_auc <- sum(camptothecin_test_kidney$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_large_intestine_exp <- camptothecin_rna_seq_test[camptothecin_large_intestine_lines, ]
camptothecin_test_large_intestine <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_large_intestine_auc <- sum(camptothecin_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_lung_exp <- camptothecin_rna_seq_test[camptothecin_lung_lines, ]
camptothecin_test_lung <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_lung_auc <- sum(camptothecin_test_lung$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_lung_NSCLC_exp <- camptothecin_rna_seq_test[camptothecin_lung_NSCLC_lines, ]
camptothecin_test_lung_NSCLC <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_lung_NSCLC_auc <- sum(camptothecin_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_lung_SCLC_exp <- camptothecin_rna_seq_test[camptothecin_lung_SCLC_lines, ]
camptothecin_test_lung_SCLC <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_lung_SCLC_auc <- sum(camptothecin_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_nervous_system_exp <- camptothecin_rna_seq_test[camptothecin_nervous_system_lines, ]
camptothecin_test_nervous_system <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_nervous_system_auc <- sum(camptothecin_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_neuroblastoma_exp <- camptothecin_rna_seq_test[camptothecin_neuroblastoma_lines, ]
camptothecin_test_neuroblastoma <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_neuroblastoma_auc <- sum(camptothecin_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_pancreas_exp <- camptothecin_rna_seq_test[camptothecin_pancreas_lines, ]
camptothecin_test_pancreas <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_pancreas_auc <- sum(camptothecin_test_pancreas$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_skin_exp <- camptothecin_rna_seq_test[camptothecin_skin_lines, ]
camptothecin_test_skin <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_skin_auc <- sum(camptothecin_test_skin$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_soft_tissue_exp <- camptothecin_rna_seq_test[camptothecin_soft_tissue_lines, ]
camptothecin_test_soft_tissue <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_soft_tissue_auc <- sum(camptothecin_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_thyroid_exp <- camptothecin_rna_seq_test[camptothecin_thyroid_lines, ]
camptothecin_test_thyroid <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_thyroid_auc <- sum(camptothecin_test_thyroid$res_sens == new_ic50)/length(new_ic50)

camptothecin_test_urogenital_system_exp <- camptothecin_rna_seq_test[camptothecin_urogenital_system_lines, ]
camptothecin_test_urogenital_system <- camptothecin_test[which(camptothecin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
camptothecin_urogenital_system_auc <- sum(camptothecin_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

camptothecin_auc <- c(camptothecin_aero_dig_tract_auc, camptothecin_bone_auc, 
                   camptothecin_breast_auc, camptothecin_digestive_system_auc, 
                   camptothecin_kidney_auc, camptothecin_large_intestine_auc, 
                   camptothecin_lung_auc, camptothecin_lung_NSCLC_auc, camptothecin_lung_SCLC_auc, 
                   camptothecin_nervous_system_auc, camptothecin_neuroblastoma_auc, 
                   camptothecin_pancreas_auc, camptothecin_skin_auc, 
                   camptothecin_soft_tissue_auc, camptothecin_thyroid_auc, 
                   camptothecin_urogenital_system_auc)


cisplatin_aero_dig_tract_lines              <- cisplatin_test$Cell_line_tissue_type == 'aero_dig_tract'
cisplatin_bone_lines                        <- cisplatin_test$Cell_line_tissue_type == 'bone'
cisplatin_breast_lines                      <- cisplatin_test$Cell_line_tissue_type == 'breast'
cisplatin_digestive_system_lines            <- cisplatin_test$Cell_line_tissue_type == 'digestive_system'
cisplatin_kidney_lines                      <- cisplatin_test$Cell_line_tissue_type == 'kidney'
cisplatin_large_intestine_lines             <- cisplatin_test$Cell_line_tissue_type == 'large_intestine'
cisplatin_lung_lines                        <- cisplatin_test$Cell_line_tissue_type == 'lung'
cisplatin_lung_NSCLC_lines                  <- cisplatin_test$Cell_line_tissue_type == 'lung_NSCLC'
cisplatin_lung_SCLC_lines                   <- cisplatin_test$Cell_line_tissue_type == 'lung_SCLC'
cisplatin_nervous_system_lines              <- cisplatin_test$Cell_line_tissue_type == 'nervous_system'
cisplatin_neuroblastoma_lines               <- cisplatin_test$Cell_line_tissue_type == 'neuroblastoma'
cisplatin_pancreas_lines                    <- cisplatin_test$Cell_line_tissue_type == 'pancreas'
cisplatin_skin_lines                        <- cisplatin_test$Cell_line_tissue_type == 'skin'
cisplatin_soft_tissue_lines                 <- cisplatin_test$Cell_line_tissue_type == 'soft_tissue'
cisplatin_thyroid_lines                     <- cisplatin_test$Cell_line_tissue_type == 'thyroid'
cisplatin_urogenital_system_lines           <- cisplatin_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
cisplatin_test_aero_dig_tract_exp <- cisplatin_rna_seq_test[cisplatin_aero_dig_tract_lines, ]
cisplatin_test_aero_dig_tract <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_aero_dig_tract_auc <- sum(cisplatin_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_bone_exp <- cisplatin_rna_seq_test[cisplatin_bone_lines, ]
cisplatin_test_bone <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_bone_auc <- sum(cisplatin_test_bone$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_breast_exp <- cisplatin_rna_seq_test[cisplatin_breast_lines, ]
cisplatin_test_breast <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_breast_auc <- sum(cisplatin_test_breast$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_digestive_system_exp <- cisplatin_rna_seq_test[cisplatin_digestive_system_lines, ]
cisplatin_test_digestive_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_digestive_system_auc <- sum(cisplatin_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_kidney_exp <- cisplatin_rna_seq_test[cisplatin_kidney_lines, ]
cisplatin_test_kidney <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_kidney_auc <- sum(cisplatin_test_kidney$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_large_intestine_exp <- cisplatin_rna_seq_test[cisplatin_large_intestine_lines, ]
cisplatin_test_large_intestine <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_large_intestine_auc <- sum(cisplatin_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_lung_exp <- cisplatin_rna_seq_test[cisplatin_lung_lines, ]
cisplatin_test_lung <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_lung_auc <- sum(cisplatin_test_lung$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_lung_NSCLC_exp <- cisplatin_rna_seq_test[cisplatin_lung_NSCLC_lines, ]
cisplatin_test_lung_NSCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_lung_NSCLC_auc <- sum(cisplatin_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_lung_SCLC_exp <- cisplatin_rna_seq_test[cisplatin_lung_SCLC_lines, ]
cisplatin_test_lung_SCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_lung_SCLC_auc <- sum(cisplatin_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_nervous_system_exp <- cisplatin_rna_seq_test[cisplatin_nervous_system_lines, ]
cisplatin_test_nervous_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_nervous_system_auc <- sum(cisplatin_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_neuroblastoma_exp <- cisplatin_rna_seq_test[cisplatin_neuroblastoma_lines, ]
cisplatin_test_neuroblastoma <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_neuroblastoma_auc <- sum(cisplatin_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_pancreas_exp <- cisplatin_rna_seq_test[cisplatin_pancreas_lines, ]
cisplatin_test_pancreas <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_pancreas_auc <- sum(cisplatin_test_pancreas$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_skin_exp <- cisplatin_rna_seq_test[cisplatin_skin_lines, ]
cisplatin_test_skin <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_skin_auc <- sum(cisplatin_test_skin$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_soft_tissue_exp <- cisplatin_rna_seq_test[cisplatin_soft_tissue_lines, ]
cisplatin_test_soft_tissue <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_soft_tissue_auc <- sum(cisplatin_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_thyroid_exp <- cisplatin_rna_seq_test[cisplatin_thyroid_lines, ]
cisplatin_test_thyroid <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_thyroid_auc <- sum(cisplatin_test_thyroid$res_sens == new_ic50)/length(new_ic50)

cisplatin_test_urogenital_system_exp <- cisplatin_rna_seq_test[cisplatin_urogenital_system_lines, ]
cisplatin_test_urogenital_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cisplatin_urogenital_system_auc <- sum(cisplatin_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

cisplatin_auc <- c(cisplatin_aero_dig_tract_auc, cisplatin_bone_auc, 
                   cisplatin_breast_auc, cisplatin_digestive_system_auc, 
                   cisplatin_kidney_auc, cisplatin_large_intestine_auc, 
                   cisplatin_lung_auc, cisplatin_lung_NSCLC_auc, cisplatin_lung_SCLC_auc, 
                   cisplatin_nervous_system_auc, cisplatin_neuroblastoma_auc, 
                   cisplatin_pancreas_auc, cisplatin_skin_auc, 
                   cisplatin_soft_tissue_auc, cisplatin_thyroid_auc, 
                   cisplatin_urogenital_system_auc)



cytarabine_aero_dig_tract_lines              <- cytarabine_test$Cell_line_tissue_type == 'aero_dig_tract'
cytarabine_bone_lines                        <- cytarabine_test$Cell_line_tissue_type == 'bone'
cytarabine_breast_lines                      <- cytarabine_test$Cell_line_tissue_type == 'breast'
cytarabine_digestive_system_lines            <- cytarabine_test$Cell_line_tissue_type == 'digestive_system'
cytarabine_kidney_lines                      <- cytarabine_test$Cell_line_tissue_type == 'kidney'
cytarabine_large_intestine_lines             <- cytarabine_test$Cell_line_tissue_type == 'large_intestine'
cytarabine_lung_lines                        <- cytarabine_test$Cell_line_tissue_type == 'lung'
cytarabine_lung_NSCLC_lines                  <- cytarabine_test$Cell_line_tissue_type == 'lung_NSCLC'
cytarabine_lung_SCLC_lines                   <- cytarabine_test$Cell_line_tissue_type == 'lung_SCLC'
cytarabine_nervous_system_lines              <- cytarabine_test$Cell_line_tissue_type == 'nervous_system'
cytarabine_neuroblastoma_lines               <- cytarabine_test$Cell_line_tissue_type == 'neuroblastoma'
cytarabine_pancreas_lines                    <- cytarabine_test$Cell_line_tissue_type == 'pancreas'
cytarabine_skin_lines                        <- cytarabine_test$Cell_line_tissue_type == 'skin'
cytarabine_soft_tissue_lines                 <- cytarabine_test$Cell_line_tissue_type == 'soft_tissue'
cytarabine_thyroid_lines                     <- cytarabine_test$Cell_line_tissue_type == 'thyroid'
cytarabine_urogenital_system_lines           <- cytarabine_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
cytarabine_test_aero_dig_tract_exp <- cytarabine_rna_seq_test[cytarabine_aero_dig_tract_lines, ]
cytarabine_test_aero_dig_tract <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_aero_dig_tract_auc <- sum(cytarabine_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_bone_exp <- cytarabine_rna_seq_test[cytarabine_bone_lines, ]
cytarabine_test_bone <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_bone_auc <- sum(cytarabine_test_bone$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_breast_exp <- cytarabine_rna_seq_test[cytarabine_breast_lines, ]
cytarabine_test_breast <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_breast_auc <- sum(cytarabine_test_breast$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_digestive_system_exp <- cytarabine_rna_seq_test[cytarabine_digestive_system_lines, ]
cytarabine_test_digestive_system <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_digestive_system_auc <- sum(cytarabine_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_kidney_exp <- cytarabine_rna_seq_test[cytarabine_kidney_lines, ]
cytarabine_test_kidney <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_kidney_auc <- sum(cytarabine_test_kidney$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_large_intestine_exp <- cytarabine_rna_seq_test[cytarabine_large_intestine_lines, ]
cytarabine_test_large_intestine <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_large_intestine_auc <- sum(cytarabine_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_lung_exp <- cytarabine_rna_seq_test[cytarabine_lung_lines, ]
cytarabine_test_lung <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_lung_auc <- sum(cytarabine_test_lung$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_lung_NSCLC_exp <- cytarabine_rna_seq_test[cytarabine_lung_NSCLC_lines, ]
cytarabine_test_lung_NSCLC <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_lung_NSCLC_auc <- sum(cytarabine_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_lung_SCLC_exp <- cytarabine_rna_seq_test[cytarabine_lung_SCLC_lines, ]
cytarabine_test_lung_SCLC <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_lung_SCLC_auc <- sum(cytarabine_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_nervous_system_exp <- cytarabine_rna_seq_test[cytarabine_nervous_system_lines, ]
cytarabine_test_nervous_system <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_nervous_system_auc <- sum(cytarabine_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_neuroblastoma_exp <- cytarabine_rna_seq_test[cytarabine_neuroblastoma_lines, ]
cytarabine_test_neuroblastoma <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_neuroblastoma_auc <- sum(cytarabine_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_pancreas_exp <- cytarabine_rna_seq_test[cytarabine_pancreas_lines, ]
cytarabine_test_pancreas <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_pancreas_auc <- sum(cytarabine_test_pancreas$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_skin_exp <- cytarabine_rna_seq_test[cytarabine_skin_lines, ]
cytarabine_test_skin <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_skin_auc <- sum(cytarabine_test_skin$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_soft_tissue_exp <- cytarabine_rna_seq_test[cytarabine_soft_tissue_lines, ]
cytarabine_test_soft_tissue <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_soft_tissue_auc <- sum(cytarabine_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_thyroid_exp <- cytarabine_rna_seq_test[cytarabine_thyroid_lines, ]
cytarabine_test_thyroid <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_thyroid_auc <- sum(cytarabine_test_thyroid$res_sens == new_ic50)/length(new_ic50)

cytarabine_test_urogenital_system_exp <- cytarabine_rna_seq_test[cytarabine_urogenital_system_lines, ]
cytarabine_test_urogenital_system <- cytarabine_test[which(cytarabine_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
cytarabine_urogenital_system_auc <- sum(cytarabine_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

cytarabine_auc <- c(cytarabine_aero_dig_tract_auc, cytarabine_bone_auc, 
                   cytarabine_breast_auc, cytarabine_digestive_system_auc, 
                   cytarabine_kidney_auc, cytarabine_large_intestine_auc, 
                   cytarabine_lung_auc, cytarabine_lung_NSCLC_auc, cytarabine_lung_SCLC_auc, 
                   cytarabine_nervous_system_auc, cytarabine_neuroblastoma_auc, 
                   cytarabine_pancreas_auc, cytarabine_skin_auc, 
                   cytarabine_soft_tissue_auc, cytarabine_thyroid_auc, 
                   cytarabine_urogenital_system_auc)


doxorubicin_aero_dig_tract_lines              <- doxorubicin_test$Cell_line_tissue_type == 'aero_dig_tract'
doxorubicin_bone_lines                        <- doxorubicin_test$Cell_line_tissue_type == 'bone'
doxorubicin_breast_lines                      <- doxorubicin_test$Cell_line_tissue_type == 'breast'
doxorubicin_digestive_system_lines            <- doxorubicin_test$Cell_line_tissue_type == 'digestive_system'
doxorubicin_kidney_lines                      <- doxorubicin_test$Cell_line_tissue_type == 'kidney'
doxorubicin_large_intestine_lines             <- doxorubicin_test$Cell_line_tissue_type == 'large_intestine'
doxorubicin_lung_lines                        <- doxorubicin_test$Cell_line_tissue_type == 'lung'
doxorubicin_lung_NSCLC_lines                  <- doxorubicin_test$Cell_line_tissue_type == 'lung_NSCLC'
doxorubicin_lung_SCLC_lines                   <- doxorubicin_test$Cell_line_tissue_type == 'lung_SCLC'
doxorubicin_nervous_system_lines              <- doxorubicin_test$Cell_line_tissue_type == 'nervous_system'
doxorubicin_neuroblastoma_lines               <- doxorubicin_test$Cell_line_tissue_type == 'neuroblastoma'
doxorubicin_pancreas_lines                    <- doxorubicin_test$Cell_line_tissue_type == 'pancreas'
doxorubicin_skin_lines                        <- doxorubicin_test$Cell_line_tissue_type == 'skin'
doxorubicin_soft_tissue_lines                 <- doxorubicin_test$Cell_line_tissue_type == 'soft_tissue'
doxorubicin_thyroid_lines                     <- doxorubicin_test$Cell_line_tissue_type == 'thyroid'
doxorubicin_urogenital_system_lines           <- doxorubicin_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
doxorubicin_test_aero_dig_tract_exp <- doxorubicin_rna_seq_test[doxorubicin_aero_dig_tract_lines, ]
doxorubicin_test_aero_dig_tract <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_aero_dig_tract_auc <- sum(doxorubicin_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_bone_exp <- doxorubicin_rna_seq_test[doxorubicin_bone_lines, ]
doxorubicin_test_bone <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_bone_auc <- sum(doxorubicin_test_bone$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_breast_exp <- doxorubicin_rna_seq_test[doxorubicin_breast_lines, ]
doxorubicin_test_breast <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_breast_auc <- sum(doxorubicin_test_breast$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_digestive_system_exp <- doxorubicin_rna_seq_test[doxorubicin_digestive_system_lines, ]
doxorubicin_test_digestive_system <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_digestive_system_auc <- sum(doxorubicin_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_kidney_exp <- doxorubicin_rna_seq_test[doxorubicin_kidney_lines, ]
doxorubicin_test_kidney <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_kidney_auc <- sum(doxorubicin_test_kidney$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_large_intestine_exp <- doxorubicin_rna_seq_test[doxorubicin_large_intestine_lines, ]
doxorubicin_test_large_intestine <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_large_intestine_auc <- sum(doxorubicin_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_lung_exp <- doxorubicin_rna_seq_test[doxorubicin_lung_lines, ]
doxorubicin_test_lung <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_lung_auc <- sum(doxorubicin_test_lung$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_lung_NSCLC_exp <- doxorubicin_rna_seq_test[doxorubicin_lung_NSCLC_lines, ]
doxorubicin_test_lung_NSCLC <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_lung_NSCLC_auc <- sum(doxorubicin_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_lung_SCLC_exp <- doxorubicin_rna_seq_test[doxorubicin_lung_SCLC_lines, ]
doxorubicin_test_lung_SCLC <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_lung_SCLC_auc <- sum(doxorubicin_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_nervous_system_exp <- doxorubicin_rna_seq_test[doxorubicin_nervous_system_lines, ]
doxorubicin_test_nervous_system <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_nervous_system_auc <- sum(doxorubicin_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_neuroblastoma_exp <- doxorubicin_rna_seq_test[doxorubicin_neuroblastoma_lines, ]
doxorubicin_test_neuroblastoma <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_neuroblastoma_auc <- sum(doxorubicin_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_pancreas_exp <- doxorubicin_rna_seq_test[doxorubicin_pancreas_lines, ]
doxorubicin_test_pancreas <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_pancreas_auc <- sum(doxorubicin_test_pancreas$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_skin_exp <- doxorubicin_rna_seq_test[doxorubicin_skin_lines, ]
doxorubicin_test_skin <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_skin_auc <- sum(doxorubicin_test_skin$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_soft_tissue_exp <- doxorubicin_rna_seq_test[doxorubicin_soft_tissue_lines, ]
doxorubicin_test_soft_tissue <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_soft_tissue_auc <- sum(doxorubicin_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_thyroid_exp <- doxorubicin_rna_seq_test[doxorubicin_thyroid_lines, ]
doxorubicin_test_thyroid <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_thyroid_auc <- sum(doxorubicin_test_thyroid$res_sens == new_ic50)/length(new_ic50)

doxorubicin_test_urogenital_system_exp <- doxorubicin_rna_seq_test[doxorubicin_urogenital_system_lines, ]
doxorubicin_test_urogenital_system <- doxorubicin_test[which(doxorubicin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
doxorubicin_urogenital_system_auc <- sum(doxorubicin_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

doxorubicin_auc <- c(doxorubicin_aero_dig_tract_auc, doxorubicin_bone_auc, 
                   doxorubicin_breast_auc, doxorubicin_digestive_system_auc, 
                   doxorubicin_kidney_auc, doxorubicin_large_intestine_auc, 
                   doxorubicin_lung_auc, doxorubicin_lung_NSCLC_auc, doxorubicin_lung_SCLC_auc, 
                   doxorubicin_nervous_system_auc, doxorubicin_neuroblastoma_auc, 
                   doxorubicin_pancreas_auc, doxorubicin_skin_auc, 
                   doxorubicin_soft_tissue_auc, doxorubicin_thyroid_auc, 
                   doxorubicin_urogenital_system_auc)


etoposide_aero_dig_tract_lines              <- etoposide_test$Cell_line_tissue_type == 'aero_dig_tract'
etoposide_bone_lines                        <- etoposide_test$Cell_line_tissue_type == 'bone'
etoposide_breast_lines                      <- etoposide_test$Cell_line_tissue_type == 'breast'
etoposide_digestive_system_lines            <- etoposide_test$Cell_line_tissue_type == 'digestive_system'
etoposide_kidney_lines                      <- etoposide_test$Cell_line_tissue_type == 'kidney'
etoposide_large_intestine_lines             <- etoposide_test$Cell_line_tissue_type == 'large_intestine'
etoposide_lung_lines                        <- etoposide_test$Cell_line_tissue_type == 'lung'
etoposide_lung_NSCLC_lines                  <- etoposide_test$Cell_line_tissue_type == 'lung_NSCLC'
etoposide_lung_SCLC_lines                   <- etoposide_test$Cell_line_tissue_type == 'lung_SCLC'
etoposide_nervous_system_lines              <- etoposide_test$Cell_line_tissue_type == 'nervous_system'
etoposide_neuroblastoma_lines               <- etoposide_test$Cell_line_tissue_type == 'neuroblastoma'
etoposide_pancreas_lines                    <- etoposide_test$Cell_line_tissue_type == 'pancreas'
etoposide_skin_lines                        <- etoposide_test$Cell_line_tissue_type == 'skin'
etoposide_soft_tissue_lines                 <- etoposide_test$Cell_line_tissue_type == 'soft_tissue'
etoposide_thyroid_lines                     <- etoposide_test$Cell_line_tissue_type == 'thyroid'
etoposide_urogenital_system_lines           <- etoposide_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
etoposide_test_aero_dig_tract_exp <- etoposide_rna_seq_test[etoposide_aero_dig_tract_lines, ]
etoposide_test_aero_dig_tract <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_aero_dig_tract_auc <- sum(etoposide_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

etoposide_test_bone_exp <- etoposide_rna_seq_test[etoposide_bone_lines, ]
etoposide_test_bone <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_bone_auc <- sum(etoposide_test_bone$res_sens == new_ic50)/length(new_ic50)

etoposide_test_breast_exp <- etoposide_rna_seq_test[etoposide_breast_lines, ]
etoposide_test_breast <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_breast_auc <- sum(etoposide_test_breast$res_sens == new_ic50)/length(new_ic50)

etoposide_test_digestive_system_exp <- etoposide_rna_seq_test[etoposide_digestive_system_lines, ]
etoposide_test_digestive_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_digestive_system_auc <- sum(etoposide_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

etoposide_test_kidney_exp <- etoposide_rna_seq_test[etoposide_kidney_lines, ]
etoposide_test_kidney <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_kidney_auc <- sum(etoposide_test_kidney$res_sens == new_ic50)/length(new_ic50)

etoposide_test_large_intestine_exp <- etoposide_rna_seq_test[etoposide_large_intestine_lines, ]
etoposide_test_large_intestine <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_large_intestine_auc <- sum(etoposide_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

etoposide_test_lung_exp <- etoposide_rna_seq_test[etoposide_lung_lines, ]
etoposide_test_lung <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_lung_auc <- sum(etoposide_test_lung$res_sens == new_ic50)/length(new_ic50)

etoposide_test_lung_NSCLC_exp <- etoposide_rna_seq_test[etoposide_lung_NSCLC_lines, ]
etoposide_test_lung_NSCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_lung_NSCLC_auc <- sum(etoposide_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

etoposide_test_lung_SCLC_exp <- etoposide_rna_seq_test[etoposide_lung_SCLC_lines, ]
etoposide_test_lung_SCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_lung_SCLC_auc <- sum(etoposide_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

etoposide_test_nervous_system_exp <- etoposide_rna_seq_test[etoposide_nervous_system_lines, ]
etoposide_test_nervous_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_nervous_system_auc <- sum(etoposide_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

etoposide_test_neuroblastoma_exp <- etoposide_rna_seq_test[etoposide_neuroblastoma_lines, ]
etoposide_test_neuroblastoma <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_neuroblastoma_auc <- sum(etoposide_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

etoposide_test_pancreas_exp <- etoposide_rna_seq_test[etoposide_pancreas_lines, ]
etoposide_test_pancreas <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_pancreas_auc <- sum(etoposide_test_pancreas$res_sens == new_ic50)/length(new_ic50)

etoposide_test_skin_exp <- etoposide_rna_seq_test[etoposide_skin_lines, ]
etoposide_test_skin <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_skin_auc <- sum(etoposide_test_skin$res_sens == new_ic50)/length(new_ic50)

etoposide_test_soft_tissue_exp <- etoposide_rna_seq_test[etoposide_soft_tissue_lines, ]
etoposide_test_soft_tissue <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_soft_tissue_auc <- sum(etoposide_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

etoposide_test_thyroid_exp <- etoposide_rna_seq_test[etoposide_thyroid_lines, ]
etoposide_test_thyroid <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_thyroid_auc <- sum(etoposide_test_thyroid$res_sens == new_ic50)/length(new_ic50)

etoposide_test_urogenital_system_exp <- etoposide_rna_seq_test[etoposide_urogenital_system_lines, ]
etoposide_test_urogenital_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
etoposide_urogenital_system_auc <- sum(etoposide_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

etoposide_auc <- c(etoposide_aero_dig_tract_auc, etoposide_bone_auc, 
                   etoposide_breast_auc, etoposide_digestive_system_auc, 
                   etoposide_kidney_auc, etoposide_large_intestine_auc, 
                   etoposide_lung_auc, etoposide_lung_NSCLC_auc, etoposide_lung_SCLC_auc, 
                   etoposide_nervous_system_auc, etoposide_neuroblastoma_auc, 
                   etoposide_pancreas_auc, etoposide_skin_auc, 
                   etoposide_soft_tissue_auc, etoposide_thyroid_auc, 
                   etoposide_urogenital_system_auc)



gemcitabine_aero_dig_tract_lines              <- gemcitabine_test$Cell_line_tissue_type == 'aero_dig_tract'
gemcitabine_bone_lines                        <- gemcitabine_test$Cell_line_tissue_type == 'bone'
gemcitabine_breast_lines                      <- gemcitabine_test$Cell_line_tissue_type == 'breast'
gemcitabine_digestive_system_lines            <- gemcitabine_test$Cell_line_tissue_type == 'digestive_system'
gemcitabine_kidney_lines                      <- gemcitabine_test$Cell_line_tissue_type == 'kidney'
gemcitabine_large_intestine_lines             <- gemcitabine_test$Cell_line_tissue_type == 'large_intestine'
gemcitabine_lung_lines                        <- gemcitabine_test$Cell_line_tissue_type == 'lung'
gemcitabine_lung_NSCLC_lines                  <- gemcitabine_test$Cell_line_tissue_type == 'lung_NSCLC'
gemcitabine_lung_SCLC_lines                   <- gemcitabine_test$Cell_line_tissue_type == 'lung_SCLC'
gemcitabine_nervous_system_lines              <- gemcitabine_test$Cell_line_tissue_type == 'nervous_system'
gemcitabine_neuroblastoma_lines               <- gemcitabine_test$Cell_line_tissue_type == 'neuroblastoma'
gemcitabine_pancreas_lines                    <- gemcitabine_test$Cell_line_tissue_type == 'pancreas'
gemcitabine_skin_lines                        <- gemcitabine_test$Cell_line_tissue_type == 'skin'
gemcitabine_soft_tissue_lines                 <- gemcitabine_test$Cell_line_tissue_type == 'soft_tissue'
gemcitabine_thyroid_lines                     <- gemcitabine_test$Cell_line_tissue_type == 'thyroid'
gemcitabine_urogenital_system_lines           <- gemcitabine_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
gemcitabine_test_aero_dig_tract_exp <- gemcitabine_rna_seq_test[gemcitabine_aero_dig_tract_lines, ]
gemcitabine_test_aero_dig_tract <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_aero_dig_tract_auc <- sum(gemcitabine_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_bone_exp <- gemcitabine_rna_seq_test[gemcitabine_bone_lines, ]
gemcitabine_test_bone <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_bone_auc <- sum(gemcitabine_test_bone$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_breast_exp <- gemcitabine_rna_seq_test[gemcitabine_breast_lines, ]
gemcitabine_test_breast <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_breast_auc <- sum(gemcitabine_test_breast$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_digestive_system_exp <- gemcitabine_rna_seq_test[gemcitabine_digestive_system_lines, ]
gemcitabine_test_digestive_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_digestive_system_auc <- sum(gemcitabine_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_kidney_exp <- gemcitabine_rna_seq_test[gemcitabine_kidney_lines, ]
gemcitabine_test_kidney <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_kidney_auc <- sum(gemcitabine_test_kidney$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_large_intestine_exp <- gemcitabine_rna_seq_test[gemcitabine_large_intestine_lines, ]
gemcitabine_test_large_intestine <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_large_intestine_auc <- sum(gemcitabine_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_lung_exp <- gemcitabine_rna_seq_test[gemcitabine_lung_lines, ]
gemcitabine_test_lung <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_lung_auc <- sum(gemcitabine_test_lung$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_lung_NSCLC_exp <- gemcitabine_rna_seq_test[gemcitabine_lung_NSCLC_lines, ]
gemcitabine_test_lung_NSCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_lung_NSCLC_auc <- sum(gemcitabine_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_lung_SCLC_exp <- gemcitabine_rna_seq_test[gemcitabine_lung_SCLC_lines, ]
gemcitabine_test_lung_SCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_lung_SCLC_auc <- sum(gemcitabine_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_nervous_system_exp <- gemcitabine_rna_seq_test[gemcitabine_nervous_system_lines, ]
gemcitabine_test_nervous_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_nervous_system_auc <- sum(gemcitabine_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_neuroblastoma_exp <- gemcitabine_rna_seq_test[gemcitabine_neuroblastoma_lines, ]
gemcitabine_test_neuroblastoma <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_neuroblastoma_auc <- sum(gemcitabine_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_pancreas_exp <- gemcitabine_rna_seq_test[gemcitabine_pancreas_lines, ]
gemcitabine_test_pancreas <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_pancreas_auc <- sum(gemcitabine_test_pancreas$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_skin_exp <- gemcitabine_rna_seq_test[gemcitabine_skin_lines, ]
gemcitabine_test_skin <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_skin_auc <- sum(gemcitabine_test_skin$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_soft_tissue_exp <- gemcitabine_rna_seq_test[gemcitabine_soft_tissue_lines, ]
gemcitabine_test_soft_tissue <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_soft_tissue_auc <- sum(gemcitabine_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_thyroid_exp <- gemcitabine_rna_seq_test[gemcitabine_thyroid_lines, ]
gemcitabine_test_thyroid <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_thyroid_auc <- sum(gemcitabine_test_thyroid$res_sens == new_ic50)/length(new_ic50)

gemcitabine_test_urogenital_system_exp <- gemcitabine_rna_seq_test[gemcitabine_urogenital_system_lines, ]
gemcitabine_test_urogenital_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
gemcitabine_urogenital_system_auc <- sum(gemcitabine_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

gemcitabine_auc <- c(gemcitabine_aero_dig_tract_auc, gemcitabine_bone_auc, 
                   gemcitabine_breast_auc, gemcitabine_digestive_system_auc, 
                   gemcitabine_kidney_auc, gemcitabine_large_intestine_auc, 
                   gemcitabine_lung_auc, gemcitabine_lung_NSCLC_auc, gemcitabine_lung_SCLC_auc, 
                   gemcitabine_nervous_system_auc, gemcitabine_neuroblastoma_auc, 
                   gemcitabine_pancreas_auc, gemcitabine_skin_auc, 
                   gemcitabine_soft_tissue_auc, gemcitabine_thyroid_auc, 
                   gemcitabine_urogenital_system_auc)



methotrexate_aero_dig_tract_lines              <- methotrexate_test$Cell_line_tissue_type == 'aero_dig_tract'
methotrexate_bone_lines                        <- methotrexate_test$Cell_line_tissue_type == 'bone'
methotrexate_breast_lines                      <- methotrexate_test$Cell_line_tissue_type == 'breast'
methotrexate_digestive_system_lines            <- methotrexate_test$Cell_line_tissue_type == 'digestive_system'
methotrexate_kidney_lines                      <- methotrexate_test$Cell_line_tissue_type == 'kidney'
methotrexate_large_intestine_lines             <- methotrexate_test$Cell_line_tissue_type == 'large_intestine'
methotrexate_lung_lines                        <- methotrexate_test$Cell_line_tissue_type == 'lung'
methotrexate_lung_NSCLC_lines                  <- methotrexate_test$Cell_line_tissue_type == 'lung_NSCLC'
methotrexate_lung_SCLC_lines                   <- methotrexate_test$Cell_line_tissue_type == 'lung_SCLC'
methotrexate_nervous_system_lines              <- methotrexate_test$Cell_line_tissue_type == 'nervous_system'
methotrexate_neuroblastoma_lines               <- methotrexate_test$Cell_line_tissue_type == 'neuroblastoma'
methotrexate_pancreas_lines                    <- methotrexate_test$Cell_line_tissue_type == 'pancreas'
methotrexate_skin_lines                        <- methotrexate_test$Cell_line_tissue_type == 'skin'
methotrexate_soft_tissue_lines                 <- methotrexate_test$Cell_line_tissue_type == 'soft_tissue'
methotrexate_thyroid_lines                     <- methotrexate_test$Cell_line_tissue_type == 'thyroid'
methotrexate_urogenital_system_lines           <- methotrexate_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
methotrexate_test_aero_dig_tract_exp <- methotrexate_rna_seq_test[methotrexate_aero_dig_tract_lines, ]
methotrexate_test_aero_dig_tract <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_aero_dig_tract_auc <- sum(methotrexate_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_bone_exp <- methotrexate_rna_seq_test[methotrexate_bone_lines, ]
methotrexate_test_bone <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_bone_auc <- sum(methotrexate_test_bone$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_breast_exp <- methotrexate_rna_seq_test[methotrexate_breast_lines, ]
methotrexate_test_breast <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_breast_auc <- sum(methotrexate_test_breast$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_digestive_system_exp <- methotrexate_rna_seq_test[methotrexate_digestive_system_lines, ]
methotrexate_test_digestive_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_digestive_system_auc <- sum(methotrexate_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_kidney_exp <- methotrexate_rna_seq_test[methotrexate_kidney_lines, ]
methotrexate_test_kidney <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_kidney_auc <- sum(methotrexate_test_kidney$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_large_intestine_exp <- methotrexate_rna_seq_test[methotrexate_large_intestine_lines, ]
methotrexate_test_large_intestine <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_large_intestine_auc <- sum(methotrexate_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_lung_exp <- methotrexate_rna_seq_test[methotrexate_lung_lines, ]
methotrexate_test_lung <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_lung_auc <- sum(methotrexate_test_lung$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_lung_NSCLC_exp <- methotrexate_rna_seq_test[methotrexate_lung_NSCLC_lines, ]
methotrexate_test_lung_NSCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_lung_NSCLC_auc <- sum(methotrexate_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_lung_SCLC_exp <- methotrexate_rna_seq_test[methotrexate_lung_SCLC_lines, ]
methotrexate_test_lung_SCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_lung_SCLC_auc <- sum(methotrexate_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_nervous_system_exp <- methotrexate_rna_seq_test[methotrexate_nervous_system_lines, ]
methotrexate_test_nervous_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_nervous_system_auc <- sum(methotrexate_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_neuroblastoma_exp <- methotrexate_rna_seq_test[methotrexate_neuroblastoma_lines, ]
methotrexate_test_neuroblastoma <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_neuroblastoma_auc <- sum(methotrexate_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_pancreas_exp <- methotrexate_rna_seq_test[methotrexate_pancreas_lines, ]
methotrexate_test_pancreas <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_pancreas_auc <- sum(methotrexate_test_pancreas$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_skin_exp <- methotrexate_rna_seq_test[methotrexate_skin_lines, ]
methotrexate_test_skin <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_skin_auc <- sum(methotrexate_test_skin$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_soft_tissue_exp <- methotrexate_rna_seq_test[methotrexate_soft_tissue_lines, ]
methotrexate_test_soft_tissue <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_soft_tissue_auc <- sum(methotrexate_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_thyroid_exp <- methotrexate_rna_seq_test[methotrexate_thyroid_lines, ]
methotrexate_test_thyroid <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_thyroid_auc <- sum(methotrexate_test_thyroid$res_sens == new_ic50)/length(new_ic50)

methotrexate_test_urogenital_system_exp <- methotrexate_rna_seq_test[methotrexate_urogenital_system_lines, ]
methotrexate_test_urogenital_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(methotrexate_fit_elnet, newx = as.matrix(methotrexate_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
methotrexate_urogenital_system_auc <- sum(methotrexate_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

methotrexate_auc <- c(methotrexate_aero_dig_tract_auc, methotrexate_bone_auc, 
                   methotrexate_breast_auc, methotrexate_digestive_system_auc, 
                   methotrexate_kidney_auc, methotrexate_large_intestine_auc, 
                   methotrexate_lung_auc, methotrexate_lung_NSCLC_auc, methotrexate_lung_SCLC_auc, 
                   methotrexate_nervous_system_auc, methotrexate_neuroblastoma_auc, 
                   methotrexate_pancreas_auc, methotrexate_skin_auc, 
                   methotrexate_soft_tissue_auc, methotrexate_thyroid_auc, 
                   methotrexate_urogenital_system_auc)



mitomycin_aero_dig_tract_lines              <- mitomycin_test$Cell_line_tissue_type == 'aero_dig_tract'
mitomycin_bone_lines                        <- mitomycin_test$Cell_line_tissue_type == 'bone'
mitomycin_breast_lines                      <- mitomycin_test$Cell_line_tissue_type == 'breast'
mitomycin_digestive_system_lines            <- mitomycin_test$Cell_line_tissue_type == 'digestive_system'
mitomycin_kidney_lines                      <- mitomycin_test$Cell_line_tissue_type == 'kidney'
mitomycin_large_intestine_lines             <- mitomycin_test$Cell_line_tissue_type == 'large_intestine'
mitomycin_lung_lines                        <- mitomycin_test$Cell_line_tissue_type == 'lung'
mitomycin_lung_NSCLC_lines                  <- mitomycin_test$Cell_line_tissue_type == 'lung_NSCLC'
mitomycin_lung_SCLC_lines                   <- mitomycin_test$Cell_line_tissue_type == 'lung_SCLC'
mitomycin_nervous_system_lines              <- mitomycin_test$Cell_line_tissue_type == 'nervous_system'
mitomycin_neuroblastoma_lines               <- mitomycin_test$Cell_line_tissue_type == 'neuroblastoma'
mitomycin_pancreas_lines                    <- mitomycin_test$Cell_line_tissue_type == 'pancreas'
mitomycin_skin_lines                        <- mitomycin_test$Cell_line_tissue_type == 'skin'
mitomycin_soft_tissue_lines                 <- mitomycin_test$Cell_line_tissue_type == 'soft_tissue'
mitomycin_thyroid_lines                     <- mitomycin_test$Cell_line_tissue_type == 'thyroid'
mitomycin_urogenital_system_lines           <- mitomycin_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
mitomycin_test_aero_dig_tract_exp <- mitomycin_rna_seq_test[mitomycin_aero_dig_tract_lines, ]
mitomycin_test_aero_dig_tract <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_aero_dig_tract_auc <- sum(mitomycin_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_bone_exp <- mitomycin_rna_seq_test[mitomycin_bone_lines, ]
mitomycin_test_bone <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_bone_auc <- sum(mitomycin_test_bone$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_breast_exp <- mitomycin_rna_seq_test[mitomycin_breast_lines, ]
mitomycin_test_breast <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_breast_auc <- sum(mitomycin_test_breast$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_digestive_system_exp <- mitomycin_rna_seq_test[mitomycin_digestive_system_lines, ]
mitomycin_test_digestive_system <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_digestive_system_auc <- sum(mitomycin_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_kidney_exp <- mitomycin_rna_seq_test[mitomycin_kidney_lines, ]
mitomycin_test_kidney <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_kidney_auc <- sum(mitomycin_test_kidney$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_large_intestine_exp <- mitomycin_rna_seq_test[mitomycin_large_intestine_lines, ]
mitomycin_test_large_intestine <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_large_intestine_auc <- sum(mitomycin_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_lung_exp <- mitomycin_rna_seq_test[mitomycin_lung_lines, ]
mitomycin_test_lung <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_lung_auc <- sum(mitomycin_test_lung$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_lung_NSCLC_exp <- mitomycin_rna_seq_test[mitomycin_lung_NSCLC_lines, ]
mitomycin_test_lung_NSCLC <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_lung_NSCLC_auc <- sum(mitomycin_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_lung_SCLC_exp <- mitomycin_rna_seq_test[mitomycin_lung_SCLC_lines, ]
mitomycin_test_lung_SCLC <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_lung_SCLC_auc <- sum(mitomycin_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_nervous_system_exp <- mitomycin_rna_seq_test[mitomycin_nervous_system_lines, ]
mitomycin_test_nervous_system <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_nervous_system_auc <- sum(mitomycin_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_neuroblastoma_exp <- mitomycin_rna_seq_test[mitomycin_neuroblastoma_lines, ]
mitomycin_test_neuroblastoma <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_neuroblastoma_auc <- sum(mitomycin_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_pancreas_exp <- mitomycin_rna_seq_test[mitomycin_pancreas_lines, ]
mitomycin_test_pancreas <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_pancreas_auc <- sum(mitomycin_test_pancreas$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_skin_exp <- mitomycin_rna_seq_test[mitomycin_skin_lines, ]
mitomycin_test_skin <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_skin_auc <- sum(mitomycin_test_skin$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_soft_tissue_exp <- mitomycin_rna_seq_test[mitomycin_soft_tissue_lines, ]
mitomycin_test_soft_tissue <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_soft_tissue_auc <- sum(mitomycin_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_thyroid_exp <- mitomycin_rna_seq_test[mitomycin_thyroid_lines, ]
mitomycin_test_thyroid <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_thyroid_auc <- sum(mitomycin_test_thyroid$res_sens == new_ic50)/length(new_ic50)

mitomycin_test_urogenital_system_exp <- mitomycin_rna_seq_test[mitomycin_urogenital_system_lines, ]
mitomycin_test_urogenital_system <- mitomycin_test[which(mitomycin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(mitomycin_fit_elnet, newx = as.matrix(mitomycin_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
mitomycin_urogenital_system_auc <- sum(mitomycin_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

mitomycin_auc <- c(mitomycin_aero_dig_tract_auc, mitomycin_bone_auc, 
                   mitomycin_breast_auc, mitomycin_digestive_system_auc, 
                   mitomycin_kidney_auc, mitomycin_large_intestine_auc, 
                   mitomycin_lung_auc, mitomycin_lung_NSCLC_auc, mitomycin_lung_SCLC_auc, 
                   mitomycin_nervous_system_auc, mitomycin_neuroblastoma_auc, 
                   mitomycin_pancreas_auc, mitomycin_skin_auc, 
                   mitomycin_soft_tissue_auc, mitomycin_thyroid_auc, 
                   mitomycin_urogenital_system_auc)



sn38_aero_dig_tract_lines              <- sn38_test$Cell_line_tissue_type == 'aero_dig_tract'
sn38_bone_lines                        <- sn38_test$Cell_line_tissue_type == 'bone'
sn38_breast_lines                      <- sn38_test$Cell_line_tissue_type == 'breast'
sn38_digestive_system_lines            <- sn38_test$Cell_line_tissue_type == 'digestive_system'
sn38_kidney_lines                      <- sn38_test$Cell_line_tissue_type == 'kidney'
sn38_large_intestine_lines             <- sn38_test$Cell_line_tissue_type == 'large_intestine'
sn38_lung_lines                        <- sn38_test$Cell_line_tissue_type == 'lung'
sn38_lung_NSCLC_lines                  <- sn38_test$Cell_line_tissue_type == 'lung_NSCLC'
sn38_lung_SCLC_lines                   <- sn38_test$Cell_line_tissue_type == 'lung_SCLC'
sn38_nervous_system_lines              <- sn38_test$Cell_line_tissue_type == 'nervous_system'
sn38_neuroblastoma_lines               <- sn38_test$Cell_line_tissue_type == 'neuroblastoma'
sn38_pancreas_lines                    <- sn38_test$Cell_line_tissue_type == 'pancreas'
sn38_skin_lines                        <- sn38_test$Cell_line_tissue_type == 'skin'
sn38_soft_tissue_lines                 <- sn38_test$Cell_line_tissue_type == 'soft_tissue'
sn38_thyroid_lines                     <- sn38_test$Cell_line_tissue_type == 'thyroid'
sn38_urogenital_system_lines           <- sn38_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
sn38_test_aero_dig_tract_exp <- sn38_rna_seq_test[sn38_aero_dig_tract_lines, ]
sn38_test_aero_dig_tract <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_aero_dig_tract_auc <- sum(sn38_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

sn38_test_bone_exp <- sn38_rna_seq_test[sn38_bone_lines, ]
sn38_test_bone <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_bone_auc <- sum(sn38_test_bone$res_sens == new_ic50)/length(new_ic50)

sn38_test_breast_exp <- sn38_rna_seq_test[sn38_breast_lines, ]
sn38_test_breast <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_breast_auc <- sum(sn38_test_breast$res_sens == new_ic50)/length(new_ic50)

sn38_test_digestive_system_exp <- sn38_rna_seq_test[sn38_digestive_system_lines, ]
sn38_test_digestive_system <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_digestive_system_auc <- sum(sn38_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

sn38_test_kidney_exp <- sn38_rna_seq_test[sn38_kidney_lines, ]
sn38_test_kidney <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_kidney_auc <- sum(sn38_test_kidney$res_sens == new_ic50)/length(new_ic50)

sn38_test_large_intestine_exp <- sn38_rna_seq_test[sn38_large_intestine_lines, ]
sn38_test_large_intestine <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_large_intestine_auc <- sum(sn38_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

sn38_test_lung_exp <- sn38_rna_seq_test[sn38_lung_lines, ]
sn38_test_lung <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_lung_auc <- sum(sn38_test_lung$res_sens == new_ic50)/length(new_ic50)

sn38_test_lung_NSCLC_exp <- sn38_rna_seq_test[sn38_lung_NSCLC_lines, ]
sn38_test_lung_NSCLC <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_lung_NSCLC_auc <- sum(sn38_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

sn38_test_lung_SCLC_exp <- sn38_rna_seq_test[sn38_lung_SCLC_lines, ]
sn38_test_lung_SCLC <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_lung_SCLC_auc <- sum(sn38_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

sn38_test_nervous_system_exp <- sn38_rna_seq_test[sn38_nervous_system_lines, ]
sn38_test_nervous_system <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_nervous_system_auc <- sum(sn38_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

sn38_test_neuroblastoma_exp <- sn38_rna_seq_test[sn38_neuroblastoma_lines, ]
sn38_test_neuroblastoma <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_neuroblastoma_auc <- sum(sn38_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

sn38_test_pancreas_exp <- sn38_rna_seq_test[sn38_pancreas_lines, ]
sn38_test_pancreas <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_pancreas_auc <- sum(sn38_test_pancreas$res_sens == new_ic50)/length(new_ic50)

sn38_test_skin_exp <- sn38_rna_seq_test[sn38_skin_lines, ]
sn38_test_skin <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_skin_auc <- sum(sn38_test_skin$res_sens == new_ic50)/length(new_ic50)

sn38_test_soft_tissue_exp <- sn38_rna_seq_test[sn38_soft_tissue_lines, ]
sn38_test_soft_tissue <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_soft_tissue_auc <- sum(sn38_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

sn38_test_thyroid_exp <- sn38_rna_seq_test[sn38_thyroid_lines, ]
sn38_test_thyroid <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_thyroid_auc <- sum(sn38_test_thyroid$res_sens == new_ic50)/length(new_ic50)

sn38_test_urogenital_system_exp <- sn38_rna_seq_test[sn38_urogenital_system_lines, ]
sn38_test_urogenital_system <- sn38_test[which(sn38_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(sn38_fit_elnet, newx = as.matrix(sn38_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
sn38_urogenital_system_auc <- sum(sn38_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

sn38_auc <- c(sn38_aero_dig_tract_auc, sn38_bone_auc, 
                   sn38_breast_auc, sn38_digestive_system_auc, 
                   sn38_kidney_auc, sn38_large_intestine_auc, 
                   sn38_lung_auc, sn38_lung_NSCLC_auc, sn38_lung_SCLC_auc, 
                   sn38_nervous_system_auc, sn38_neuroblastoma_auc, 
                   sn38_pancreas_auc, sn38_skin_auc, 
                   sn38_soft_tissue_auc, sn38_thyroid_auc, 
                   sn38_urogenital_system_auc)



temozolomide_aero_dig_tract_lines              <- temozolomide_test$Cell_line_tissue_type == 'aero_dig_tract'
temozolomide_bone_lines                        <- temozolomide_test$Cell_line_tissue_type == 'bone'
temozolomide_breast_lines                      <- temozolomide_test$Cell_line_tissue_type == 'breast'
temozolomide_digestive_system_lines            <- temozolomide_test$Cell_line_tissue_type == 'digestive_system'
temozolomide_kidney_lines                      <- temozolomide_test$Cell_line_tissue_type == 'kidney'
temozolomide_large_intestine_lines             <- temozolomide_test$Cell_line_tissue_type == 'large_intestine'
temozolomide_lung_lines                        <- temozolomide_test$Cell_line_tissue_type == 'lung'
temozolomide_lung_NSCLC_lines                  <- temozolomide_test$Cell_line_tissue_type == 'lung_NSCLC'
temozolomide_lung_SCLC_lines                   <- temozolomide_test$Cell_line_tissue_type == 'lung_SCLC'
temozolomide_nervous_system_lines              <- temozolomide_test$Cell_line_tissue_type == 'nervous_system'
temozolomide_neuroblastoma_lines               <- temozolomide_test$Cell_line_tissue_type == 'neuroblastoma'
temozolomide_pancreas_lines                    <- temozolomide_test$Cell_line_tissue_type == 'pancreas'
temozolomide_skin_lines                        <- temozolomide_test$Cell_line_tissue_type == 'skin'
temozolomide_soft_tissue_lines                 <- temozolomide_test$Cell_line_tissue_type == 'soft_tissue'
temozolomide_thyroid_lines                     <- temozolomide_test$Cell_line_tissue_type == 'thyroid'
temozolomide_urogenital_system_lines           <- temozolomide_test$Cell_line_tissue_type == 'urogenital_system'

#test pan-cancer models against individual cancer types
temozolomide_test_aero_dig_tract_exp <- temozolomide_rna_seq_test[temozolomide_aero_dig_tract_lines, ]
temozolomide_test_aero_dig_tract <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_aero_dig_tract_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_aero_dig_tract_auc <- sum(temozolomide_test_aero_dig_tract$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_bone_exp <- temozolomide_rna_seq_test[temozolomide_bone_lines, ]
temozolomide_test_bone <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_bone_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_bone_auc <- sum(temozolomide_test_bone$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_breast_exp <- temozolomide_rna_seq_test[temozolomide_breast_lines, ]
temozolomide_test_breast <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_breast_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_breast_auc <- sum(temozolomide_test_breast$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_digestive_system_exp <- temozolomide_rna_seq_test[temozolomide_digestive_system_lines, ]
temozolomide_test_digestive_system <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_digestive_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_digestive_system_auc <- sum(temozolomide_test_digestive_system$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_kidney_exp <- temozolomide_rna_seq_test[temozolomide_kidney_lines, ]
temozolomide_test_kidney <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_kidney_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_kidney_auc <- sum(temozolomide_test_kidney$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_large_intestine_exp <- temozolomide_rna_seq_test[temozolomide_large_intestine_lines, ]
temozolomide_test_large_intestine <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_large_intestine_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_large_intestine_auc <- sum(temozolomide_test_large_intestine$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_lung_exp <- temozolomide_rna_seq_test[temozolomide_lung_lines, ]
temozolomide_test_lung <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_lung_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_lung_auc <- sum(temozolomide_test_lung$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_lung_NSCLC_exp <- temozolomide_rna_seq_test[temozolomide_lung_NSCLC_lines, ]
temozolomide_test_lung_NSCLC <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_lung_NSCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_lung_NSCLC_auc <- sum(temozolomide_test_lung_NSCLC$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_lung_SCLC_exp <- temozolomide_rna_seq_test[temozolomide_lung_SCLC_lines, ]
temozolomide_test_lung_SCLC <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_lung_SCLC_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_lung_SCLC_auc <- sum(temozolomide_test_lung_SCLC$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_nervous_system_exp <- temozolomide_rna_seq_test[temozolomide_nervous_system_lines, ]
temozolomide_test_nervous_system <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_nervous_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_nervous_system_auc <- sum(temozolomide_test_nervous_system$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_neuroblastoma_exp <- temozolomide_rna_seq_test[temozolomide_neuroblastoma_lines, ]
temozolomide_test_neuroblastoma <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_neuroblastoma_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_neuroblastoma_auc <- sum(temozolomide_test_neuroblastoma$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_pancreas_exp <- temozolomide_rna_seq_test[temozolomide_pancreas_lines, ]
temozolomide_test_pancreas <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_pancreas_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_pancreas_auc <- sum(temozolomide_test_pancreas$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_skin_exp <- temozolomide_rna_seq_test[temozolomide_skin_lines, ]
temozolomide_test_skin <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_skin_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_skin_auc <- sum(temozolomide_test_skin$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_soft_tissue_exp <- temozolomide_rna_seq_test[temozolomide_soft_tissue_lines, ]
temozolomide_test_soft_tissue <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_soft_tissue_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_soft_tissue_auc <- sum(temozolomide_test_soft_tissue$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_thyroid_exp <- temozolomide_rna_seq_test[temozolomide_thyroid_lines, ]
temozolomide_test_thyroid <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_thyroid_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_thyroid_auc <- sum(temozolomide_test_thyroid$res_sens == new_ic50)/length(new_ic50)

temozolomide_test_urogenital_system_exp <- temozolomide_rna_seq_test[temozolomide_urogenital_system_lines, ]
temozolomide_test_urogenital_system <- temozolomide_test[which(temozolomide_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(temozolomide_fit_elnet, newx = as.matrix(temozolomide_test_urogenital_system_exp), s = 'lambda.1se', interval = 'conf', type = 'class')
temozolomide_urogenital_system_auc <- sum(temozolomide_test_urogenital_system$res_sens == new_ic50)/length(new_ic50)

temozolomide_auc <- c(temozolomide_aero_dig_tract_auc, temozolomide_bone_auc, 
                   temozolomide_breast_auc, temozolomide_digestive_system_auc, 
                   temozolomide_kidney_auc, temozolomide_large_intestine_auc, 
                   temozolomide_lung_auc, temozolomide_lung_NSCLC_auc, temozolomide_lung_SCLC_auc, 
                   temozolomide_nervous_system_auc, temozolomide_neuroblastoma_auc, 
                   temozolomide_pancreas_auc, temozolomide_skin_auc, 
                   temozolomide_soft_tissue_auc, temozolomide_thyroid_auc, 
                   temozolomide_urogenital_system_auc)




#put everything together
all_gdsc_auc_by_type <- data.frame(bleomycin_auc, camptothecin_auc, cisplatin_auc, cytarabine_auc, doxorubicin_auc, 
                                   etoposide_auc, gemcitabine_auc, methotrexate_auc, mitomycin_auc, 
                                   sn38_auc, temozolomide_auc)

#all_gdsc_auc_by_type <- data.frame(t(all_gdsc_auc_by_type))
all_gdsc_auc_by_type <- rbind(all_gdsc_auc_by_type, overall_auc)
#all_gdsc_auc_by_type <- all_gdsc_auc_by_type[c(9, 1:8), ]
colnames(all_gdsc_auc_by_type) <- c('bleomycin', 'camptothecin', 'cisplatin', 'doxorubicin', 
                                    'etoposide', 'gemcitabine', 'methotrexate', 
                                    'mitomycin', 'sn38', 'temozolomide')
rownames(all_gdsc_auc_by_type) <- c('aero_dig_tract', 'bone', 'breast', 'digestive_system', 'kidney', 
                                    'large_intestine', 'lung', 'lung_NSCLC', 'lung_SCLC', 
                                    'nervous_system', 'neuroblastoma', 'pancreas', 'skin', 
                                    'soft_tissue', 'thyroid', 'urogenital_system', 'overall')

all_gdsc_auc_by_type <- round(all_gdsc_auc_by_type, digits = 2)
all_gdsc_auc_by_type <- data.frame(t(all_gdsc_auc_by_type))
all_gdsc_auc_by_type <- all_gdsc_auc_by_type[, c(17, 1:16)]
colors <- colorRampPalette(c("dodgerblue", "white", "red"))(100)

png(filename = 'Images/GDSC_AUC_heatmap.png', width = 960)
heatmap.2(as.matrix(all_gdsc_auc_by_type), trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, density.info = 'none', key.xlab = 'AUC', key.title = '', cexRow = 0.9, cexCol = 0.9, cellnote = all_gdsc_auc_by_type, notecol = 'black', colsep = 1, sepwidth = c(0.1,0.1), srtCol = 45, margins = c(5,8))
dev.off()


### testing on TCGA classes treated with drug ----
# BLCA W CISPLATIN (49)
blca_clinical <- read.csv('Processed_Clinical_Data/blca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(blca_clinical$most_sensitive)
blca_clinical <- blca_clinical[!na_idx, ]
table(blca_clinical$drug_name)
blca_clinical_cisplatin <- blca_clinical[which(blca_clinical$drug_name == 'cisplatin' | blca_clinical$drug_name == 'Cisplatin' | 
                                                 blca_clinical$drug_name == 'Cisplatnin'), ]

blca_gene <- read.csv('Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(blca_gene) <- gsub('\\.', '-', colnames(blca_gene))
blca_matching_idx <- blca_clinical_cisplatin$submitter_id.samples %in% colnames(blca_gene)
blca_clinical_cisplatin_short <- blca_clinical_cisplatin[blca_matching_idx, ]
blca_matching_idx <- colnames(blca_gene) %in% blca_clinical_cisplatin_short$submitter_id.samples
blca_gene_short <- blca_gene[, blca_matching_idx]
blca_gene_short <- t(blca_gene_short)
blca_gene_short_scaled <- apply(blca_gene_short, 2, scale)

new_blca_tcga_cisplatin <- predict(cisplatin_fit_elnet, newx = as.matrix(blca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

blca_surv_times <- blca_clinical_cisplatin_short$PFS
blca_status <- ifelse(blca_clinical_cisplatin_short$PFS == blca_clinical_cisplatin_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_cisplatin)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_cisplatin,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_cisplatin,
                data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)



# BLCA W GEMCITABINE (30)
blca_clinical <- read.csv('Processed_Clinical_Data/blca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(blca_clinical$most_sensitive)
blca_clinical <- blca_clinical[!na_idx, ]
table(blca_clinical$drug_name)
blca_clinical_gemcitabine <- blca_clinical[which(blca_clinical$drug_name == 'gemcitabine' | blca_clinical$drug_name == 'Gemcitabine' | 
                                                 blca_clinical$drug_name == 'Gemzar'), ]

blca_gene <- read.csv('Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(blca_gene) <- gsub('\\.', '-', colnames(blca_gene))
blca_matching_idx <- blca_clinical_gemcitabine$submitter_id.samples %in% colnames(blca_gene)
blca_clinical_gemcitabine_short <- blca_clinical_gemcitabine[blca_matching_idx, ]
blca_matching_idx <- colnames(blca_gene) %in% blca_clinical_gemcitabine_short$submitter_id.samples
blca_gene_short <- blca_gene[, blca_matching_idx]
blca_gene_short <- t(blca_gene_short)
blca_gene_short_scaled <- apply(blca_gene_short, 2, scale)

new_blca_tcga_gemcitabine <- predict(gemcitabine_fit_elnet, newx = as.matrix(blca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

blca_surv_times <- blca_clinical_gemcitabine_short$PFS
blca_status <- ifelse(blca_clinical_gemcitabine_short$PFS == blca_clinical_gemcitabine_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_gemcitabine)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
                data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)


# BRCA W DOXORUBICIN (47)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_doxorubicin <- brca_clinical[which(brca_clinical$drug_name == 'doxorubicin' | brca_clinical$drug_name == 'Doxorubicin' |
                                                 brca_clinical$drug_name == 'DOXORUBICIN'), ]

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_doxorubicin$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_doxorubicin_short <- brca_clinical_doxorubicin[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_doxorubicin_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

new_brca_tcga_doxorubicin <- predict(doxorubicin_fit_elnet, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

brca_surv_times <- brca_clinical_doxorubicin_short$PFS
brca_status <- ifelse(brca_clinical_doxorubicin_short$PFS == brca_clinical_doxorubicin_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_doxorubicin)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_doxorubicin,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_doxorubicin,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)

# BRCA W GEMCITABINE (4)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_gemcitabine <- brca_clinical[which(brca_clinical$drug_name == 'Gemcitabine' | brca_clinical$drug_name == 'GEMZAR'), ]

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_gemcitabine$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_gemcitabine_short <- brca_clinical_gemcitabine[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_gemcitabine_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)
new_blca_tcga_gemcitabine <- predict(gemcitabine_fit_elnet, newx = as.matrix(blca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

new_brca_tcga_gemcitabine <- predict(gemcitabine_fit_elnet, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

brca_surv_times <- brca_clinical_gemcitabine_short$PFS
brca_status <- ifelse(brca_clinical_gemcitabine_short$PFS == brca_clinical_gemcitabine_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_gemcitabine)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_gemcitabine,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_gemcitabine,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)


# BRCA W METHOTREXATE (9)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_methotrexate <- brca_clinical[which(brca_clinical$drug_name == 'Methotrexate' | brca_clinical$drug_name == 'METHOTREXATE'), ]

brca_clinical_methotrexate$most_sensitive  <- ifelse(brca_clinical_methotrexate$PFS < quantile(brca_clinical_methotrexate$PFS, probs = .20), 1, 0)
brca_clinical_methotrexate$least_sensitive <- ifelse(brca_clinical_methotrexate$PFS > quantile(brca_clinical_methotrexate$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_methotrexate$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_methotrexate_short <- brca_clinical_methotrexate[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_methotrexate_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

new_brca_tcga_methotrexate <- predict(methotrexate_fit_elnet, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

brca_surv_times <- brca_clinical_methotrexate_short$PFS
brca_status <- ifelse(brca_clinical_methotrexate_short$PFS == brca_clinical_methotrexate_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_methotrexate)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_methotrexate,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_methotrexate,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)



# CESC W CISPLATIN (104)
cesc_clinical <- read.csv('Processed_Clinical_Data/cesc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(cesc_clinical$most_sensitive)
cesc_clinical <- cesc_clinical[!na_idx, ]
table(cesc_clinical$drug_name)
cesc_clinical_cisplatin <- cesc_clinical[which(cesc_clinical$drug_name == 'cisplatin' | cesc_clinical$drug_name == 'Cisplatin' | 
                                                 cesc_clinical$drug_name == 'Cisplatin-xrt' | cesc_clinical$drug_name == 'Cisplatinum'), ]

cesc_clinical_cisplatin$most_sensitive  <- ifelse(cesc_clinical_cisplatin$PFS < quantile(cesc_clinical_cisplatin$PFS, probs = .20), 1, 0)
cesc_clinical_cisplatin$least_sensitive <- ifelse(cesc_clinical_cisplatin$PFS > quantile(cesc_clinical_cisplatin$PFS, probs = .80), 1, 0)

cesc_gene <- read.csv('Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(cesc_gene) <- gsub('\\.', '-', colnames(cesc_gene))
cesc_matching_idx <- cesc_clinical_cisplatin$submitter_id.samples %in% colnames(cesc_gene)
cesc_clinical_cisplatin_short <- cesc_clinical_cisplatin[cesc_matching_idx, ]
cesc_matching_idx <- colnames(cesc_gene) %in% cesc_clinical_cisplatin_short$submitter_id.samples
cesc_gene_short <- cesc_gene[, cesc_matching_idx]
cesc_gene_short <- t(cesc_gene_short)
cesc_gene_short_scaled <- apply(cesc_gene_short, 2, scale)

new_cesc_tcga_cisplatin <- predict(cisplatin_fit_elnet, newx = as.matrix(cesc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

cesc_surv_times <- cesc_clinical_cisplatin_short$PFS
cesc_status <- ifelse(cesc_clinical_cisplatin_short$PFS == cesc_clinical_cisplatin_short$OS, 0, 1)

cesc_surv_df <- data.frame(cesc_surv_times, cesc_status, new_cesc_tcga_cisplatin)
fit <- survfit(Surv(cesc_surv_times, cesc_status) ~ new_cesc_tcga_cisplatin,
               data = cesc_surv_df)
fit2 <- survfit(Surv(cesc_surv_times, cesc_status) ~ new_cesc_tcga_cisplatin,
                data = cesc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)

# CHOL W GEMCITABINE (8)
chol_clinical <- read.csv('Processed_Clinical_Data/chol_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(chol_clinical$most_sensitive)
chol_clinical <- chol_clinical[!na_idx, ]
table(chol_clinical$drug_name)
chol_clinical_gemcitabine <- chol_clinical[which(chol_clinical$drug_name == 'gemcitabine' | chol_clinical$drug_name == 'Gemcitabine' | 
                                                 chol_clinical$drug_name == 'GEMCITABINE' | chol_clinical$drug_name == 'Gemzar'), ]

chol_clinical_gemcitabine$most_sensitive  <- ifelse(chol_clinical_gemcitabine$PFS < quantile(chol_clinical_gemcitabine$PFS, probs = .20), 1, 0)
chol_clinical_gemcitabine$least_sensitive <- ifelse(chol_clinical_gemcitabine$PFS > quantile(chol_clinical_gemcitabine$PFS, probs = .80), 1, 0)

chol_gene <- read.csv('Processed_Gene_Expression/chol_tcga_rna_seq_processed.csv', row.names = 1)
colnames(chol_gene) <- gsub('\\.', '-', colnames(chol_gene))
chol_matching_idx <- chol_clinical_gemcitabine$submitter_id.samples %in% colnames(chol_gene)
chol_clinical_gemcitabine_short <- chol_clinical_gemcitabine[chol_matching_idx, ]
chol_matching_idx <- colnames(chol_gene) %in% chol_clinical_gemcitabine_short$submitter_id.samples
chol_gene_short <- chol_gene[, chol_matching_idx]
chol_gene_short <- t(chol_gene_short)
chol_gene_short_scaled <- apply(chol_gene_short, 2, scale)

rm(chol_tcga_most_min_auc)
rm(chol_tcga_most_1se_auc)
rm(chol_tcga_least_min_auc)
rm(chol_tcga_least_1se_auc)

new_chol_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
chol_tcga_most_min_auc <- auc(chol_clinical_gemcitabine_short$most_sensitive, new_chol_tcga_gemcitabine_most_sensitive_min)
chol_tcga_most_min_auc <- round(chol_tcga_most_min_auc, digits = 2)

new_chol_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
chol_tcga_most_1se_auc <- auc(chol_clinical_gemcitabine_short$most_sensitive, new_chol_tcga_gemcitabine_most_sensitive_1se)
chol_tcga_most_1se_auc <- round(chol_tcga_most_1se_auc, digits = 2)

new_chol_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
chol_tcga_least_min_auc <- auc(chol_clinical_gemcitabine_short$least_sensitive, new_chol_tcga_gemcitabine_least_sensitive_min)
chol_tcga_least_min_auc <- round(chol_tcga_least_min_auc, digits = 2)

new_chol_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
chol_tcga_least_1se_auc <- auc(chol_clinical_gemcitabine_short$least_sensitive, new_chol_tcga_gemcitabine_least_sensitive_1se)
chol_tcga_least_1se_auc <- round(chol_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_chol_tcga_gemcitabine_most_sensitive_min, chol_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_chol_tcga_gemcitabine_most_sensitive_1se, chol_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_chol_tcga_gemcitabine_least_sensitive_min, chol_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_chol_tcga_gemcitabine_least_sensitive_1se, chol_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/chol_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# ESCA W CISPLATIN (15)
esca_clinical <- read.csv('Processed_Clinical_Data/esca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(esca_clinical$most_sensitive)
esca_clinical <- esca_clinical[!na_idx, ]
table(esca_clinical$drug_name)
esca_clinical_cisplatin <- esca_clinical[which(esca_clinical$drug_name == 'Cisplatin'), ]

esca_clinical_cisplatin$most_sensitive  <- ifelse(esca_clinical_cisplatin$PFS < quantile(esca_clinical_cisplatin$PFS, probs = .20), 1, 0)
esca_clinical_cisplatin$least_sensitive <- ifelse(esca_clinical_cisplatin$PFS > quantile(esca_clinical_cisplatin$PFS, probs = .80), 1, 0)

esca_gene <- read.csv('Processed_Gene_Expression/esca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(esca_gene) <- gsub('\\.', '-', colnames(esca_gene))
esca_matching_idx <- esca_clinical_cisplatin$submitter_id.samples %in% colnames(esca_gene)
esca_clinical_cisplatin_short <- esca_clinical_cisplatin[esca_matching_idx, ]
esca_matching_idx <- colnames(esca_gene) %in% esca_clinical_cisplatin_short$submitter_id.samples
esca_gene_short <- esca_gene[, esca_matching_idx]
esca_gene_short <- t(esca_gene_short)
esca_gene_short_scaled <- apply(esca_gene_short, 2, scale)

rm(esca_tcga_most_min_auc)
rm(esca_tcga_most_1se_auc)
rm(esca_tcga_least_min_auc)
rm(esca_tcga_least_1se_auc)

new_esca_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
esca_tcga_most_min_auc <- auc(esca_clinical_cisplatin_short$most_sensitive, new_esca_tcga_cisplatin_most_sensitive_min)
esca_tcga_most_min_auc <- round(esca_tcga_most_min_auc, digits = 2)

new_esca_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
esca_tcga_most_1se_auc <- auc(esca_clinical_cisplatin_short$most_sensitive, new_esca_tcga_cisplatin_most_sensitive_1se)
esca_tcga_most_1se_auc <- round(esca_tcga_most_1se_auc, digits = 2)

new_esca_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
esca_tcga_least_min_auc <- auc(esca_clinical_cisplatin_short$least_sensitive, new_esca_tcga_cisplatin_least_sensitive_min)
esca_tcga_least_min_auc <- round(esca_tcga_least_min_auc, digits = 2)

new_esca_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
esca_tcga_least_1se_auc <- auc(esca_clinical_cisplatin_short$least_sensitive, new_esca_tcga_cisplatin_least_sensitive_1se)
esca_tcga_least_1se_auc <- round(esca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_esca_tcga_cisplatin_most_sensitive_min, esca_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_esca_tcga_cisplatin_most_sensitive_1se, esca_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_esca_tcga_cisplatin_least_sensitive_min, esca_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_esca_tcga_cisplatin_least_sensitive_1se, esca_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/esca_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[2], lwd = 2, lty = 2, main = 'ESCA treated with cisplatin')
plot(perf3, col = colors_i_need[2], lwd = 2, add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.50)', 'resistant model (AUC = 0.67)'), lwd = 2, lty = c(1,2), col = colors_i_need[2], cex = 0.8, bty = 'n')
dev.off()

# HNSC W CISPLATIN (71)
hnsc_clinical <- read.csv('Processed_Clinical_Data/hnsc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(hnsc_clinical$most_sensitive)
hnsc_clinical <- hnsc_clinical[!na_idx, ]
table(hnsc_clinical$drug_name)
hnsc_clinical_cisplatin <- hnsc_clinical[which(hnsc_clinical$drug_name == 'cisplatin' | hnsc_clinical$drug_name == 'Cisplatin' | 
                                                 hnsc_clinical$drug_name == 'CISPLATIN'), ]

hnsc_clinical_cisplatin$most_sensitive  <- ifelse(hnsc_clinical_cisplatin$PFS < quantile(hnsc_clinical_cisplatin$PFS, probs = .20), 1, 0)
hnsc_clinical_cisplatin$least_sensitive <- ifelse(hnsc_clinical_cisplatin$PFS > quantile(hnsc_clinical_cisplatin$PFS, probs = .80), 1, 0)

hnsc_gene <- read.csv('Processed_Gene_Expression/hnsc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(hnsc_gene) <- gsub('\\.', '-', colnames(hnsc_gene))
hnsc_matching_idx <- hnsc_clinical_cisplatin$submitter_id.samples %in% colnames(hnsc_gene)
hnsc_clinical_cisplatin_short <- hnsc_clinical_cisplatin[hnsc_matching_idx, ]
hnsc_matching_idx <- colnames(hnsc_gene) %in% hnsc_clinical_cisplatin_short$submitter_id.samples
hnsc_gene_short <- hnsc_gene[, hnsc_matching_idx]
hnsc_gene_short <- t(hnsc_gene_short)
hnsc_gene_short_scaled <- apply(hnsc_gene_short, 2, scale)

rm(hnsc_tcga_most_min_auc)
rm(hnsc_tcga_most_1se_auc)
rm(hnsc_tcga_least_min_auc)
rm(hnsc_tcga_least_1se_auc)

new_hnsc_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_min_auc <- auc(hnsc_clinical_cisplatin_short$most_sensitive, new_hnsc_tcga_cisplatin_most_sensitive_min)
hnsc_tcga_most_min_auc <- round(hnsc_tcga_most_min_auc, digits = 2)

new_hnsc_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_1se_auc <- auc(hnsc_clinical_cisplatin_short$most_sensitive, new_hnsc_tcga_cisplatin_most_sensitive_1se)
hnsc_tcga_most_1se_auc <- round(hnsc_tcga_most_1se_auc, digits = 2)

new_hnsc_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_min_auc <- auc(hnsc_clinical_cisplatin_short$least_sensitive, new_hnsc_tcga_cisplatin_least_sensitive_min)
hnsc_tcga_least_min_auc <- round(hnsc_tcga_least_min_auc, digits = 2)

new_hnsc_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_1se_auc <- auc(hnsc_clinical_cisplatin_short$least_sensitive, new_hnsc_tcga_cisplatin_least_sensitive_1se)
hnsc_tcga_least_1se_auc <- round(hnsc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_hnsc_tcga_cisplatin_most_sensitive_min, hnsc_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_hnsc_tcga_cisplatin_most_sensitive_1se, hnsc_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_hnsc_tcga_cisplatin_least_sensitive_min, hnsc_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_hnsc_tcga_cisplatin_least_sensitive_1se, hnsc_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/hnsc_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[2], lwd = 2, lty = 2, main = 'HNSC treated with cisplatin')
plot(perf3, col = colors_i_need[2], lwd = 2, add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.53)', 'resistant model (AUC = 0.61)'), lty = c(1,2), col = colors_i_need[2], cex = 0.8, bty = 'n', lwd = 2)
dev.off()

# KIRC W GEMCITABINE (6)
kirc_clinical <- read.csv('Processed_Clinical_Data/kirc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(kirc_clinical$most_sensitive)
kirc_clinical <- kirc_clinical[!na_idx, ]
table(kirc_clinical$drug_name)
kirc_clinical_gemcitabine <- kirc_clinical[which(kirc_clinical$drug_name == 'gemcitabine' | kirc_clinical$drug_name == 'Gemcitabine'), ]

kirc_clinical_gemcitabine$most_sensitive  <- ifelse(kirc_clinical_gemcitabine$PFS < quantile(kirc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
kirc_clinical_gemcitabine$least_sensitive <- ifelse(kirc_clinical_gemcitabine$PFS > quantile(kirc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

kirc_gene <- read.csv('Processed_Gene_Expression/kirc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(kirc_gene) <- gsub('\\.', '-', colnames(kirc_gene))
kirc_matching_idx <- kirc_clinical_gemcitabine$submitter_id.samples %in% colnames(kirc_gene)
kirc_clinical_gemcitabine_short <- kirc_clinical_gemcitabine[kirc_matching_idx, ]
kirc_matching_idx <- colnames(kirc_gene) %in% kirc_clinical_gemcitabine_short$submitter_id.samples
kirc_gene_short <- kirc_gene[, kirc_matching_idx]
kirc_gene_short <- t(kirc_gene_short)
kirc_gene_short_scaled <- apply(kirc_gene_short, 2, scale)

rm(kirc_tcga_most_min_auc)
rm(kirc_tcga_most_1se_auc)
rm(kirc_tcga_least_min_auc)
rm(kirc_tcga_least_1se_auc)

new_kirc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = kirc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
kirc_tcga_most_min_auc <- auc(kirc_clinical_gemcitabine_short$most_sensitive, new_kirc_tcga_gemcitabine_most_sensitive_min)
kirc_tcga_most_min_auc <- round(kirc_tcga_most_min_auc, digits = 2)

new_kirc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = kirc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
kirc_tcga_most_1se_auc <- auc(kirc_clinical_gemcitabine_short$most_sensitive, new_kirc_tcga_gemcitabine_most_sensitive_1se)
kirc_tcga_most_1se_auc <- round(kirc_tcga_most_1se_auc, digits = 2)

new_kirc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = kirc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
kirc_tcga_least_min_auc <- auc(kirc_clinical_gemcitabine_short$least_sensitive, new_kirc_tcga_gemcitabine_least_sensitive_min)
kirc_tcga_least_min_auc <- round(kirc_tcga_least_min_auc, digits = 2)

new_kirc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = kirc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
kirc_tcga_least_1se_auc <- auc(kirc_clinical_gemcitabine_short$least_sensitive, new_kirc_tcga_gemcitabine_least_sensitive_1se)
kirc_tcga_least_1se_auc <- round(kirc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_kirc_tcga_gemcitabine_most_sensitive_min, kirc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_kirc_tcga_gemcitabine_most_sensitive_1se, kirc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_kirc_tcga_gemcitabine_least_sensitive_min, kirc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_kirc_tcga_gemcitabine_least_sensitive_1se, kirc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/kirc_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LIHC W GEMCITABINE (8)
lihc_clinical <- read.csv('Processed_Clinical_Data/lihc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lihc_clinical$most_sensitive)
lihc_clinical <- lihc_clinical[!na_idx, ]
table(lihc_clinical$drug_name)
lihc_clinical_gemcitabine <- lihc_clinical[which(lihc_clinical$drug_name == 'Gemcitabine'), ]

lihc_clinical_gemcitabine$most_sensitive  <- ifelse(lihc_clinical_gemcitabine$PFS < quantile(lihc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
lihc_clinical_gemcitabine$least_sensitive <- ifelse(lihc_clinical_gemcitabine$PFS > quantile(lihc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

lihc_gene <- read.csv('Processed_Gene_Expression/lihc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lihc_gene) <- gsub('\\.', '-', colnames(lihc_gene))
lihc_matching_idx <- lihc_clinical_gemcitabine$submitter_id.samples %in% colnames(lihc_gene)
lihc_clinical_gemcitabine_short <- lihc_clinical_gemcitabine[lihc_matching_idx, ]
lihc_matching_idx <- colnames(lihc_gene) %in% lihc_clinical_gemcitabine_short$submitter_id.samples
lihc_gene_short <- lihc_gene[, lihc_matching_idx]
lihc_gene_short <- t(lihc_gene_short)
lihc_gene_short_scaled <- apply(lihc_gene_short, 2, scale)

rm(lihc_tcga_most_min_auc)
rm(lihc_tcga_most_1se_auc)
rm(lihc_tcga_least_min_auc)
rm(lihc_tcga_least_1se_auc)

new_lihc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_most_min_auc <- auc(lihc_clinical_gemcitabine_short$most_sensitive, new_lihc_tcga_gemcitabine_most_sensitive_min)
lihc_tcga_most_min_auc <- round(lihc_tcga_most_min_auc, digits = 2)

new_lihc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_most_1se_auc <- auc(lihc_clinical_gemcitabine_short$most_sensitive, new_lihc_tcga_gemcitabine_most_sensitive_1se)
lihc_tcga_most_1se_auc <- round(lihc_tcga_most_1se_auc, digits = 2)

new_lihc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_least_min_auc <- auc(lihc_clinical_gemcitabine_short$least_sensitive, new_lihc_tcga_gemcitabine_least_sensitive_min)
lihc_tcga_least_min_auc <- round(lihc_tcga_least_min_auc, digits = 2)

new_lihc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_least_1se_auc <- auc(lihc_clinical_gemcitabine_short$least_sensitive, new_lihc_tcga_gemcitabine_least_sensitive_1se)
lihc_tcga_least_1se_auc <- round(lihc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lihc_tcga_gemcitabine_most_sensitive_min, lihc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_lihc_tcga_gemcitabine_most_sensitive_1se, lihc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_lihc_tcga_gemcitabine_least_sensitive_min, lihc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_lihc_tcga_gemcitabine_least_sensitive_1se, lihc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lihc_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LIHC W TEMOZOLOMIDE (4)
lihc_clinical <- read.csv('Processed_Clinical_Data/lihc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lihc_clinical$most_sensitive)
lihc_clinical <- lihc_clinical[!na_idx, ]
table(lihc_clinical$drug_name)
lihc_clinical_temozolomide <- lihc_clinical[which(lihc_clinical$drug_name == 'Temozolomide'), ]

lihc_clinical_temozolomide$most_sensitive  <- ifelse(lihc_clinical_temozolomide$PFS < quantile(lihc_clinical_temozolomide$PFS, probs = .20), 1, 0)
lihc_clinical_temozolomide$least_sensitive <- ifelse(lihc_clinical_temozolomide$PFS > quantile(lihc_clinical_temozolomide$PFS, probs = .80), 1, 0)

lihc_gene <- read.csv('Processed_Gene_Expression/lihc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lihc_gene) <- gsub('\\.', '-', colnames(lihc_gene))
lihc_matching_idx <- lihc_clinical_temozolomide$submitter_id.samples %in% colnames(lihc_gene)
lihc_clinical_temozolomide_short <- lihc_clinical_temozolomide[lihc_matching_idx, ]
lihc_matching_idx <- colnames(lihc_gene) %in% lihc_clinical_temozolomide_short$submitter_id.samples
lihc_gene_short <- lihc_gene[, lihc_matching_idx]
lihc_gene_short <- t(lihc_gene_short)
lihc_gene_short_scaled <- apply(lihc_gene_short, 2, scale)

rm(lihc_tcga_most_min_auc)
rm(lihc_tcga_most_1se_auc)
rm(lihc_tcga_least_min_auc)
rm(lihc_tcga_least_1se_auc)

new_lihc_tcga_temozolomide_most_sensitive_min <- predict(temozolomide_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_most_min_auc <- auc(lihc_clinical_temozolomide_short$most_sensitive, new_lihc_tcga_temozolomide_most_sensitive_min)
lihc_tcga_most_min_auc <- round(lihc_tcga_most_min_auc, digits = 2)

new_lihc_tcga_temozolomide_most_sensitive_1se <- predict(temozolomide_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_most_1se_auc <- auc(lihc_clinical_temozolomide_short$most_sensitive, new_lihc_tcga_temozolomide_most_sensitive_1se)
lihc_tcga_most_1se_auc <- round(lihc_tcga_most_1se_auc, digits = 2)

new_lihc_tcga_temozolomide_least_sensitive_min <- predict(temozolomide_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_least_min_auc <- auc(lihc_clinical_temozolomide_short$least_sensitive, new_lihc_tcga_temozolomide_least_sensitive_min)
lihc_tcga_least_min_auc <- round(lihc_tcga_least_min_auc, digits = 2)

new_lihc_tcga_temozolomide_least_sensitive_1se <- predict(temozolomide_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_least_1se_auc <- auc(lihc_clinical_temozolomide_short$least_sensitive, new_lihc_tcga_temozolomide_least_sensitive_1se)
lihc_tcga_least_1se_auc <- round(lihc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lihc_tcga_temozolomide_most_sensitive_min, lihc_clinical_temozolomide_short$most_sensitive)
pred2 <- prediction(new_lihc_tcga_temozolomide_most_sensitive_1se, lihc_clinical_temozolomide_short$most_sensitive)
pred3 <- prediction(new_lihc_tcga_temozolomide_least_sensitive_min, lihc_clinical_temozolomide_short$least_sensitive)
pred4 <- prediction(new_lihc_tcga_temozolomide_least_sensitive_1se, lihc_clinical_temozolomide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lihc_temozolomide_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUAD W CISPLATIN (30)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_cisplatin <- luad_clinical[which(luad_clinical$drug_name == 'cisplatin' | luad_clinical$drug_name == 'Cisplatin' | 
                                                 luad_clinical$drug_name == 'CISPLATIN' | luad_clinical$drug_name == 'CISplatinum'), ]

luad_clinical_cisplatin$most_sensitive  <- ifelse(luad_clinical_cisplatin$PFS < quantile(luad_clinical_cisplatin$PFS, probs = .20), 1, 0)
luad_clinical_cisplatin$least_sensitive <- ifelse(luad_clinical_cisplatin$PFS > quantile(luad_clinical_cisplatin$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_cisplatin$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_cisplatin_short <- luad_clinical_cisplatin[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_cisplatin_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_cisplatin_short$most_sensitive, new_luad_tcga_cisplatin_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_cisplatin_short$most_sensitive, new_luad_tcga_cisplatin_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_cisplatin_short$least_sensitive, new_luad_tcga_cisplatin_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_cisplatin_short$least_sensitive, new_luad_tcga_cisplatin_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_cisplatin_most_sensitive_min, luad_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_cisplatin_most_sensitive_1se, luad_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_cisplatin_least_sensitive_min, luad_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_cisplatin_least_sensitive_1se, luad_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 1, main = 'LUAD treated with cisplatin')
plot(perf3, col = 'blue', add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.56)', 'least_sensitive (AUC = 0.86)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
dev.off()

# LUAD W ETOPOSIDE (7)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_etoposide <- luad_clinical[which(luad_clinical$drug_name == 'Etoposide'), ]

luad_clinical_etoposide$most_sensitive  <- ifelse(luad_clinical_etoposide$PFS < quantile(luad_clinical_etoposide$PFS, probs = .20), 1, 0)
luad_clinical_etoposide$least_sensitive <- ifelse(luad_clinical_etoposide$PFS > quantile(luad_clinical_etoposide$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_etoposide$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_etoposide_short <- luad_clinical_etoposide[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_etoposide_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)
for(i in 1:nrow(luad_gene_short_scaled)) {
  for (j in 1:ncol(luad_gene_short_scaled)) {
    if (is.na(luad_gene_short_scaled[i,j]) == TRUE) {
      luad_gene_short_scaled[i,j] <- 0
    }
  }
}

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_etoposide_short$most_sensitive, new_luad_tcga_etoposide_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_etoposide_most_sensitive_1se <- predict(etoposide_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_etoposide_short$most_sensitive, new_luad_tcga_etoposide_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_etoposide_short$least_sensitive, new_luad_tcga_etoposide_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_etoposide_least_sensitive_1se <- predict(etoposide_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_etoposide_short$least_sensitive, new_luad_tcga_etoposide_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_etoposide_most_sensitive_min, luad_clinical_etoposide_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_etoposide_most_sensitive_1se, luad_clinical_etoposide_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_etoposide_least_sensitive_min, luad_clinical_etoposide_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_etoposide_least_sensitive_1se, luad_clinical_etoposide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_etoposide_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = colors_i_need[5], lwd = 2, lty = 2, main = 'LUAD treated with etoposide (GDSC)')
#plot(perf2, col = 'red', lty = 1, main = 'LUAD treated with etoposide (GDSC)')
plot(perf3, col = colors_i_need[5], lwd = 2, lty =1, add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.5, y = 0.3, legend = c('sensitive model (AUC = 1.0)', 'resistant model (AUC = 0.50'), lwd = 2, cex = 0.8, lty = c(1,2), col = colors_i_need[5], bty = 'n')
dev.off()

# LUAD W GEMCITABINE (18)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_gemcitabine <- luad_clinical[which(luad_clinical$drug_name == 'gemcitabine' | luad_clinical$drug_name == 'Gemcitabine' | 
                                                 luad_clinical$drug_name == 'Gemzar'), ]

luad_clinical_gemcitabine$most_sensitive  <- ifelse(luad_clinical_gemcitabine$PFS < quantile(luad_clinical_gemcitabine$PFS, probs = .20), 1, 0)
luad_clinical_gemcitabine$least_sensitive <- ifelse(luad_clinical_gemcitabine$PFS > quantile(luad_clinical_gemcitabine$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_gemcitabine$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_gemcitabine_short <- luad_clinical_gemcitabine[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_gemcitabine_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_gemcitabine_short$most_sensitive, new_luad_tcga_gemcitabine_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_gemcitabine_short$most_sensitive, new_luad_tcga_gemcitabine_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_gemcitabine_short$least_sensitive, new_luad_tcga_gemcitabine_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_gemcitabine_short$least_sensitive, new_luad_tcga_gemcitabine_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_gemcitabine_most_sensitive_min, luad_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_gemcitabine_most_sensitive_1se, luad_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_gemcitabine_least_sensitive_min, luad_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_gemcitabine_least_sensitive_1se, luad_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUSC W CISPLATIN (24)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_cisplatin <- lusc_clinical[which(lusc_clinical$drug_name == 'cisplatin' | lusc_clinical$drug_name == 'Cisplatin'), ]

lusc_clinical_cisplatin$most_sensitive  <- ifelse(lusc_clinical_cisplatin$PFS < quantile(lusc_clinical_cisplatin$PFS, probs = .20), 1, 0)
lusc_clinical_cisplatin$least_sensitive <- ifelse(lusc_clinical_cisplatin$PFS > quantile(lusc_clinical_cisplatin$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_cisplatin$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_cisplatin_short <- lusc_clinical_cisplatin[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_cisplatin_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_cisplatin_short$most_sensitive, new_lusc_tcga_cisplatin_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_cisplatin_short$most_sensitive, new_lusc_tcga_cisplatin_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_cisplatin_short$least_sensitive, new_lusc_tcga_cisplatin_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_cisplatin_short$least_sensitive, new_lusc_tcga_cisplatin_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_cisplatin_most_sensitive_min, lusc_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_cisplatin_most_sensitive_1se, lusc_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_cisplatin_least_sensitive_min, lusc_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_cisplatin_least_sensitive_1se, lusc_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUSC W ETOPOSIDE (6)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_etoposide <- lusc_clinical[which(lusc_clinical$drug_name == 'etoposide' | lusc_clinical$drug_name == 'Etoposide'), ]

lusc_clinical_etoposide$most_sensitive  <- ifelse(lusc_clinical_etoposide$PFS < quantile(lusc_clinical_etoposide$PFS, probs = .20), 1, 0)
lusc_clinical_etoposide$least_sensitive <- ifelse(lusc_clinical_etoposide$PFS > quantile(lusc_clinical_etoposide$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_etoposide$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_etoposide_short <- lusc_clinical_etoposide[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_etoposide_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_etoposide_short$most_sensitive, new_lusc_tcga_etoposide_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_etoposide_most_sensitive_1se <- predict(etoposide_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_etoposide_short$most_sensitive, new_lusc_tcga_etoposide_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_etoposide_short$least_sensitive, new_lusc_tcga_etoposide_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_etoposide_least_sensitive_1se <- predict(etoposide_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_etoposide_short$least_sensitive, new_lusc_tcga_etoposide_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_etoposide_most_sensitive_min, lusc_clinical_etoposide_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_etoposide_most_sensitive_1se, lusc_clinical_etoposide_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_etoposide_least_sensitive_min, lusc_clinical_etoposide_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_etoposide_least_sensitive_1se, lusc_clinical_etoposide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_etoposide_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUSC W GEMCITABINE (30)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_gemcitabine <- lusc_clinical[which(lusc_clinical$drug_name == 'gemcitabine' | lusc_clinical$drug_name == 'Gemcitabine' | 
                                                 lusc_clinical$drug_name == 'Gemzar' | lusc_clinical$drug_name == 'GEMZAR' | 
                                                   lusc_clinical$drug_name == 'Gemzar (Gemcitabine'), ]

lusc_clinical_gemcitabine$most_sensitive  <- ifelse(lusc_clinical_gemcitabine$PFS < quantile(lusc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
lusc_clinical_gemcitabine$least_sensitive <- ifelse(lusc_clinical_gemcitabine$PFS > quantile(lusc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_gemcitabine$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_gemcitabine_short <- lusc_clinical_gemcitabine[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_gemcitabine_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_gemcitabine_short$most_sensitive, new_lusc_tcga_gemcitabine_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_gemcitabine_short$most_sensitive, new_lusc_tcga_gemcitabine_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_gemcitabine_short$least_sensitive, new_lusc_tcga_gemcitabine_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_gemcitabine_short$least_sensitive, new_lusc_tcga_gemcitabine_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_gemcitabine_most_sensitive_min, lusc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_gemcitabine_most_sensitive_1se, lusc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_gemcitabine_least_sensitive_min, lusc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_gemcitabine_least_sensitive_1se, lusc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# MESO W CISPLATIN (7)
meso_clinical <- read.csv('Processed_Clinical_Data/meso_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(meso_clinical$most_sensitive)
meso_clinical <- meso_clinical[!na_idx, ]
table(meso_clinical$drug_name)
meso_clinical_cisplatin <- meso_clinical[which(meso_clinical$drug_name == 'Cisplatin' | meso_clinical$drug_name == 'Cisplatinum'), ]

meso_clinical_cisplatin$most_sensitive  <- ifelse(meso_clinical_cisplatin$PFS < quantile(meso_clinical_cisplatin$PFS, probs = .20), 1, 0)
meso_clinical_cisplatin$least_sensitive <- ifelse(meso_clinical_cisplatin$PFS > quantile(meso_clinical_cisplatin$PFS, probs = .80), 1, 0)

meso_gene <- read.csv('Processed_Gene_Expression/meso_tcga_rna_seq_processed.csv', row.names = 1)
colnames(meso_gene) <- gsub('\\.', '-', colnames(meso_gene))
meso_matching_idx <- meso_clinical_cisplatin$submitter_id.samples %in% colnames(meso_gene)
meso_clinical_cisplatin_short <- meso_clinical_cisplatin[meso_matching_idx, ]
meso_matching_idx <- colnames(meso_gene) %in% meso_clinical_cisplatin_short$submitter_id.samples
meso_gene_short <- meso_gene[, meso_matching_idx]
meso_gene_short <- t(meso_gene_short)
meso_gene_short_scaled <- apply(meso_gene_short, 2, scale)

rm(meso_tcga_most_min_auc)
rm(meso_tcga_most_1se_auc)
rm(meso_tcga_least_min_auc)
rm(meso_tcga_least_1se_auc)

new_meso_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = meso_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
meso_tcga_most_min_auc <- auc(meso_clinical_cisplatin_short$most_sensitive, new_meso_tcga_cisplatin_most_sensitive_min)
meso_tcga_most_min_auc <- round(meso_tcga_most_min_auc, digits = 2)

new_meso_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = meso_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
meso_tcga_most_1se_auc <- auc(meso_clinical_cisplatin_short$most_sensitive, new_meso_tcga_cisplatin_most_sensitive_1se)
meso_tcga_most_1se_auc <- round(meso_tcga_most_1se_auc, digits = 2)

new_meso_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = meso_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
meso_tcga_least_min_auc <- auc(meso_clinical_cisplatin_short$least_sensitive, new_meso_tcga_cisplatin_least_sensitive_min)
meso_tcga_least_min_auc <- round(meso_tcga_least_min_auc, digits = 2)

new_meso_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = meso_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
meso_tcga_least_1se_auc <- auc(meso_clinical_cisplatin_short$least_sensitive, new_meso_tcga_cisplatin_least_sensitive_1se)
meso_tcga_least_1se_auc <- round(meso_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_meso_tcga_cisplatin_most_sensitive_min, meso_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_meso_tcga_cisplatin_most_sensitive_1se, meso_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_meso_tcga_cisplatin_least_sensitive_min, meso_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_meso_tcga_cisplatin_least_sensitive_1se, meso_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/meso_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W CISPLATIN (28)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_cisplatin <- ov_clinical[which(ov_clinical$drug_name == 'Ciplastin' | ov_clinical$drug_name == 'Cisplatin' | 
                                                 ov_clinical$drug_name == 'Cisplatin #2-7'), ]

ov_clinical_cisplatin$most_sensitive  <- ifelse(ov_clinical_cisplatin$PFS < quantile(ov_clinical_cisplatin$PFS, probs = .20), 1, 0)
ov_clinical_cisplatin$least_sensitive <- ifelse(ov_clinical_cisplatin$PFS > quantile(ov_clinical_cisplatin$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_cisplatin$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_cisplatin_short <- ov_clinical_cisplatin[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_cisplatin_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_cisplatin_short$most_sensitive, new_ov_tcga_cisplatin_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_cisplatin_short$most_sensitive, new_ov_tcga_cisplatin_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_cisplatin_short$least_sensitive, new_ov_tcga_cisplatin_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_cisplatin_short$least_sensitive, new_ov_tcga_cisplatin_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_cisplatin_most_sensitive_min, ov_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_cisplatin_most_sensitive_1se, ov_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_cisplatin_least_sensitive_min, ov_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_cisplatin_least_sensitive_1se, ov_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[2], lty = 2, main = 'OV treated with cisplatin', lwd = 2)
plot(perf3, col = colors_i_need[2], add = TRUE, lwd = 2)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.4, y = 0.4, legend = c('sensitive model (AUC = 0.80)', 'resistant model (AUC = 0.83)'), lty = c(1,2), lwd = 2, col = c(colors_i_need[2], colors_i_need[2]), bty = 'n', cex = 0.8)
dev.off()

# OV W DOXORUBICIN (5)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_doxorubicin <- ov_clinical[which(ov_clinical$drug_name == 'Doxorubicin' | ov_clinical$drug_name == 'Doxorubicin HCL'), ]

ov_clinical_doxorubicin$most_sensitive  <- ifelse(ov_clinical_doxorubicin$PFS < quantile(ov_clinical_doxorubicin$PFS, probs = .20), 1, 0)
ov_clinical_doxorubicin$least_sensitive <- ifelse(ov_clinical_doxorubicin$PFS > quantile(ov_clinical_doxorubicin$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_doxorubicin$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_doxorubicin_short <- ov_clinical_doxorubicin[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_doxorubicin_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_doxorubicin_most_sensitive_min <- predict(dox_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_doxorubicin_short$most_sensitive, new_ov_tcga_doxorubicin_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_doxorubicin_most_sensitive_1se <- predict(dox_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_doxorubicin_short$most_sensitive, new_ov_tcga_doxorubicin_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_doxorubicin_least_sensitive_min <- predict(dox_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_doxorubicin_short$least_sensitive, new_ov_tcga_doxorubicin_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_doxorubicin_least_sensitive_1se <- predict(dox_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_doxorubicin_short$least_sensitive, new_ov_tcga_doxorubicin_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_doxorubicin_most_sensitive_min, ov_clinical_doxorubicin_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_doxorubicin_most_sensitive_1se, ov_clinical_doxorubicin_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_doxorubicin_least_sensitive_min, ov_clinical_doxorubicin_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_doxorubicin_least_sensitive_1se, ov_clinical_doxorubicin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_doxorubicin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W GEMCITABINE (41)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_gemcitabine <- ov_clinical[which(ov_clinical$drug_name == 'Gamzar' | ov_clinical$drug_name == 'gemcitabine' | 
                                                 ov_clinical$drug_name == 'Gemcitabine' | ov_clinical$drug_name == 'Gemcitabine HCL' | 
                                               ov_clinical$drug_name == 'Gemcitibine' | ov_clinical$drug_name == 'Gemzar'), ]

ov_clinical_gemcitabine$most_sensitive  <- ifelse(ov_clinical_gemcitabine$PFS < quantile(ov_clinical_gemcitabine$PFS, probs = .20), 1, 0)
ov_clinical_gemcitabine$least_sensitive <- ifelse(ov_clinical_gemcitabine$PFS > quantile(ov_clinical_gemcitabine$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_gemcitabine$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_gemcitabine_short <- ov_clinical_gemcitabine[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_gemcitabine_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_gemcitabine_short$most_sensitive, new_ov_tcga_gemcitabine_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_gemcitabine_short$most_sensitive, new_ov_tcga_gemcitabine_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_gemcitabine_short$least_sensitive, new_ov_tcga_gemcitabine_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_gemcitabine_short$least_sensitive, new_ov_tcga_gemcitabine_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_gemcitabine_most_sensitive_min, ov_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_gemcitabine_most_sensitive_1se, ov_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_gemcitabine_least_sensitive_min, ov_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_gemcitabine_least_sensitive_1se, ov_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# PAAD W GEMCITABINE (39)
paad_clinical <- read.csv('Processed_Clinical_Data/paad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(paad_clinical$most_sensitive)
paad_clinical <- paad_clinical[!na_idx, ]
table(paad_clinical$drug_name)
paad_clinical_gemcitabine <- paad_clinical[which(paad_clinical$drug_name == 'Gemcitabine' | paad_clinical$drug_name == 'Gemcitabine Injection' | 
                                                 paad_clinical$drug_name == 'gemzar' | paad_clinical$drug_name == 'Gemzar'), ]

paad_clinical_gemcitabine$most_sensitive  <- ifelse(paad_clinical_gemcitabine$PFS < quantile(paad_clinical_gemcitabine$PFS, probs = .20), 1, 0)
paad_clinical_gemcitabine$least_sensitive <- ifelse(paad_clinical_gemcitabine$PFS > quantile(paad_clinical_gemcitabine$PFS, probs = .80), 1, 0)

paad_gene <- read.csv('Processed_Gene_Expression/paad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(paad_gene) <- gsub('\\.', '-', colnames(paad_gene))
paad_matching_idx <- paad_clinical_gemcitabine$submitter_id.samples %in% colnames(paad_gene)
paad_clinical_gemcitabine_short <- paad_clinical_gemcitabine[paad_matching_idx, ]
paad_matching_idx <- colnames(paad_gene) %in% paad_clinical_gemcitabine_short$submitter_id.samples
paad_gene_short <- paad_gene[, paad_matching_idx]
paad_gene_short <- t(paad_gene_short)
paad_gene_short_scaled <- apply(paad_gene_short, 2, scale)

rm(paad_tcga_most_min_auc)
rm(paad_tcga_most_1se_auc)
rm(paad_tcga_least_min_auc)
rm(paad_tcga_least_1se_auc)

new_paad_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_most_min_auc <- auc(paad_clinical_gemcitabine_short$most_sensitive, new_paad_tcga_gemcitabine_most_sensitive_min)
paad_tcga_most_min_auc <- round(paad_tcga_most_min_auc, digits = 2)

new_paad_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_most_1se_auc <- auc(paad_clinical_gemcitabine_short$most_sensitive, new_paad_tcga_gemcitabine_most_sensitive_1se)
paad_tcga_most_1se_auc <- round(paad_tcga_most_1se_auc, digits = 2)

new_paad_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_least_min_auc <- auc(paad_clinical_gemcitabine_short$least_sensitive, new_paad_tcga_gemcitabine_least_sensitive_min)
paad_tcga_least_min_auc <- round(paad_tcga_least_min_auc, digits = 2)

new_paad_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_least_1se_auc <- auc(paad_clinical_gemcitabine_short$least_sensitive, new_paad_tcga_gemcitabine_least_sensitive_1se)
paad_tcga_least_1se_auc <- round(paad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_paad_tcga_gemcitabine_most_sensitive_min, paad_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_paad_tcga_gemcitabine_most_sensitive_1se, paad_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_paad_tcga_gemcitabine_least_sensitive_min, paad_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_paad_tcga_gemcitabine_least_sensitive_1se, paad_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/paad_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# SARC W GEMCITABINE (17)
sarc_clinical <- read.csv('Processed_Clinical_Data/sarc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(sarc_clinical$most_sensitive)
sarc_clinical <- sarc_clinical[!na_idx, ]
table(sarc_clinical$drug_name)
sarc_clinical_gemcitabine <- sarc_clinical[which(sarc_clinical$drug_name == 'gemcitabine' | sarc_clinical$drug_name == 'Gemcitabine' | 
                                                 sarc_clinical$drug_name == 'GEMCITABINE' | sarc_clinical$drug_name == 'Gemzar'), ]

sarc_clinical_gemcitabine$most_sensitive  <- ifelse(sarc_clinical_gemcitabine$PFS < quantile(sarc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
sarc_clinical_gemcitabine$least_sensitive <- ifelse(sarc_clinical_gemcitabine$PFS > quantile(sarc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

sarc_gene <- read.csv('Processed_Gene_Expression/sarc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(sarc_gene) <- gsub('\\.', '-', colnames(sarc_gene))
sarc_matching_idx <- sarc_clinical_gemcitabine$submitter_id.samples %in% colnames(sarc_gene)
sarc_clinical_gemcitabine_short <- sarc_clinical_gemcitabine[sarc_matching_idx, ]
sarc_matching_idx <- colnames(sarc_gene) %in% sarc_clinical_gemcitabine_short$submitter_id.samples
sarc_gene_short <- sarc_gene[, sarc_matching_idx]
sarc_gene_short <- t(sarc_gene_short)
sarc_gene_short_scaled <- apply(sarc_gene_short, 2, scale)

rm(sarc_tcga_most_min_auc)
rm(sarc_tcga_most_1se_auc)
rm(sarc_tcga_least_min_auc)
rm(sarc_tcga_least_1se_auc)

new_sarc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_most_min_auc <- auc(sarc_clinical_gemcitabine_short$most_sensitive, new_sarc_tcga_gemcitabine_most_sensitive_min)
sarc_tcga_most_min_auc <- round(sarc_tcga_most_min_auc, digits = 2)

new_sarc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_most_1se_auc <- auc(sarc_clinical_gemcitabine_short$most_sensitive, new_sarc_tcga_gemcitabine_most_sensitive_1se)
sarc_tcga_most_1se_auc <- round(sarc_tcga_most_1se_auc, digits = 2)

new_sarc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_least_min_auc <- auc(sarc_clinical_gemcitabine_short$least_sensitive, new_sarc_tcga_gemcitabine_least_sensitive_min)
sarc_tcga_least_min_auc <- round(sarc_tcga_least_min_auc, digits = 2)

new_sarc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_least_1se_auc <- auc(sarc_clinical_gemcitabine_short$least_sensitive, new_sarc_tcga_gemcitabine_least_sensitive_1se)
sarc_tcga_least_1se_auc <- round(sarc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_sarc_tcga_gemcitabine_most_sensitive_min, sarc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_sarc_tcga_gemcitabine_most_sensitive_1se, sarc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_sarc_tcga_gemcitabine_least_sensitive_min, sarc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_sarc_tcga_gemcitabine_least_sensitive_1se, sarc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/sarc_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# SKCM W TEMOZOLOMIDE (4)
skcm_clinical <- read.csv('Processed_Clinical_Data/skcm_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(skcm_clinical$most_sensitive)
skcm_clinical <- skcm_clinical[!na_idx, ]
table(skcm_clinical$drug_name)
skcm_clinical_temozolomide <- skcm_clinical[which(skcm_clinical$drug_name == 'Temodal' | skcm_clinical$drug_name == 'Temodar' | 
                                                 skcm_clinical$drug_name == 'Temozolomide'), ]

skcm_clinical_temozolomide$most_sensitive  <- ifelse(skcm_clinical_temozolomide$PFS < quantile(skcm_clinical_temozolomide$PFS, probs = .20), 1, 0)
skcm_clinical_temozolomide$least_sensitive <- ifelse(skcm_clinical_temozolomide$PFS > quantile(skcm_clinical_temozolomide$PFS, probs = .80), 1, 0)

skcm_gene <- read.csv('Processed_Gene_Expression/skcm_tcga_rna_seq_processed.csv', row.names = 1)
colnames(skcm_gene) <- gsub('\\.', '-', colnames(skcm_gene))
skcm_matching_idx <- skcm_clinical_temozolomide$submitter_id.samples %in% colnames(skcm_gene)
skcm_clinical_temozolomide_short <- skcm_clinical_temozolomide[skcm_matching_idx, ]
skcm_matching_idx <- colnames(skcm_gene) %in% skcm_clinical_temozolomide_short$submitter_id.samples
skcm_gene_short <- skcm_gene[, skcm_matching_idx]
skcm_gene_short <- t(skcm_gene_short)
skcm_gene_short_scaled <- apply(skcm_gene_short, 2, scale)

rm(skcm_tcga_most_min_auc)
rm(skcm_tcga_most_1se_auc)
rm(skcm_tcga_least_min_auc)
rm(skcm_tcga_least_1se_auc)

new_skcm_tcga_temozolomide_most_sensitive_min <- predict(temozolomide_most_fit_elnet, newx = skcm_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
skcm_tcga_most_min_auc <- auc(skcm_clinical_temozolomide_short$most_sensitive, new_skcm_tcga_temozolomide_most_sensitive_min)
skcm_tcga_most_min_auc <- round(skcm_tcga_most_min_auc, digits = 2)

new_skcm_tcga_temozolomide_most_sensitive_1se <- predict(temozolomide_most_fit_elnet, newx = skcm_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
skcm_tcga_most_1se_auc <- auc(skcm_clinical_temozolomide_short$most_sensitive, new_skcm_tcga_temozolomide_most_sensitive_1se)
skcm_tcga_most_1se_auc <- round(skcm_tcga_most_1se_auc, digits = 2)

new_skcm_tcga_temozolomide_least_sensitive_min <- predict(temozolomide_least_fit_elnet, newx = skcm_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
skcm_tcga_least_min_auc <- auc(skcm_clinical_temozolomide_short$least_sensitive, new_skcm_tcga_temozolomide_least_sensitive_min)
skcm_tcga_least_min_auc <- round(skcm_tcga_least_min_auc, digits = 2)

new_skcm_tcga_temozolomide_least_sensitive_1se <- predict(temozolomide_least_fit_elnet, newx = skcm_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
skcm_tcga_least_1se_auc <- auc(skcm_clinical_temozolomide_short$least_sensitive, new_skcm_tcga_temozolomide_least_sensitive_1se)
skcm_tcga_least_1se_auc <- round(skcm_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_skcm_tcga_temozolomide_most_sensitive_min, skcm_clinical_temozolomide_short$most_sensitive)
pred2 <- prediction(new_skcm_tcga_temozolomide_most_sensitive_1se, skcm_clinical_temozolomide_short$most_sensitive)
pred3 <- prediction(new_skcm_tcga_temozolomide_least_sensitive_min, skcm_clinical_temozolomide_short$least_sensitive)
pred4 <- prediction(new_skcm_tcga_temozolomide_least_sensitive_1se, skcm_clinical_temozolomide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/skcm_temozolomide_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# STAD W CISPLATIN (13)
stad_clinical <- read.csv('Processed_Clinical_Data/stad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(stad_clinical$most_sensitive)
stad_clinical <- stad_clinical[!na_idx, ]
table(stad_clinical$drug_name)
stad_clinical_cisplatin <- stad_clinical[which(stad_clinical$drug_name == 'Cisplatin'), ]

stad_clinical_cisplatin$most_sensitive  <- ifelse(stad_clinical_cisplatin$PFS < quantile(stad_clinical_cisplatin$PFS, probs = .20), 1, 0)
stad_clinical_cisplatin$least_sensitive <- ifelse(stad_clinical_cisplatin$PFS > quantile(stad_clinical_cisplatin$PFS, probs = .80), 1, 0)

stad_gene <- read.csv('Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(stad_gene) <- gsub('\\.', '-', colnames(stad_gene))
stad_matching_idx <- stad_clinical_cisplatin$submitter_id.samples %in% colnames(stad_gene)
stad_clinical_cisplatin_short <- stad_clinical_cisplatin[stad_matching_idx, ]
stad_matching_idx <- colnames(stad_gene) %in% stad_clinical_cisplatin_short$submitter_id.samples
stad_gene_short <- stad_gene[, stad_matching_idx]
stad_gene_short <- t(stad_gene_short)
stad_gene_short_scaled <- apply(stad_gene_short, 2, scale)

rm(stad_tcga_most_min_auc)
rm(stad_tcga_most_1se_auc)
rm(stad_tcga_least_min_auc)
rm(stad_tcga_least_1se_auc)

new_stad_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_most_min_auc <- auc(stad_clinical_cisplatin_short$most_sensitive, new_stad_tcga_cisplatin_most_sensitive_min)
stad_tcga_most_min_auc <- round(stad_tcga_most_min_auc, digits = 2)

new_stad_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_most_1se_auc <- auc(stad_clinical_cisplatin_short$most_sensitive, new_stad_tcga_cisplatin_most_sensitive_1se)
stad_tcga_most_1se_auc <- round(stad_tcga_most_1se_auc, digits = 2)

new_stad_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_least_min_auc <- auc(stad_clinical_cisplatin_short$least_sensitive, new_stad_tcga_cisplatin_least_sensitive_min)
stad_tcga_least_min_auc <- round(stad_tcga_least_min_auc, digits = 2)

new_stad_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_least_1se_auc <- auc(stad_clinical_cisplatin_short$least_sensitive, new_stad_tcga_cisplatin_least_sensitive_1se)
stad_tcga_least_1se_auc <- round(stad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_stad_tcga_cisplatin_most_sensitive_min, stad_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_stad_tcga_cisplatin_most_sensitive_1se, stad_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_stad_tcga_cisplatin_least_sensitive_min, stad_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_stad_tcga_cisplatin_least_sensitive_1se, stad_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/stad_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# STAD W ETOPOSIDE (7)
stad_clinical <- read.csv('Processed_Clinical_Data/stad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(stad_clinical$most_sensitive)
stad_clinical <- stad_clinical[!na_idx, ]
table(stad_clinical$drug_name)
stad_clinical_etoposide <- stad_clinical[which(stad_clinical$drug_name == 'Etoposide'), ]

stad_clinical_etoposide$most_sensitive  <- ifelse(stad_clinical_etoposide$PFS < quantile(stad_clinical_etoposide$PFS, probs = .20), 1, 0)
stad_clinical_etoposide$least_sensitive <- ifelse(stad_clinical_etoposide$PFS > quantile(stad_clinical_etoposide$PFS, probs = .80), 1, 0)

stad_gene <- read.csv('Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(stad_gene) <- gsub('\\.', '-', colnames(stad_gene))
stad_matching_idx <- stad_clinical_etoposide$submitter_id.samples %in% colnames(stad_gene)
stad_clinical_etoposide_short <- stad_clinical_etoposide[stad_matching_idx, ]
stad_matching_idx <- colnames(stad_gene) %in% stad_clinical_etoposide_short$submitter_id.samples
stad_gene_short <- stad_gene[, stad_matching_idx]
stad_gene_short <- t(stad_gene_short)
stad_gene_short_scaled <- apply(stad_gene_short, 2, scale)
for(i in 1:nrow(stad_gene_short_scaled)) {
  for (j in 1:ncol(stad_gene_short_scaled)) {
    if (is.na(stad_gene_short_scaled[i,j]) == TRUE) {
      stad_gene_short_scaled[i,j] <- 0
    }
  }
}

rm(stad_tcga_most_min_auc)
rm(stad_tcga_most_1se_auc)
rm(stad_tcga_least_min_auc)
rm(stad_tcga_least_1se_auc)

new_stad_tcga_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_most_min_auc <- auc(stad_clinical_etoposide_short$most_sensitive, new_stad_tcga_etoposide_most_sensitive_min)
stad_tcga_most_min_auc <- round(stad_tcga_most_min_auc, digits = 2)

new_stad_tcga_etoposide_most_sensitive_1se <- predict(etoposide_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_most_1se_auc <- auc(stad_clinical_etoposide_short$most_sensitive, new_stad_tcga_etoposide_most_sensitive_1se)
stad_tcga_most_1se_auc <- round(stad_tcga_most_1se_auc, digits = 2)

new_stad_tcga_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_least_min_auc <- auc(stad_clinical_etoposide_short$least_sensitive, new_stad_tcga_etoposide_least_sensitive_min)
stad_tcga_least_min_auc <- round(stad_tcga_least_min_auc, digits = 2)

new_stad_tcga_etoposide_least_sensitive_1se <- predict(etoposide_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_least_1se_auc <- auc(stad_clinical_etoposide_short$least_sensitive, new_stad_tcga_etoposide_least_sensitive_1se)
stad_tcga_least_1se_auc <- round(stad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_stad_tcga_etoposide_most_sensitive_min, stad_clinical_etoposide_short$most_sensitive)
pred2 <- prediction(new_stad_tcga_etoposide_most_sensitive_1se, stad_clinical_etoposide_short$most_sensitive)
pred3 <- prediction(new_stad_tcga_etoposide_least_sensitive_min, stad_clinical_etoposide_short$least_sensitive)
pred4 <- prediction(new_stad_tcga_etoposide_least_sensitive_1se, stad_clinical_etoposide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/stad_etoposide_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
#plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# TGCT W CISPLATIN (48)
tgct_clinical <- read.csv('Processed_Clinical_Data/tgct_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(tgct_clinical$most_sensitive)
tgct_clinical <- tgct_clinical[!na_idx, ]
table(tgct_clinical$drug_name)
tgct_clinical_cisplatin <- tgct_clinical[which(tgct_clinical$drug_name == 'cisplatin' | tgct_clinical$drug_name == 'Cisplatin'), ]

tgct_clinical_cisplatin$most_sensitive  <- ifelse(tgct_clinical_cisplatin$PFS < quantile(tgct_clinical_cisplatin$PFS, probs = .20), 1, 0)
tgct_clinical_cisplatin$least_sensitive <- ifelse(tgct_clinical_cisplatin$PFS > quantile(tgct_clinical_cisplatin$PFS, probs = .80), 1, 0)

tgct_gene <- read.csv('Processed_Gene_Expression/tgct_tcga_rna_seq_processed.csv', row.names = 1)
colnames(tgct_gene) <- gsub('\\.', '-', colnames(tgct_gene))
tgct_matching_idx <- tgct_clinical_cisplatin$submitter_id.samples %in% colnames(tgct_gene)
tgct_clinical_cisplatin_short <- tgct_clinical_cisplatin[tgct_matching_idx, ]
tgct_matching_idx <- colnames(tgct_gene) %in% tgct_clinical_cisplatin_short$submitter_id.samples
tgct_gene_short <- tgct_gene[, tgct_matching_idx]
tgct_gene_short <- t(tgct_gene_short)
tgct_gene_short_scaled <- apply(tgct_gene_short, 2, scale)

rm(tgct_tcga_most_min_auc)
rm(tgct_tcga_most_1se_auc)
rm(tgct_tcga_least_min_auc)
rm(tgct_tcga_least_1se_auc)

new_tgct_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = tgct_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
tgct_tcga_most_min_auc <- auc(tgct_clinical_cisplatin_short$most_sensitive, new_tgct_tcga_cisplatin_most_sensitive_min)
tgct_tcga_most_min_auc <- round(tgct_tcga_most_min_auc, digits = 2)

new_tgct_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = tgct_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
tgct_tcga_most_1se_auc <- auc(tgct_clinical_cisplatin_short$most_sensitive, new_tgct_tcga_cisplatin_most_sensitive_1se)
tgct_tcga_most_1se_auc <- round(tgct_tcga_most_1se_auc, digits = 2)

new_tgct_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = tgct_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
tgct_tcga_least_min_auc <- auc(tgct_clinical_cisplatin_short$least_sensitive, new_tgct_tcga_cisplatin_least_sensitive_min)
tgct_tcga_least_min_auc <- round(tgct_tcga_least_min_auc, digits = 2)

new_tgct_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = tgct_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
tgct_tcga_least_1se_auc <- auc(tgct_clinical_cisplatin_short$least_sensitive, new_tgct_tcga_cisplatin_least_sensitive_1se)
tgct_tcga_least_1se_auc <- round(tgct_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_tgct_tcga_cisplatin_most_sensitive_min, tgct_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_tgct_tcga_cisplatin_most_sensitive_1se, tgct_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_tgct_tcga_cisplatin_least_sensitive_min, tgct_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_tgct_tcga_cisplatin_least_sensitive_1se, tgct_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/tgct_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCEC W CISPLATIN (5)
ucec_clinical <- read.csv('Processed_Clinical_Data/ucec_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucec_clinical$most_sensitive)
ucec_clinical <- ucec_clinical[!na_idx, ]
table(ucec_clinical$drug_name)
ucec_clinical_cisplatin <- ucec_clinical[which(ucec_clinical$drug_name == 'cisplatin' | ucec_clinical$drug_name == 'Cisplatin'), ]

ucec_clinical_cisplatin$most_sensitive  <- ifelse(ucec_clinical_cisplatin$PFS < quantile(ucec_clinical_cisplatin$PFS, probs = .20), 1, 0)
ucec_clinical_cisplatin$least_sensitive <- ifelse(ucec_clinical_cisplatin$PFS > quantile(ucec_clinical_cisplatin$PFS, probs = .80), 1, 0)

ucec_gene <- read.csv('Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucec_gene) <- gsub('\\.', '-', colnames(ucec_gene))
ucec_matching_idx <- ucec_clinical_cisplatin$submitter_id.samples %in% colnames(ucec_gene)
ucec_clinical_cisplatin_short <- ucec_clinical_cisplatin[ucec_matching_idx, ]
ucec_matching_idx <- colnames(ucec_gene) %in% ucec_clinical_cisplatin_short$submitter_id.samples
ucec_gene_short <- ucec_gene[, ucec_matching_idx]
ucec_gene_short <- t(ucec_gene_short)
ucec_gene_short_scaled <- apply(ucec_gene_short, 2, scale)

rm(ucec_tcga_most_min_auc)
rm(ucec_tcga_most_1se_auc)
rm(ucec_tcga_least_min_auc)
rm(ucec_tcga_least_1se_auc)

new_ucec_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_most_min_auc <- auc(ucec_clinical_cisplatin_short$most_sensitive, new_ucec_tcga_cisplatin_most_sensitive_min)
ucec_tcga_most_min_auc <- round(ucec_tcga_most_min_auc, digits = 2)

new_ucec_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_most_1se_auc <- auc(ucec_clinical_cisplatin_short$most_sensitive, new_ucec_tcga_cisplatin_most_sensitive_1se)
ucec_tcga_most_1se_auc <- round(ucec_tcga_most_1se_auc, digits = 2)

new_ucec_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_least_min_auc <- auc(ucec_clinical_cisplatin_short$least_sensitive, new_ucec_tcga_cisplatin_least_sensitive_min)
ucec_tcga_least_min_auc <- round(ucec_tcga_least_min_auc, digits = 2)

new_ucec_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_least_1se_auc <- auc(ucec_clinical_cisplatin_short$least_sensitive, new_ucec_tcga_cisplatin_least_sensitive_1se)
ucec_tcga_least_1se_auc <- round(ucec_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucec_tcga_cisplatin_most_sensitive_min, ucec_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_ucec_tcga_cisplatin_most_sensitive_1se, ucec_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_ucec_tcga_cisplatin_least_sensitive_min, ucec_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_ucec_tcga_cisplatin_least_sensitive_1se, ucec_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucec_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCEC W DOXORUBICIN (7)
ucec_clinical <- read.csv('Processed_Clinical_Data/ucec_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucec_clinical$most_sensitive)
ucec_clinical <- ucec_clinical[!na_idx, ]
table(ucec_clinical$drug_name)
ucec_clinical_doxorubicin <- ucec_clinical[which(ucec_clinical$drug_name == 'doxorubicin' | ucec_clinical$drug_name == 'Doxorubicin'), ]

ucec_clinical_doxorubicin$most_sensitive  <- ifelse(ucec_clinical_doxorubicin$PFS < quantile(ucec_clinical_doxorubicin$PFS, probs = .20), 1, 0)
ucec_clinical_doxorubicin$least_sensitive <- ifelse(ucec_clinical_doxorubicin$PFS > quantile(ucec_clinical_doxorubicin$PFS, probs = .80), 1, 0)

ucec_gene <- read.csv('Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucec_gene) <- gsub('\\.', '-', colnames(ucec_gene))
ucec_matching_idx <- ucec_clinical_doxorubicin$submitter_id.samples %in% colnames(ucec_gene)
ucec_clinical_doxorubicin_short <- ucec_clinical_doxorubicin[ucec_matching_idx, ]
ucec_matching_idx <- colnames(ucec_gene) %in% ucec_clinical_doxorubicin_short$submitter_id.samples
ucec_gene_short <- ucec_gene[, ucec_matching_idx]
ucec_gene_short <- t(ucec_gene_short)
ucec_gene_short_scaled <- apply(ucec_gene_short, 2, scale)

rm(ucec_tcga_most_min_auc)
rm(ucec_tcga_most_1se_auc)
rm(ucec_tcga_least_min_auc)
rm(ucec_tcga_least_1se_auc)

new_ucec_tcga_doxorubicin_most_sensitive_min <- predict(dox_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_most_min_auc <- auc(ucec_clinical_doxorubicin_short$most_sensitive, new_ucec_tcga_doxorubicin_most_sensitive_min)
ucec_tcga_most_min_auc <- round(ucec_tcga_most_min_auc, digits = 2)

new_ucec_tcga_doxorubicin_most_sensitive_1se <- predict(dox_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_most_1se_auc <- auc(ucec_clinical_doxorubicin_short$most_sensitive, new_ucec_tcga_doxorubicin_most_sensitive_1se)
ucec_tcga_most_1se_auc <- round(ucec_tcga_most_1se_auc, digits = 2)

new_ucec_tcga_doxorubicin_least_sensitive_min <- predict(dox_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_least_min_auc <- auc(ucec_clinical_doxorubicin_short$least_sensitive, new_ucec_tcga_doxorubicin_least_sensitive_min)
ucec_tcga_least_min_auc <- round(ucec_tcga_least_min_auc, digits = 2)

new_ucec_tcga_doxorubicin_least_sensitive_1se <- predict(dox_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_least_1se_auc <- auc(ucec_clinical_doxorubicin_short$least_sensitive, new_ucec_tcga_doxorubicin_least_sensitive_1se)
ucec_tcga_least_1se_auc <- round(ucec_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucec_tcga_doxorubicin_most_sensitive_min, ucec_clinical_doxorubicin_short$most_sensitive)
pred2 <- prediction(new_ucec_tcga_doxorubicin_most_sensitive_1se, ucec_clinical_doxorubicin_short$most_sensitive)
pred3 <- prediction(new_ucec_tcga_doxorubicin_least_sensitive_min, ucec_clinical_doxorubicin_short$least_sensitive)
pred4 <- prediction(new_ucec_tcga_doxorubicin_least_sensitive_1se, ucec_clinical_doxorubicin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucec_doxorubicin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCS W CISPLATIN (4)
ucs_clinical <- read.csv('Processed_Clinical_Data/ucs_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucs_clinical$most_sensitive)
ucs_clinical <- ucs_clinical[!na_idx, ]
table(ucs_clinical$drug_name)
ucs_clinical_cisplatin <- ucs_clinical[which(ucs_clinical$drug_name == 'Cisplatin'), ]

ucs_clinical_cisplatin$most_sensitive  <- ifelse(ucs_clinical_cisplatin$PFS < quantile(ucs_clinical_cisplatin$PFS, probs = .20), 1, 0)
ucs_clinical_cisplatin$least_sensitive <- ifelse(ucs_clinical_cisplatin$PFS > quantile(ucs_clinical_cisplatin$PFS, probs = .80), 1, 0)

ucs_gene <- read.csv('Processed_Gene_Expression/ucs_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucs_gene) <- gsub('\\.', '-', colnames(ucs_gene))
ucs_matching_idx <- ucs_clinical_cisplatin$submitter_id.samples %in% colnames(ucs_gene)
ucs_clinical_cisplatin_short <- ucs_clinical_cisplatin[ucs_matching_idx, ]
ucs_matching_idx <- colnames(ucs_gene) %in% ucs_clinical_cisplatin_short$submitter_id.samples
ucs_gene_short <- ucs_gene[, ucs_matching_idx]
ucs_gene_short <- t(ucs_gene_short)
ucs_gene_short_scaled <- apply(ucs_gene_short, 2, scale)

rm(ucs_tcga_most_min_auc)
rm(ucs_tcga_most_1se_auc)
rm(ucs_tcga_least_min_auc)
rm(ucs_tcga_least_1se_auc)

new_ucs_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_most_min_auc <- auc(ucs_clinical_cisplatin_short$most_sensitive, new_ucs_tcga_cisplatin_most_sensitive_min)
ucs_tcga_most_min_auc <- round(ucs_tcga_most_min_auc, digits = 2)

new_ucs_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_most_1se_auc <- auc(ucs_clinical_cisplatin_short$most_sensitive, new_ucs_tcga_cisplatin_most_sensitive_1se)
ucs_tcga_most_1se_auc <- round(ucs_tcga_most_1se_auc, digits = 2)

new_ucs_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_least_min_auc <- auc(ucs_clinical_cisplatin_short$least_sensitive, new_ucs_tcga_cisplatin_least_sensitive_min)
ucs_tcga_least_min_auc <- round(ucs_tcga_least_min_auc, digits = 2)

new_ucs_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_least_1se_auc <- auc(ucs_clinical_cisplatin_short$least_sensitive, new_ucs_tcga_cisplatin_least_sensitive_1se)
ucs_tcga_least_1se_auc <- round(ucs_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucs_tcga_cisplatin_most_sensitive_min, ucs_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_ucs_tcga_cisplatin_most_sensitive_1se, ucs_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_ucs_tcga_cisplatin_least_sensitive_min, ucs_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_ucs_tcga_cisplatin_least_sensitive_1se, ucs_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucs_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

### CCLE ----
## load clinical data ----
carboplatin_ccle <- read.csv('Processed_Clinical_Data/carboplatin_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
cyclophosphamide_ccle <- read.csv('Processed_Clinical_Data/cyclophosphamide_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
docetaxel_ccle <- read.csv('Processed_Clinical_Data/docetaxel_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
fluorouracil_ccle <- read.csv('Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
gemcitabine_ccle <- read.csv('Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)

## load gene expression data ----
ccle_microarray <- read.csv('Processed_Gene_Expression/ccle_microarray_processed.csv', row.names = 1)
colnames(ccle_microarray) <- gsub('_.*$', '', colnames(ccle_microarray))
colnames(ccle_microarray) <- gsub('^X', '', colnames(ccle_microarray))


## set up data for model building ----
carboplatin_lines <- carboplatin_ccle$Cell.line.name #650
carboplatin_lines <- as.character(carboplatin_lines)
cyclophosphamide_lines <- cyclophosphamide_ccle$Cell.line.name #335
cyclophosphamide_lines <- as.character(cyclophosphamide_lines)
docetaxel_lines <- docetaxel_ccle$Cell.line.name #302
docetaxel_lines <- as.character(docetaxel_lines)
fluorouracil_lines <- fluorouracil_ccle$Cell.line.name #652
fluorouracil_lines <- as.character(fluorouracil_lines)
gemcitabine_lines <- gemcitabine_ccle$Cell.line.name #589
gemcitabine_lines <- as.character(gemcitabine_lines)


ccle <- data.frame(t(ccle_microarray))
#855 x 14209

set.seed(5)

random_sample <- sample(x = rownames(ccle), size = nrow(ccle)/2)


ccle_train         <- ccle[random_sample, ] #427 x 14209

ccle_test          <- ccle[which(rownames(ccle) %ni% random_sample), ] #481 x 14209

intersect(rownames(ccle_train), rownames(ccle_test))

carboplatin_rna_seq_train     <- ccle_train[intersect(carboplatin_lines, rownames(ccle_train)), ]
# 312 x 14209
carboplatin_rna_seq_test      <- ccle_test[intersect(carboplatin_lines, rownames(ccle_test)), ]
# 306 x 14209

cyclophosphamide_rna_seq_train     <- ccle_train[intersect(cyclophosphamide_lines, rownames(ccle_train)), ]
# 165 x 14209
cyclophosphamide_rna_seq_test      <- ccle_test[intersect(cyclophosphamide_lines, rownames(ccle_test)), ]
# 151 x 14209

docetaxel_rna_seq_train     <- ccle_train[intersect(docetaxel_lines, rownames(ccle_train)), ]
# 142 x 14209
docetaxel_rna_seq_test      <- ccle_test[intersect(docetaxel_lines, rownames(ccle_test)), ]
# 155 x 14209

fluorouracil_rna_seq_train     <- ccle_train[intersect(fluorouracil_lines, rownames(ccle_train)), ]
# 313 x 14209
fluorouracil_rna_seq_test      <- ccle_test[intersect(fluorouracil_lines, rownames(ccle_test)), ]
# 307 x 14209

gemcitabine_rna_seq_train     <- ccle_train[intersect(gemcitabine_lines, rownames(ccle_train)), ]
# 276 x 14209
gemcitabine_rna_seq_test      <- ccle_test[intersect(gemcitabine_lines, rownames(ccle_test)), ]
# 286 x 14209


carboplatin_train        <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_train)), ]
carboplatin_test         <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_test)), ]

cyclophosphamide_train        <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_train)), ]
cyclophosphamide_test         <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_test)), ]

docetaxel_train        <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_train)), ]
docetaxel_test         <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_test)), ]

fluorouracil_train        <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_train)), ]
fluorouracil_test         <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_test)), ]

gemcitabine_train        <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test         <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_test)), ]



carboplatin_rna_seq_train_scaled          <- apply(carboplatin_rna_seq_train, 2, scale)
carboplatin_rna_seq_test_scaled           <- as.data.frame(apply(carboplatin_rna_seq_test, 2, scale))

cyclophosphamide_rna_seq_train_scaled          <- apply(cyclophosphamide_rna_seq_train, 2, scale)
cyclophosphamide_rna_seq_test_scaled           <- as.data.frame(apply(cyclophosphamide_rna_seq_test, 2, scale))

docetaxel_rna_seq_train_scaled          <- apply(docetaxel_rna_seq_train, 2, scale)
docetaxel_rna_seq_test_scaled           <- as.data.frame(apply(docetaxel_rna_seq_test, 2, scale))

fluorouracil_rna_seq_train_scaled          <- apply(fluorouracil_rna_seq_train, 2, scale)
fluorouracil_rna_seq_test_scaled           <- as.data.frame(apply(fluorouracil_rna_seq_test, 2, scale))

gemcitabine_rna_seq_train_scaled          <- apply(gemcitabine_rna_seq_train, 2, scale)
gemcitabine_rna_seq_test_scaled           <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))



### load models ----
carboplatin_ccle_most_fit_elnet      <- readRDS('GLM_Models/carboplatin_ccle_most_model.rds')
carboplatin_ccle_least_fit_elnet     <- readRDS('GLM_Models/carboplatin_ccle_least_model.rds')

cyclophosphamide_ccle_most_fit_elnet      <- readRDS('GLM_Models/cyclophosphamide_ccle_most_model.rds')
cyclophosphamide_ccle_least_fit_elnet     <- readRDS('GLM_Models/cyclophosphamide_ccle_least_model.rds')

docetaxel_ccle_most_fit_elnet   <- readRDS('GLM_Models/docetaxel_ccle_most_model.rds')
docetaxel_ccle_least_fit_elnet  <- readRDS('GLM_Models/docetaxel_ccle_least_model.rds')

fluorouracil_ccle_most_fit_elnet <- readRDS('GLM_Models/fluorouracil_ccle_most_model.rds')
fluorouracil_ccle_least_fit_elnet <- readRDS('GLM_Models/fluorouracil_ccle_least_model.rds')

gemcitabine_ccle_most_fit_elnet <- readRDS('GLM_Models/gemcitabine_ccle_most_model.rds')
gemcitabine_ccle_least_fit_elnet <- readRDS('GLM_Models/gemcitabine_ccle_least_model.rds')


### get accuracy on testing data --------

## CARBOPLATIN
new_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

carboplatin_most_test_gdsc_auc_1se <- auc(carboplatin_test$most_sensitive, new_carboplatin_most_sensitive_1se)
carboplatin_most_test_gdsc_auc_1se <- round(carboplatin_most_test_gdsc_auc_1se, digits = 2)


new_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

carboplatin_least_test_gdsc_auc_1se <- auc(carboplatin_test$least_sensitive, new_carboplatin_least_sensitive_1se)
carboplatin_least_test_gdsc_auc_1se <- round(carboplatin_least_test_gdsc_auc_1se, digits = 2)

## cyclophosphamide
new_cyclophosphamide_most_sensitive_1se <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cyclophosphamide_most_test_gdsc_auc_1se <- auc(cyclophosphamide_test$most_sensitive, new_cyclophosphamide_most_sensitive_1se)
cyclophosphamide_most_test_gdsc_auc_1se <- round(cyclophosphamide_most_test_gdsc_auc_1se, digits = 2)


new_cyclophosphamide_least_sensitive_1se <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cyclophosphamide_least_test_gdsc_auc_1se <- auc(cyclophosphamide_test$least_sensitive, new_cyclophosphamide_least_sensitive_1se)
cyclophosphamide_least_test_gdsc_auc_1se <- round(cyclophosphamide_least_test_gdsc_auc_1se, digits = 2)

## docetaxel
new_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

docetaxel_most_test_gdsc_auc_1se <- auc(docetaxel_test$most_sensitive, new_docetaxel_most_sensitive_1se)
docetaxel_most_test_gdsc_auc_1se <- round(docetaxel_most_test_gdsc_auc_1se, digits = 2)


new_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

docetaxel_least_test_gdsc_auc_min <- auc(docetaxel_test$least_sensitive, new_docetaxel_least_sensitive_min)
docetaxel_least_test_gdsc_auc_min <- round(docetaxel_least_test_gdsc_auc_min, digits = 2)

## fluorouracil
new_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

fluorouracil_most_test_gdsc_auc_1se <- auc(fluorouracil_test$most_sensitive, new_fluorouracil_most_sensitive_1se)
fluorouracil_most_test_gdsc_auc_1se <- round(fluorouracil_most_test_gdsc_auc_1se, digits = 2)


new_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

fluorouracil_least_test_gdsc_auc_1se <- auc(fluorouracil_test$least_sensitive, new_fluorouracil_least_sensitive_1se)
fluorouracil_least_test_gdsc_auc_1se <- round(fluorouracil_least_test_gdsc_auc_1se, digits = 2)

## gemcitabine
new_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_most_test_gdsc_auc_1se <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_1se)
gemcitabine_most_test_gdsc_auc_1se <- round(gemcitabine_most_test_gdsc_auc_1se, digits = 2)


new_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_least_test_gdsc_auc_1se <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_1se)
gemcitabine_least_test_gdsc_auc_1se <- round(gemcitabine_least_test_gdsc_auc_1se, digits = 2)


overall_auc <- c(carboplatin_most_test_gdsc_auc_1se, carboplatin_least_test_gdsc_auc_1se, 
                 cyclophosphamide_most_test_gdsc_auc_1se, cyclophosphamide_least_test_gdsc_auc_1se, 
                 docetaxel_most_test_gdsc_auc_1se, docetaxel_least_test_gdsc_auc_min, 
                 fluorouracil_most_test_gdsc_auc_1se, fluorouracil_least_test_gdsc_auc_1se, 
                 gemcitabine_most_test_gdsc_auc_1se, gemcitabine_least_test_gdsc_auc_1se)

# subsetting lines by cancer type
carboplatin_autonomic_ganglia_lines              <- carboplatin_test$Tissue == 'autonomic_ganglia'
carboplatin_biliary_tract_lines                  <- carboplatin_test$Tissue == 'biliary_tract'
carboplatin_bone_lines                           <- carboplatin_test$Tissue == 'bone'
carboplatin_breast_lines                         <- carboplatin_test$Tissue == 'breast'
carboplatin_central_nervous_system_lines         <- carboplatin_test$Tissue == 'central_nervous_system'
carboplatin_endometrium_lines                    <- carboplatin_test$Tissue == 'endometrium'
carboplatin_kidney_lines                         <- carboplatin_test$Tissue == 'kidney'
carboplatin_large_intestine_lines                <- carboplatin_test$Tissue == 'large_intestine'
carboplatin_liver_lines                          <- carboplatin_test$Tissue == 'liver'
carboplatin_lung_lines                           <- carboplatin_test$Tissue == 'lung'
carboplatin_oesaphagus_lines                     <- carboplatin_test$Tissue == 'oesphagus'
carboplatin_ovary_lines                          <- carboplatin_test$Tissue == 'ovary'
carboplatin_pancreas_lines                       <- carboplatin_test$Tissue == 'pancreas'
carboplatin_pleura_lines                         <- carboplatin_test$Tissue == 'pleura'
carboplatin_prostate_lines                       <- carboplatin_test$Tissue == 'prostate'
carboplatin_salivary_gland_lines                 <- carboplatin_test$Tissue == 'salivary_gland'
carboplatin_skin_lines                           <- carboplatin_test$Tissue == 'skin'
carboplatin_soft_tissue_lines                    <- carboplatin_test$Tissue == 'soft_tissue'
carboplatin_stomach_lines                        <- carboplatin_test$Tissue == 'stomach'
carboplatin_thyroid_lines                        <- carboplatin_test$Tissue == 'thyroid'
carboplatin_upper_aerodigestive_tract_lines      <- carboplatin_test$Tissue == 'upper_aerodigestive_tract'
carboplatin_urinary_tract_lines                  <- carboplatin_test$Tissue == 'urinary_tract'

carboplatin_test_scaled_autonomic_ganglia <- carboplatin_rna_seq_test_scaled[carboplatin_autonomic_ganglia_lines, ]
carboplatin_test_autonomic_ganglia <- carboplatin_test[which(carboplatin_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
carboplatin_autonomic_ganglia_most_auc <- auc(carboplatin_test_autonomic_ganglia$most_sensitive, new_ic50)

carboplatin_test_scaled_biliary_tract <- carboplatin_rna_seq_test_scaled[carboplatin_biliary_tract_lines, ]
carboplatin_test_biliary_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_biliary_tract_most_auc <- auc(carboplatin_test_biliary_tract$most_sensitive, new_ic50)

carboplatin_test_scaled_bone <- carboplatin_rna_seq_test_scaled[carboplatin_bone_lines, ]
carboplatin_test_bone <- carboplatin_test[which(carboplatin_test$Tissue == 'bone'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
carboplatin_bone_most_auc <- auc(carboplatin_test_bone$most_sensitive, new_ic50)

carboplatin_test_scaled_breast <- carboplatin_rna_seq_test_scaled[carboplatin_breast_lines, ]
carboplatin_test_breast <- carboplatin_test[which(carboplatin_test$Tissue == 'breast'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
carboplatin_breast_most_auc <- auc(carboplatin_test_breast$most_sensitive, new_ic50)

carboplatin_test_scaled_central_nervous_system <- carboplatin_rna_seq_test_scaled[carboplatin_central_nervous_system_lines, ]
carboplatin_test_central_nervous_system <- carboplatin_test[which(carboplatin_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
carboplatin_central_nervous_system_most_auc <- auc(carboplatin_test_central_nervous_system$most_sensitive, new_ic50)

carboplatin_test_scaled_endometrium <- carboplatin_rna_seq_test_scaled[carboplatin_endometrium_lines, ]
carboplatin_test_endometrium <- carboplatin_test[which(carboplatin_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
carboplatin_endometrium_most_auc <- auc(carboplatin_test_endometrium$most_sensitive, new_ic50)

carboplatin_test_scaled_kidney <- carboplatin_rna_seq_test_scaled[carboplatin_kidney_lines, ]
carboplatin_test_kidney <- carboplatin_test[which(carboplatin_test$Tissue == 'kidney'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
carboplatin_kidney_most_auc <- auc(carboplatin_test_kidney$most_sensitive, new_ic50)

carboplatin_test_scaled_large_intestine <- carboplatin_rna_seq_test_scaled[carboplatin_large_intestine_lines, ]
carboplatin_test_large_intestine <- carboplatin_test[which(carboplatin_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
carboplatin_large_intestine_most_auc <- auc(carboplatin_test_large_intestine$most_sensitive, new_ic50)

carboplatin_test_scaled_liver <- carboplatin_rna_seq_test_scaled[carboplatin_liver_lines, ]
carboplatin_test_liver <- carboplatin_test[which(carboplatin_test$Tissue == 'liver'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
carboplatin_liver_most_auc <- auc(carboplatin_test_liver$most_sensitive, new_ic50)

carboplatin_test_scaled_lung <- carboplatin_rna_seq_test_scaled[carboplatin_lung_lines, ]
carboplatin_test_lung <- carboplatin_test[which(carboplatin_test$Tissue == 'lung'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
carboplatin_lung_most_auc <- auc(carboplatin_test_lung$most_sensitive, new_ic50)

carboplatin_test_scaled_oesaphagus <- carboplatin_rna_seq_test_scaled[carboplatin_oesaphagus_lines, ]
carboplatin_test_oesaphagus <- carboplatin_test[which(carboplatin_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
carboplatin_oesaphagus_most_auc <- auc(carboplatin_test_oesaphagus$most_sensitive, new_ic50)

carboplatin_test_scaled_ovary <- carboplatin_rna_seq_test_scaled[carboplatin_ovary_lines, ]
carboplatin_test_ovary <- carboplatin_test[which(carboplatin_test$Tissue == 'ovary'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
carboplatin_ovary_most_auc <- auc(carboplatin_test_ovary$most_sensitive, new_ic50)

carboplatin_test_scaled_pancreas <- carboplatin_rna_seq_test_scaled[carboplatin_pancreas_lines, ]
carboplatin_test_pancreas <- carboplatin_test[which(carboplatin_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
carboplatin_pancreas_most_auc <- auc(carboplatin_test_pancreas$most_sensitive, new_ic50)

carboplatin_test_scaled_pleura <- carboplatin_rna_seq_test_scaled[carboplatin_pleura_lines, ]
carboplatin_test_pleura <- carboplatin_test[which(carboplatin_test$Tissue == 'pleura'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
carboplatin_pleura_most_auc <- auc(carboplatin_test_pleura$most_sensitive, new_ic50)

carboplatin_test_scaled_prostate <- carboplatin_rna_seq_test_scaled[carboplatin_prostate_lines, ]
carboplatin_test_prostate <- carboplatin_test[which(carboplatin_test$Tissue == 'prostate'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
carboplatin_prostate_most_auc <- auc(carboplatin_test_prostate$most_sensitive, new_ic50)

carboplatin_test_scaled_salivary_gland <- carboplatin_rna_seq_test_scaled[carboplatin_salivary_gland_lines, ]
carboplatin_test_salivary_gland <- carboplatin_test[which(carboplatin_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
carboplatin_salivary_gland_most_auc <- auc(carboplatin_test_salivary_gland$most_sensitive, new_ic50)

carboplatin_test_scaled_skin <- carboplatin_rna_seq_test_scaled[carboplatin_skin_lines, ]
carboplatin_test_skin <- carboplatin_test[which(carboplatin_test$Tissue == 'skin'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
carboplatin_skin_most_auc <- auc(carboplatin_test_skin$most_sensitive, new_ic50)

carboplatin_test_scaled_soft_tissue <- carboplatin_rna_seq_test_scaled[carboplatin_soft_tissue_lines, ]
carboplatin_test_soft_tissue <- carboplatin_test[which(carboplatin_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
carboplatin_soft_tissue_most_auc <- auc(carboplatin_test_soft_tissue$most_sensitive, new_ic50)

carboplatin_test_scaled_stomach <- carboplatin_rna_seq_test_scaled[carboplatin_stomach_lines, ]
carboplatin_test_stomach <- carboplatin_test[which(carboplatin_test$Tissue == 'stomach'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
carboplatin_stomach_most_auc <- auc(carboplatin_test_stomach$most_sensitive, new_ic50)

carboplatin_test_scaled_thyroid <- carboplatin_rna_seq_test_scaled[carboplatin_thyroid_lines, ]
carboplatin_test_thyroid <- carboplatin_test[which(carboplatin_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
carboplatin_thyroid_most_auc <- auc(carboplatin_test_thyroid$most_sensitive, new_ic50)

carboplatin_test_scaled_upper_aerodigestive_tract <- carboplatin_rna_seq_test_scaled[carboplatin_upper_aerodigestive_tract_lines, ]
carboplatin_test_upper_aerodigestive_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_upper_aerodigestive_tract_most_auc <- auc(carboplatin_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

carboplatin_test_scaled_urinary_tract <- carboplatin_rna_seq_test_scaled[carboplatin_urinary_tract_lines, ]
carboplatin_test_urinary_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_urinary_tract_most_auc <- auc(carboplatin_test_urinary_tract$most_sensitive, new_ic50)

carboplatin_most_auc <- c(carboplatin_autonomic_ganglia_most_auc, carboplatin_biliary_tract_most_auc, 
                          carboplatin_bone_most_auc, carboplatin_breast_most_auc, carboplatin_central_nervous_system_most_auc, 
                          carboplatin_endometrium_most_auc, carboplatin_kidney_most_auc, 
                          carboplatin_large_intestine_most_auc, carboplatin_liver_most_auc, 
                          carboplatin_lung_most_auc, carboplatin_oesaphagus_most_auc, 
                          carboplatin_ovary_most_auc, carboplatin_pancreas_most_auc, 
                          carboplatin_pleura_most_auc, carboplatin_prostate_most_auc, 
                          carboplatin_salivary_gland_most_auc, carboplatin_skin_most_auc, 
                          carboplatin_soft_tissue_most_auc, carboplatin_stomach_most_auc, 
                          carboplatin_thyroid_most_auc, carboplatin_upper_aerodigestive_tract_most_auc, 
                          carboplatin_urinary_tract_most_auc)

carboplatin_test_scaled_autonomic_ganglia <- carboplatin_rna_seq_test_scaled[carboplatin_autonomic_ganglia_lines, ]
carboplatin_test_autonomic_ganglia <- carboplatin_test[which(carboplatin_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
carboplatin_autonomic_ganglia_least_auc <- auc(carboplatin_test_autonomic_ganglia$least_sensitive, new_ic50)

carboplatin_test_scaled_biliary_tract <- carboplatin_rna_seq_test_scaled[carboplatin_biliary_tract_lines, ]
carboplatin_test_biliary_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_biliary_tract_least_auc <- auc(carboplatin_test_biliary_tract$least_sensitive, new_ic50)

carboplatin_test_scaled_bone <- carboplatin_rna_seq_test_scaled[carboplatin_bone_lines, ]
carboplatin_test_bone <- carboplatin_test[which(carboplatin_test$Tissue == 'bone'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
carboplatin_bone_least_auc <- auc(carboplatin_test_bone$least_sensitive, new_ic50)

carboplatin_test_scaled_breast <- carboplatin_rna_seq_test_scaled[carboplatin_breast_lines, ]
carboplatin_test_breast <- carboplatin_test[which(carboplatin_test$Tissue == 'breast'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
carboplatin_breast_least_auc <- auc(carboplatin_test_breast$least_sensitive, new_ic50)

carboplatin_test_scaled_central_nervous_system <- carboplatin_rna_seq_test_scaled[carboplatin_central_nervous_system_lines, ]
carboplatin_test_central_nervous_system <- carboplatin_test[which(carboplatin_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
carboplatin_central_nervous_system_least_auc <- auc(carboplatin_test_central_nervous_system$least_sensitive, new_ic50)

carboplatin_test_scaled_endometrium <- carboplatin_rna_seq_test_scaled[carboplatin_endometrium_lines, ]
carboplatin_test_endometrium <- carboplatin_test[which(carboplatin_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
carboplatin_endometrium_least_auc <- auc(carboplatin_test_endometrium$least_sensitive, new_ic50)

carboplatin_test_scaled_kidney <- carboplatin_rna_seq_test_scaled[carboplatin_kidney_lines, ]
carboplatin_test_kidney <- carboplatin_test[which(carboplatin_test$Tissue == 'kidney'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
carboplatin_kidney_least_auc <- auc(carboplatin_test_kidney$least_sensitive, new_ic50)

carboplatin_test_scaled_large_intestine <- carboplatin_rna_seq_test_scaled[carboplatin_large_intestine_lines, ]
carboplatin_test_large_intestine <- carboplatin_test[which(carboplatin_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
carboplatin_large_intestine_least_auc <- auc(carboplatin_test_large_intestine$least_sensitive, new_ic50)

carboplatin_test_scaled_liver <- carboplatin_rna_seq_test_scaled[carboplatin_liver_lines, ]
carboplatin_test_liver <- carboplatin_test[which(carboplatin_test$Tissue == 'liver'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
carboplatin_liver_least_auc <- auc(carboplatin_test_liver$least_sensitive, new_ic50)

carboplatin_test_scaled_lung <- carboplatin_rna_seq_test_scaled[carboplatin_lung_lines, ]
carboplatin_test_lung <- carboplatin_test[which(carboplatin_test$Tissue == 'lung'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
carboplatin_lung_least_auc <- auc(carboplatin_test_lung$least_sensitive, new_ic50)

carboplatin_test_scaled_oesaphagus <- carboplatin_rna_seq_test_scaled[carboplatin_oesaphagus_lines, ]
carboplatin_test_oesaphagus <- carboplatin_test[which(carboplatin_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
carboplatin_oesaphagus_least_auc <- auc(carboplatin_test_oesaphagus$least_sensitive, new_ic50)

carboplatin_test_scaled_ovary <- carboplatin_rna_seq_test_scaled[carboplatin_ovary_lines, ]
carboplatin_test_ovary <- carboplatin_test[which(carboplatin_test$Tissue == 'ovary'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
carboplatin_ovary_least_auc <- auc(carboplatin_test_ovary$least_sensitive, new_ic50)

carboplatin_test_scaled_pancreas <- carboplatin_rna_seq_test_scaled[carboplatin_pancreas_lines, ]
carboplatin_test_pancreas <- carboplatin_test[which(carboplatin_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
carboplatin_pancreas_least_auc <- auc(carboplatin_test_pancreas$least_sensitive, new_ic50)

carboplatin_test_scaled_pleura <- carboplatin_rna_seq_test_scaled[carboplatin_pleura_lines, ]
carboplatin_test_pleura <- carboplatin_test[which(carboplatin_test$Tissue == 'pleura'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
carboplatin_pleura_least_auc <- auc(carboplatin_test_pleura$least_sensitive, new_ic50)

carboplatin_test_scaled_prostate <- carboplatin_rna_seq_test_scaled[carboplatin_prostate_lines, ]
carboplatin_test_prostate <- carboplatin_test[which(carboplatin_test$Tissue == 'prostate'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
carboplatin_prostate_least_auc <- auc(carboplatin_test_prostate$least_sensitive, new_ic50)

carboplatin_test_scaled_salivary_gland <- carboplatin_rna_seq_test_scaled[carboplatin_salivary_gland_lines, ]
carboplatin_test_salivary_gland <- carboplatin_test[which(carboplatin_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
carboplatin_salivary_gland_least_auc <- auc(carboplatin_test_salivary_gland$least_sensitive, new_ic50)

carboplatin_test_scaled_skin <- carboplatin_rna_seq_test_scaled[carboplatin_skin_lines, ]
carboplatin_test_skin <- carboplatin_test[which(carboplatin_test$Tissue == 'skin'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
carboplatin_skin_least_auc <- auc(carboplatin_test_skin$least_sensitive, new_ic50)

carboplatin_test_scaled_soft_tissue <- carboplatin_rna_seq_test_scaled[carboplatin_soft_tissue_lines, ]
carboplatin_test_soft_tissue <- carboplatin_test[which(carboplatin_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
carboplatin_soft_tissue_least_auc <- auc(carboplatin_test_soft_tissue$least_sensitive, new_ic50)

carboplatin_test_scaled_stomach <- carboplatin_rna_seq_test_scaled[carboplatin_stomach_lines, ]
carboplatin_test_stomach <- carboplatin_test[which(carboplatin_test$Tissue == 'stomach'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
carboplatin_stomach_least_auc <- auc(carboplatin_test_stomach$least_sensitive, new_ic50)

carboplatin_test_scaled_thyroid <- carboplatin_rna_seq_test_scaled[carboplatin_thyroid_lines, ]
carboplatin_test_thyroid <- carboplatin_test[which(carboplatin_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
carboplatin_thyroid_least_auc <- auc(carboplatin_test_thyroid$least_sensitive, new_ic50)

carboplatin_test_scaled_upper_aerodigestive_tract <- carboplatin_rna_seq_test_scaled[carboplatin_upper_aerodigestive_tract_lines, ]
carboplatin_test_upper_aerodigestive_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_upper_aerodigestive_tract_least_auc <- auc(carboplatin_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

carboplatin_test_scaled_urinary_tract <- carboplatin_rna_seq_test_scaled[carboplatin_urinary_tract_lines, ]
carboplatin_test_urinary_tract <- carboplatin_test[which(carboplatin_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
carboplatin_urinary_tract_least_auc <- auc(carboplatin_test_urinary_tract$least_sensitive, new_ic50)

carboplatin_least_auc <- c(carboplatin_autonomic_ganglia_least_auc, carboplatin_biliary_tract_least_auc, 
                          carboplatin_bone_least_auc, carboplatin_breast_least_auc, carboplatin_central_nervous_system_least_auc, 
                          carboplatin_endometrium_least_auc, carboplatin_kidney_least_auc, 
                          carboplatin_large_intestine_least_auc, carboplatin_liver_least_auc, 
                          carboplatin_lung_least_auc, carboplatin_oesaphagus_least_auc, 
                          carboplatin_ovary_least_auc, carboplatin_pancreas_least_auc, 
                          carboplatin_pleura_least_auc, carboplatin_prostate_least_auc, 
                          carboplatin_salivary_gland_least_auc, carboplatin_skin_least_auc, 
                          carboplatin_soft_tissue_least_auc, carboplatin_stomach_least_auc, 
                          carboplatin_thyroid_least_auc, carboplatin_upper_aerodigestive_tract_least_auc, 
                          carboplatin_urinary_tract_least_auc)

cyclophosphamide_autonomic_ganglia_lines              <- cyclophosphamide_test$Tissue == 'autonomic_ganglia'
cyclophosphamide_biliary_tract_lines                  <- cyclophosphamide_test$Tissue == 'biliary_tract'
cyclophosphamide_bone_lines                           <- cyclophosphamide_test$Tissue == 'bone'
cyclophosphamide_breast_lines                         <- cyclophosphamide_test$Tissue == 'breast'
cyclophosphamide_central_nervous_system_lines         <- cyclophosphamide_test$Tissue == 'central_nervous_system'
cyclophosphamide_endometrium_lines                    <- cyclophosphamide_test$Tissue == 'endometrium'
cyclophosphamide_kidney_lines                         <- cyclophosphamide_test$Tissue == 'kidney'
cyclophosphamide_large_intestine_lines                <- cyclophosphamide_test$Tissue == 'large_intestine'
cyclophosphamide_liver_lines                          <- cyclophosphamide_test$Tissue == 'liver'
cyclophosphamide_lung_lines                           <- cyclophosphamide_test$Tissue == 'lung'
cyclophosphamide_oesaphagus_lines                     <- cyclophosphamide_test$Tissue == 'oesphagus'
cyclophosphamide_ovary_lines                          <- cyclophosphamide_test$Tissue == 'ovary'
cyclophosphamide_pancreas_lines                       <- cyclophosphamide_test$Tissue == 'pancreas'
cyclophosphamide_pleura_lines                         <- cyclophosphamide_test$Tissue == 'pleura'
cyclophosphamide_prostate_lines                       <- cyclophosphamide_test$Tissue == 'prostate'
cyclophosphamide_salivary_gland_lines                 <- cyclophosphamide_test$Tissue == 'salivary_gland'
cyclophosphamide_skin_lines                           <- cyclophosphamide_test$Tissue == 'skin'
cyclophosphamide_soft_tissue_lines                    <- cyclophosphamide_test$Tissue == 'soft_tissue'
cyclophosphamide_stomach_lines                        <- cyclophosphamide_test$Tissue == 'stomach'
cyclophosphamide_thyroid_lines                        <- cyclophosphamide_test$Tissue == 'thyroid'
cyclophosphamide_upper_aerodigestive_tract_lines      <- cyclophosphamide_test$Tissue == 'upper_aerodigestive_tract'
cyclophosphamide_urinary_tract_lines                  <- cyclophosphamide_test$Tissue == 'urinary_tract'

cyclophosphamide_test_scaled_autonomic_ganglia <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_autonomic_ganglia_lines, ]
cyclophosphamide_test_autonomic_ganglia <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_autonomic_ganglia_most_auc <- auc(cyclophosphamide_test_autonomic_ganglia$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_biliary_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_biliary_tract_lines, ]
cyclophosphamide_test_biliary_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_biliary_tract_most_auc <- auc(cyclophosphamide_test_biliary_tract$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_bone <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_bone_lines, ]
cyclophosphamide_test_bone <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'bone'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_bone_most_auc <- auc(cyclophosphamide_test_bone$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_breast <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_breast_lines, ]
cyclophosphamide_test_breast <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'breast'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_breast_most_auc <- auc(cyclophosphamide_test_breast$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_central_nervous_system <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_central_nervous_system_lines, ]
cyclophosphamide_test_central_nervous_system <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_central_nervous_system_most_auc <- auc(cyclophosphamide_test_central_nervous_system$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_endometrium <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_endometrium_lines, ]
cyclophosphamide_test_endometrium <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_endometrium_most_auc <- auc(cyclophosphamide_test_endometrium$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_kidney <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_kidney_lines, ]
cyclophosphamide_test_kidney <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'kidney'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_kidney_most_auc <- auc(cyclophosphamide_test_kidney$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_large_intestine <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_large_intestine_lines, ]
cyclophosphamide_test_large_intestine <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_large_intestine_most_auc <- auc(cyclophosphamide_test_large_intestine$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_liver <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_liver_lines, ]
cyclophosphamide_test_liver <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'liver'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_liver_most_auc <- auc(cyclophosphamide_test_liver$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_lung <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_lung_lines, ]
cyclophosphamide_test_lung <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'lung'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_lung_most_auc <- auc(cyclophosphamide_test_lung$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_oesaphagus <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_oesaphagus_lines, ]
cyclophosphamide_test_oesaphagus <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_oesaphagus_most_auc <- auc(cyclophosphamide_test_oesaphagus$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_ovary <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_ovary_lines, ]
cyclophosphamide_test_ovary <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'ovary'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_ovary_most_auc <- auc(cyclophosphamide_test_ovary$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_pancreas <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_pancreas_lines, ]
cyclophosphamide_test_pancreas <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_pancreas_most_auc <- auc(cyclophosphamide_test_pancreas$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_pleura <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_pleura_lines, ]
cyclophosphamide_test_pleura <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'pleura'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_pleura_most_auc <- auc(cyclophosphamide_test_pleura$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_prostate <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_prostate_lines, ]
cyclophosphamide_test_prostate <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'prostate'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_prostate_most_auc <- auc(cyclophosphamide_test_prostate$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_salivary_gland <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_salivary_gland_lines, ]
cyclophosphamide_test_salivary_gland <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_salivary_gland_most_auc <- auc(cyclophosphamide_test_salivary_gland$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_skin <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_skin_lines, ]
cyclophosphamide_test_skin <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'skin'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_skin_most_auc <- auc(cyclophosphamide_test_skin$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_soft_tissue <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_soft_tissue_lines, ]
cyclophosphamide_test_soft_tissue <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_soft_tissue_most_auc <- auc(cyclophosphamide_test_soft_tissue$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_stomach <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_stomach_lines, ]
cyclophosphamide_test_stomach <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'stomach'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_stomach_most_auc <- auc(cyclophosphamide_test_stomach$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_thyroid <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_thyroid_lines, ]
cyclophosphamide_test_thyroid <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_thyroid_most_auc <- auc(cyclophosphamide_test_thyroid$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_upper_aerodigestive_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_upper_aerodigestive_tract_lines, ]
cyclophosphamide_test_upper_aerodigestive_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_upper_aerodigestive_tract_most_auc <- auc(cyclophosphamide_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

cyclophosphamide_test_scaled_urinary_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_urinary_tract_lines, ]
cyclophosphamide_test_urinary_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_urinary_tract_most_auc <- auc(cyclophosphamide_test_urinary_tract$most_sensitive, new_ic50)

cyclophosphamide_most_auc <- c(cyclophosphamide_autonomic_ganglia_most_auc, cyclophosphamide_biliary_tract_most_auc, 
                          cyclophosphamide_bone_most_auc, cyclophosphamide_breast_most_auc, cyclophosphamide_central_nervous_system_most_auc, 
                          cyclophosphamide_endometrium_most_auc, cyclophosphamide_kidney_most_auc, 
                          cyclophosphamide_large_intestine_most_auc, cyclophosphamide_liver_most_auc, 
                          cyclophosphamide_lung_most_auc, cyclophosphamide_oesaphagus_most_auc, 
                          cyclophosphamide_ovary_most_auc, cyclophosphamide_pancreas_most_auc, 
                          cyclophosphamide_pleura_most_auc, cyclophosphamide_prostate_most_auc, 
                          cyclophosphamide_salivary_gland_most_auc, cyclophosphamide_skin_most_auc, 
                          cyclophosphamide_soft_tissue_most_auc, cyclophosphamide_stomach_most_auc, 
                          cyclophosphamide_thyroid_most_auc, cyclophosphamide_upper_aerodigestive_tract_most_auc, 
                          cyclophosphamide_urinary_tract_most_auc)

cyclophosphamide_test_scaled_autonomic_ganglia <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_autonomic_ganglia_lines, ]
cyclophosphamide_test_autonomic_ganglia <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_autonomic_ganglia_least_auc <- auc(cyclophosphamide_test_autonomic_ganglia$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_biliary_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_biliary_tract_lines, ]
cyclophosphamide_test_biliary_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_biliary_tract_least_auc <- auc(cyclophosphamide_test_biliary_tract$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_bone <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_bone_lines, ]
cyclophosphamide_test_bone <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'bone'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_bone_least_auc <- auc(cyclophosphamide_test_bone$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_breast <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_breast_lines, ]
cyclophosphamide_test_breast <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'breast'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_breast_least_auc <- auc(cyclophosphamide_test_breast$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_central_nervous_system <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_central_nervous_system_lines, ]
cyclophosphamide_test_central_nervous_system <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_central_nervous_system_least_auc <- auc(cyclophosphamide_test_central_nervous_system$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_endometrium <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_endometrium_lines, ]
cyclophosphamide_test_endometrium <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_endometrium_least_auc <- auc(cyclophosphamide_test_endometrium$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_kidney <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_kidney_lines, ]
cyclophosphamide_test_kidney <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'kidney'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_kidney_least_auc <- auc(cyclophosphamide_test_kidney$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_large_intestine <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_large_intestine_lines, ]
cyclophosphamide_test_large_intestine <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_large_intestine_least_auc <- auc(cyclophosphamide_test_large_intestine$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_liver <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_liver_lines, ]
cyclophosphamide_test_liver <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'liver'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_liver_least_auc <- auc(cyclophosphamide_test_liver$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_lung <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_lung_lines, ]
cyclophosphamide_test_lung <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'lung'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_lung_least_auc <- auc(cyclophosphamide_test_lung$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_oesaphagus <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_oesaphagus_lines, ]
cyclophosphamide_test_oesaphagus <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_oesaphagus_least_auc <- auc(cyclophosphamide_test_oesaphagus$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_ovary <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_ovary_lines, ]
cyclophosphamide_test_ovary <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'ovary'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_ovary_least_auc <- auc(cyclophosphamide_test_ovary$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_pancreas <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_pancreas_lines, ]
cyclophosphamide_test_pancreas <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_pancreas_least_auc <- auc(cyclophosphamide_test_pancreas$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_pleura <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_pleura_lines, ]
cyclophosphamide_test_pleura <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'pleura'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_pleura_least_auc <- auc(cyclophosphamide_test_pleura$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_prostate <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_prostate_lines, ]
cyclophosphamide_test_prostate <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'prostate'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_prostate_least_auc <- auc(cyclophosphamide_test_prostate$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_salivary_gland <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_salivary_gland_lines, ]
cyclophosphamide_test_salivary_gland <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_salivary_gland_least_auc <- auc(cyclophosphamide_test_salivary_gland$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_skin <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_skin_lines, ]
cyclophosphamide_test_skin <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'skin'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_skin_least_auc <- auc(cyclophosphamide_test_skin$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_soft_tissue <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_soft_tissue_lines, ]
cyclophosphamide_test_soft_tissue <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_soft_tissue_least_auc <- auc(cyclophosphamide_test_soft_tissue$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_stomach <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_stomach_lines, ]
cyclophosphamide_test_stomach <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'stomach'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_stomach_least_auc <- auc(cyclophosphamide_test_stomach$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_thyroid <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_thyroid_lines, ]
cyclophosphamide_test_thyroid <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_thyroid_least_auc <- auc(cyclophosphamide_test_thyroid$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_upper_aerodigestive_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_upper_aerodigestive_tract_lines, ]
cyclophosphamide_test_upper_aerodigestive_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_upper_aerodigestive_tract_least_auc <- auc(cyclophosphamide_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

cyclophosphamide_test_scaled_urinary_tract <- cyclophosphamide_rna_seq_test_scaled[cyclophosphamide_urinary_tract_lines, ]
cyclophosphamide_test_urinary_tract <- cyclophosphamide_test[which(cyclophosphamide_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
cyclophosphamide_urinary_tract_least_auc <- auc(cyclophosphamide_test_urinary_tract$least_sensitive, new_ic50)

cyclophosphamide_least_auc <- c(cyclophosphamide_autonomic_ganglia_least_auc, cyclophosphamide_biliary_tract_least_auc, 
                           cyclophosphamide_bone_least_auc, cyclophosphamide_breast_least_auc, cyclophosphamide_central_nervous_system_least_auc, 
                           cyclophosphamide_endometrium_least_auc, cyclophosphamide_kidney_least_auc, 
                           cyclophosphamide_large_intestine_least_auc, cyclophosphamide_liver_least_auc, 
                           cyclophosphamide_lung_least_auc, cyclophosphamide_oesaphagus_least_auc, 
                           cyclophosphamide_ovary_least_auc, cyclophosphamide_pancreas_least_auc, 
                           cyclophosphamide_pleura_least_auc, cyclophosphamide_prostate_least_auc, 
                           cyclophosphamide_salivary_gland_least_auc, cyclophosphamide_skin_least_auc, 
                           cyclophosphamide_soft_tissue_least_auc, cyclophosphamide_stomach_least_auc, 
                           cyclophosphamide_thyroid_least_auc, cyclophosphamide_upper_aerodigestive_tract_least_auc, 
                           cyclophosphamide_urinary_tract_least_auc)

docetaxel_autonomic_ganglia_lines              <- docetaxel_test$Tissue == 'autonomic_ganglia'
docetaxel_biliary_tract_lines                  <- docetaxel_test$Tissue == 'biliary_tract'
docetaxel_bone_lines                           <- docetaxel_test$Tissue == 'bone'
docetaxel_breast_lines                         <- docetaxel_test$Tissue == 'breast'
docetaxel_central_nervous_system_lines         <- docetaxel_test$Tissue == 'central_nervous_system'
docetaxel_endometrium_lines                    <- docetaxel_test$Tissue == 'endometrium'
docetaxel_kidney_lines                         <- docetaxel_test$Tissue == 'kidney'
docetaxel_large_intestine_lines                <- docetaxel_test$Tissue == 'large_intestine'
docetaxel_liver_lines                          <- docetaxel_test$Tissue == 'liver'
docetaxel_lung_lines                           <- docetaxel_test$Tissue == 'lung'
docetaxel_oesaphagus_lines                     <- docetaxel_test$Tissue == 'oesphagus'
docetaxel_ovary_lines                          <- docetaxel_test$Tissue == 'ovary'
docetaxel_pancreas_lines                       <- docetaxel_test$Tissue == 'pancreas'
docetaxel_pleura_lines                         <- docetaxel_test$Tissue == 'pleura'
docetaxel_prostate_lines                       <- docetaxel_test$Tissue == 'prostate'
docetaxel_salivary_gland_lines                 <- docetaxel_test$Tissue == 'salivary_gland'
docetaxel_skin_lines                           <- docetaxel_test$Tissue == 'skin'
docetaxel_soft_tissue_lines                    <- docetaxel_test$Tissue == 'soft_tissue'
docetaxel_stomach_lines                        <- docetaxel_test$Tissue == 'stomach'
docetaxel_thyroid_lines                        <- docetaxel_test$Tissue == 'thyroid'
docetaxel_upper_aerodigestive_tract_lines      <- docetaxel_test$Tissue == 'upper_aerodigestive_tract'
docetaxel_urinary_tract_lines                  <- docetaxel_test$Tissue == 'urinary_tract'

docetaxel_test_scaled_autonomic_ganglia <- docetaxel_rna_seq_test_scaled[docetaxel_autonomic_ganglia_lines, ]
docetaxel_test_autonomic_ganglia <- docetaxel_test[which(docetaxel_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
docetaxel_autonomic_ganglia_most_auc <- auc(docetaxel_test_autonomic_ganglia$most_sensitive, new_ic50)

docetaxel_test_scaled_biliary_tract <- docetaxel_rna_seq_test_scaled[docetaxel_biliary_tract_lines, ]
docetaxel_test_biliary_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
docetaxel_biliary_tract_most_auc <- auc(docetaxel_test_biliary_tract$most_sensitive, new_ic50)

docetaxel_test_scaled_bone <- docetaxel_rna_seq_test_scaled[docetaxel_bone_lines, ]
docetaxel_test_bone <- docetaxel_test[which(docetaxel_test$Tissue == 'bone'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
docetaxel_bone_most_auc <- auc(docetaxel_test_bone$most_sensitive, new_ic50)

docetaxel_test_scaled_breast <- docetaxel_rna_seq_test_scaled[docetaxel_breast_lines, ]
docetaxel_test_breast <- docetaxel_test[which(docetaxel_test$Tissue == 'breast'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
docetaxel_breast_most_auc <- auc(docetaxel_test_breast$most_sensitive, new_ic50)

docetaxel_test_scaled_central_nervous_system <- docetaxel_rna_seq_test_scaled[docetaxel_central_nervous_system_lines, ]
docetaxel_test_central_nervous_system <- docetaxel_test[which(docetaxel_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
docetaxel_central_nervous_system_most_auc <- auc(docetaxel_test_central_nervous_system$most_sensitive, new_ic50)

docetaxel_test_scaled_endometrium <- docetaxel_rna_seq_test_scaled[docetaxel_endometrium_lines, ]
docetaxel_test_endometrium <- docetaxel_test[which(docetaxel_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
docetaxel_endometrium_most_auc <- auc(docetaxel_test_endometrium$most_sensitive, new_ic50)

docetaxel_test_scaled_kidney <- docetaxel_rna_seq_test_scaled[docetaxel_kidney_lines, ]
docetaxel_test_kidney <- docetaxel_test[which(docetaxel_test$Tissue == 'kidney'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
docetaxel_kidney_most_auc <- auc(docetaxel_test_kidney$most_sensitive, new_ic50)

docetaxel_test_scaled_large_intestine <- docetaxel_rna_seq_test_scaled[docetaxel_large_intestine_lines, ]
docetaxel_test_large_intestine <- docetaxel_test[which(docetaxel_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
docetaxel_large_intestine_most_auc <- auc(docetaxel_test_large_intestine$most_sensitive, new_ic50)

docetaxel_test_scaled_liver <- docetaxel_rna_seq_test_scaled[docetaxel_liver_lines, ]
docetaxel_test_liver <- docetaxel_test[which(docetaxel_test$Tissue == 'liver'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
docetaxel_liver_most_auc <- auc(docetaxel_test_liver$most_sensitive, new_ic50)

docetaxel_test_scaled_lung <- docetaxel_rna_seq_test_scaled[docetaxel_lung_lines, ]
docetaxel_test_lung <- docetaxel_test[which(docetaxel_test$Tissue == 'lung'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
docetaxel_lung_most_auc <- auc(docetaxel_test_lung$most_sensitive, new_ic50)

docetaxel_test_scaled_oesaphagus <- docetaxel_rna_seq_test_scaled[docetaxel_oesaphagus_lines, ]
docetaxel_test_oesaphagus <- docetaxel_test[which(docetaxel_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
docetaxel_oesaphagus_most_auc <- auc(docetaxel_test_oesaphagus$most_sensitive, new_ic50)

docetaxel_test_scaled_ovary <- docetaxel_rna_seq_test_scaled[docetaxel_ovary_lines, ]
docetaxel_test_ovary <- docetaxel_test[which(docetaxel_test$Tissue == 'ovary'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
docetaxel_ovary_most_auc <- auc(docetaxel_test_ovary$most_sensitive, new_ic50)

docetaxel_test_scaled_pancreas <- docetaxel_rna_seq_test_scaled[docetaxel_pancreas_lines, ]
docetaxel_test_pancreas <- docetaxel_test[which(docetaxel_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
docetaxel_pancreas_most_auc <- auc(docetaxel_test_pancreas$most_sensitive, new_ic50)

docetaxel_test_scaled_pleura <- docetaxel_rna_seq_test_scaled[docetaxel_pleura_lines, ]
docetaxel_test_pleura <- docetaxel_test[which(docetaxel_test$Tissue == 'pleura'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
docetaxel_pleura_most_auc <- auc(docetaxel_test_pleura$most_sensitive, new_ic50)

docetaxel_test_scaled_prostate <- docetaxel_rna_seq_test_scaled[docetaxel_prostate_lines, ]
docetaxel_test_prostate <- docetaxel_test[which(docetaxel_test$Tissue == 'prostate'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
docetaxel_prostate_most_auc <- auc(docetaxel_test_prostate$most_sensitive, new_ic50)

docetaxel_test_scaled_salivary_gland <- docetaxel_rna_seq_test_scaled[docetaxel_salivary_gland_lines, ]
docetaxel_test_salivary_gland <- docetaxel_test[which(docetaxel_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
docetaxel_salivary_gland_most_auc <- auc(docetaxel_test_salivary_gland$most_sensitive, new_ic50)

docetaxel_test_scaled_skin <- docetaxel_rna_seq_test_scaled[docetaxel_skin_lines, ]
docetaxel_test_skin <- docetaxel_test[which(docetaxel_test$Tissue == 'skin'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
docetaxel_skin_most_auc <- auc(docetaxel_test_skin$most_sensitive, new_ic50)

docetaxel_test_scaled_soft_tissue <- docetaxel_rna_seq_test_scaled[docetaxel_soft_tissue_lines, ]
docetaxel_test_soft_tissue <- docetaxel_test[which(docetaxel_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
docetaxel_soft_tissue_most_auc <- auc(docetaxel_test_soft_tissue$most_sensitive, new_ic50)

docetaxel_test_scaled_stomach <- docetaxel_rna_seq_test_scaled[docetaxel_stomach_lines, ]
docetaxel_test_stomach <- docetaxel_test[which(docetaxel_test$Tissue == 'stomach'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
docetaxel_stomach_most_auc <- auc(docetaxel_test_stomach$most_sensitive, new_ic50)

docetaxel_test_scaled_thyroid <- docetaxel_rna_seq_test_scaled[docetaxel_thyroid_lines, ]
docetaxel_test_thyroid <- docetaxel_test[which(docetaxel_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
docetaxel_thyroid_most_auc <- auc(docetaxel_test_thyroid$most_sensitive, new_ic50)

docetaxel_test_scaled_upper_aerodigestive_tract <- docetaxel_rna_seq_test_scaled[docetaxel_upper_aerodigestive_tract_lines, ]
docetaxel_test_upper_aerodigestive_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
docetaxel_upper_aerodigestive_tract_most_auc <- auc(docetaxel_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

docetaxel_test_scaled_urinary_tract <- docetaxel_rna_seq_test_scaled[docetaxel_urinary_tract_lines, ]
docetaxel_test_urinary_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
docetaxel_urinary_tract_most_auc <- auc(docetaxel_test_urinary_tract$most_sensitive, new_ic50)

docetaxel_most_auc <- c(docetaxel_autonomic_ganglia_most_auc, docetaxel_biliary_tract_most_auc, 
                          docetaxel_bone_most_auc, docetaxel_breast_most_auc, docetaxel_central_nervous_system_most_auc, 
                          docetaxel_endometrium_most_auc, docetaxel_kidney_most_auc, 
                          docetaxel_large_intestine_most_auc, docetaxel_liver_most_auc, 
                          docetaxel_lung_most_auc, docetaxel_oesaphagus_most_auc, 
                          docetaxel_ovary_most_auc, docetaxel_pancreas_most_auc, 
                          docetaxel_pleura_most_auc, docetaxel_prostate_most_auc, 
                          docetaxel_salivary_gland_most_auc, docetaxel_skin_most_auc, 
                          docetaxel_soft_tissue_most_auc, docetaxel_stomach_most_auc, 
                          docetaxel_thyroid_most_auc, docetaxel_upper_aerodigestive_tract_most_auc, 
                          docetaxel_urinary_tract_most_auc)

docetaxel_test_scaled_autonomic_ganglia <- docetaxel_rna_seq_test_scaled[docetaxel_autonomic_ganglia_lines, ]
docetaxel_test_autonomic_ganglia <- docetaxel_test[which(docetaxel_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_autonomic_ganglia), s = 'lambda.min', interval = 'conf')
docetaxel_autonomic_ganglia_least_auc <- auc(docetaxel_test_autonomic_ganglia$least_sensitive, new_ic50)

docetaxel_test_scaled_biliary_tract <- docetaxel_rna_seq_test_scaled[docetaxel_biliary_tract_lines, ]
docetaxel_test_biliary_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_biliary_tract), s = 'lambda.min', interval = 'conf')
docetaxel_biliary_tract_least_auc <- auc(docetaxel_test_biliary_tract$least_sensitive, new_ic50)

docetaxel_test_scaled_bone <- docetaxel_rna_seq_test_scaled[docetaxel_bone_lines, ]
docetaxel_test_bone <- docetaxel_test[which(docetaxel_test$Tissue == 'bone'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_bone), s = 'lambda.min', interval = 'conf')
docetaxel_bone_least_auc <- auc(docetaxel_test_bone$least_sensitive, new_ic50)

docetaxel_test_scaled_breast <- docetaxel_rna_seq_test_scaled[docetaxel_breast_lines, ]
docetaxel_test_breast <- docetaxel_test[which(docetaxel_test$Tissue == 'breast'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_breast), s = 'lambda.min', interval = 'conf')
docetaxel_breast_least_auc <- auc(docetaxel_test_breast$least_sensitive, new_ic50)

docetaxel_test_scaled_central_nervous_system <- docetaxel_rna_seq_test_scaled[docetaxel_central_nervous_system_lines, ]
docetaxel_test_central_nervous_system <- docetaxel_test[which(docetaxel_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_central_nervous_system), s = 'lambda.min', interval = 'conf')
docetaxel_central_nervous_system_least_auc <- auc(docetaxel_test_central_nervous_system$least_sensitive, new_ic50)

docetaxel_test_scaled_endometrium <- docetaxel_rna_seq_test_scaled[docetaxel_endometrium_lines, ]
docetaxel_test_endometrium <- docetaxel_test[which(docetaxel_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_endometrium), s = 'lambda.min', interval = 'conf')
docetaxel_endometrium_least_auc <- auc(docetaxel_test_endometrium$least_sensitive, new_ic50)

docetaxel_test_scaled_kidney <- docetaxel_rna_seq_test_scaled[docetaxel_kidney_lines, ]
docetaxel_test_kidney <- docetaxel_test[which(docetaxel_test$Tissue == 'kidney'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_kidney), s = 'lambda.min', interval = 'conf')
docetaxel_kidney_least_auc <- auc(docetaxel_test_kidney$least_sensitive, new_ic50)

docetaxel_test_scaled_large_intestine <- docetaxel_rna_seq_test_scaled[docetaxel_large_intestine_lines, ]
docetaxel_test_large_intestine <- docetaxel_test[which(docetaxel_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_large_intestine), s = 'lambda.min', interval = 'conf')
docetaxel_large_intestine_least_auc <- auc(docetaxel_test_large_intestine$least_sensitive, new_ic50)

docetaxel_test_scaled_liver <- docetaxel_rna_seq_test_scaled[docetaxel_liver_lines, ]
docetaxel_test_liver <- docetaxel_test[which(docetaxel_test$Tissue == 'liver'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_liver), s = 'lambda.min', interval = 'conf')
docetaxel_liver_least_auc <- auc(docetaxel_test_liver$least_sensitive, new_ic50)

docetaxel_test_scaled_lung <- docetaxel_rna_seq_test_scaled[docetaxel_lung_lines, ]
docetaxel_test_lung <- docetaxel_test[which(docetaxel_test$Tissue == 'lung'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_lung), s = 'lambda.min', interval = 'conf')
docetaxel_lung_least_auc <- auc(docetaxel_test_lung$least_sensitive, new_ic50)

docetaxel_test_scaled_oesaphagus <- docetaxel_rna_seq_test_scaled[docetaxel_oesaphagus_lines, ]
docetaxel_test_oesaphagus <- docetaxel_test[which(docetaxel_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_oesaphagus), s = 'lambda.min', interval = 'conf')
docetaxel_oesaphagus_least_auc <- auc(docetaxel_test_oesaphagus$least_sensitive, new_ic50)

docetaxel_test_scaled_ovary <- docetaxel_rna_seq_test_scaled[docetaxel_ovary_lines, ]
docetaxel_test_ovary <- docetaxel_test[which(docetaxel_test$Tissue == 'ovary'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_ovary), s = 'lambda.min', interval = 'conf')
docetaxel_ovary_least_auc <- auc(docetaxel_test_ovary$least_sensitive, new_ic50)

docetaxel_test_scaled_pancreas <- docetaxel_rna_seq_test_scaled[docetaxel_pancreas_lines, ]
docetaxel_test_pancreas <- docetaxel_test[which(docetaxel_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_pancreas), s = 'lambda.min', interval = 'conf')
docetaxel_pancreas_least_auc <- auc(docetaxel_test_pancreas$least_sensitive, new_ic50)

docetaxel_test_scaled_pleura <- docetaxel_rna_seq_test_scaled[docetaxel_pleura_lines, ]
docetaxel_test_pleura <- docetaxel_test[which(docetaxel_test$Tissue == 'pleura'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_pleura), s = 'lambda.min', interval = 'conf')
docetaxel_pleura_least_auc <- auc(docetaxel_test_pleura$least_sensitive, new_ic50)

docetaxel_test_scaled_prostate <- docetaxel_rna_seq_test_scaled[docetaxel_prostate_lines, ]
docetaxel_test_prostate <- docetaxel_test[which(docetaxel_test$Tissue == 'prostate'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_prostate), s = 'lambda.min', interval = 'conf')
docetaxel_prostate_least_auc <- auc(docetaxel_test_prostate$least_sensitive, new_ic50)

docetaxel_test_scaled_salivary_gland <- docetaxel_rna_seq_test_scaled[docetaxel_salivary_gland_lines, ]
docetaxel_test_salivary_gland <- docetaxel_test[which(docetaxel_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_salivary_gland), s = 'lambda.min', interval = 'conf')
docetaxel_salivary_gland_least_auc <- auc(docetaxel_test_salivary_gland$least_sensitive, new_ic50)

docetaxel_test_scaled_skin <- docetaxel_rna_seq_test_scaled[docetaxel_skin_lines, ]
docetaxel_test_skin <- docetaxel_test[which(docetaxel_test$Tissue == 'skin'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_skin), s = 'lambda.min', interval = 'conf')
docetaxel_skin_least_auc <- auc(docetaxel_test_skin$least_sensitive, new_ic50)

docetaxel_test_scaled_soft_tissue <- docetaxel_rna_seq_test_scaled[docetaxel_soft_tissue_lines, ]
docetaxel_test_soft_tissue <- docetaxel_test[which(docetaxel_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_soft_tissue), s = 'lambda.min', interval = 'conf')
docetaxel_soft_tissue_least_auc <- auc(docetaxel_test_soft_tissue$least_sensitive, new_ic50)

docetaxel_test_scaled_stomach <- docetaxel_rna_seq_test_scaled[docetaxel_stomach_lines, ]
docetaxel_test_stomach <- docetaxel_test[which(docetaxel_test$Tissue == 'stomach'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_stomach), s = 'lambda.min', interval = 'conf')
docetaxel_stomach_least_auc <- auc(docetaxel_test_stomach$least_sensitive, new_ic50)

docetaxel_test_scaled_thyroid <- docetaxel_rna_seq_test_scaled[docetaxel_thyroid_lines, ]
docetaxel_test_thyroid <- docetaxel_test[which(docetaxel_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_thyroid), s = 'lambda.min', interval = 'conf')
docetaxel_thyroid_least_auc <- auc(docetaxel_test_thyroid$least_sensitive, new_ic50)

docetaxel_test_scaled_upper_aerodigestive_tract <- docetaxel_rna_seq_test_scaled[docetaxel_upper_aerodigestive_tract_lines, ]
docetaxel_test_upper_aerodigestive_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_upper_aerodigestive_tract), s = 'lambda.min', interval = 'conf')
docetaxel_upper_aerodigestive_tract_least_auc <- auc(docetaxel_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

docetaxel_test_scaled_urinary_tract <- docetaxel_rna_seq_test_scaled[docetaxel_urinary_tract_lines, ]
docetaxel_test_urinary_tract <- docetaxel_test[which(docetaxel_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_test_scaled_urinary_tract), s = 'lambda.min', interval = 'conf')
docetaxel_urinary_tract_least_auc <- auc(docetaxel_test_urinary_tract$least_sensitive, new_ic50)

docetaxel_least_auc <- c(docetaxel_autonomic_ganglia_least_auc, docetaxel_biliary_tract_least_auc, 
                           docetaxel_bone_least_auc, docetaxel_breast_least_auc, docetaxel_central_nervous_system_least_auc, 
                           docetaxel_endometrium_least_auc, docetaxel_kidney_least_auc, 
                           docetaxel_large_intestine_least_auc, docetaxel_liver_least_auc, 
                           docetaxel_lung_least_auc, docetaxel_oesaphagus_least_auc, 
                           docetaxel_ovary_least_auc, docetaxel_pancreas_least_auc, 
                           docetaxel_pleura_least_auc, docetaxel_prostate_least_auc, 
                           docetaxel_salivary_gland_least_auc, docetaxel_skin_least_auc, 
                           docetaxel_soft_tissue_least_auc, docetaxel_stomach_least_auc, 
                           docetaxel_thyroid_least_auc, docetaxel_upper_aerodigestive_tract_least_auc, 
                           docetaxel_urinary_tract_least_auc)

fluorouracil_autonomic_ganglia_lines              <- fluorouracil_test$Tissue == 'autonomic_ganglia'
fluorouracil_biliary_tract_lines                  <- fluorouracil_test$Tissue == 'biliary_tract'
fluorouracil_bone_lines                           <- fluorouracil_test$Tissue == 'bone'
fluorouracil_breast_lines                         <- fluorouracil_test$Tissue == 'breast'
fluorouracil_central_nervous_system_lines         <- fluorouracil_test$Tissue == 'central_nervous_system'
fluorouracil_endometrium_lines                    <- fluorouracil_test$Tissue == 'endometrium'
fluorouracil_kidney_lines                         <- fluorouracil_test$Tissue == 'kidney'
fluorouracil_large_intestine_lines                <- fluorouracil_test$Tissue == 'large_intestine'
fluorouracil_liver_lines                          <- fluorouracil_test$Tissue == 'liver'
fluorouracil_lung_lines                           <- fluorouracil_test$Tissue == 'lung'
fluorouracil_oesaphagus_lines                     <- fluorouracil_test$Tissue == 'oesphagus'
fluorouracil_ovary_lines                          <- fluorouracil_test$Tissue == 'ovary'
fluorouracil_pancreas_lines                       <- fluorouracil_test$Tissue == 'pancreas'
fluorouracil_pleura_lines                         <- fluorouracil_test$Tissue == 'pleura'
fluorouracil_prostate_lines                       <- fluorouracil_test$Tissue == 'prostate'
fluorouracil_salivary_gland_lines                 <- fluorouracil_test$Tissue == 'salivary_gland'
fluorouracil_skin_lines                           <- fluorouracil_test$Tissue == 'skin'
fluorouracil_soft_tissue_lines                    <- fluorouracil_test$Tissue == 'soft_tissue'
fluorouracil_stomach_lines                        <- fluorouracil_test$Tissue == 'stomach'
fluorouracil_thyroid_lines                        <- fluorouracil_test$Tissue == 'thyroid'
fluorouracil_upper_aerodigestive_tract_lines      <- fluorouracil_test$Tissue == 'upper_aerodigestive_tract'
fluorouracil_urinary_tract_lines                  <- fluorouracil_test$Tissue == 'urinary_tract'

fluorouracil_test_scaled_autonomic_ganglia <- fluorouracil_rna_seq_test_scaled[fluorouracil_autonomic_ganglia_lines, ]
fluorouracil_test_autonomic_ganglia <- fluorouracil_test[which(fluorouracil_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
fluorouracil_autonomic_ganglia_most_auc <- auc(fluorouracil_test_autonomic_ganglia$most_sensitive, new_ic50)

fluorouracil_test_scaled_biliary_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_biliary_tract_lines, ]
fluorouracil_test_biliary_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_biliary_tract_most_auc <- auc(fluorouracil_test_biliary_tract$most_sensitive, new_ic50)

fluorouracil_test_scaled_bone <- fluorouracil_rna_seq_test_scaled[fluorouracil_bone_lines, ]
fluorouracil_test_bone <- fluorouracil_test[which(fluorouracil_test$Tissue == 'bone'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
fluorouracil_bone_most_auc <- auc(fluorouracil_test_bone$most_sensitive, new_ic50)

fluorouracil_test_scaled_breast <- fluorouracil_rna_seq_test_scaled[fluorouracil_breast_lines, ]
fluorouracil_test_breast <- fluorouracil_test[which(fluorouracil_test$Tissue == 'breast'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
fluorouracil_breast_most_auc <- auc(fluorouracil_test_breast$most_sensitive, new_ic50)

fluorouracil_test_scaled_central_nervous_system <- fluorouracil_rna_seq_test_scaled[fluorouracil_central_nervous_system_lines, ]
fluorouracil_test_central_nervous_system <- fluorouracil_test[which(fluorouracil_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
fluorouracil_central_nervous_system_most_auc <- auc(fluorouracil_test_central_nervous_system$most_sensitive, new_ic50)

fluorouracil_test_scaled_endometrium <- fluorouracil_rna_seq_test_scaled[fluorouracil_endometrium_lines, ]
fluorouracil_test_endometrium <- fluorouracil_test[which(fluorouracil_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
fluorouracil_endometrium_most_auc <- auc(fluorouracil_test_endometrium$most_sensitive, new_ic50)

fluorouracil_test_scaled_kidney <- fluorouracil_rna_seq_test_scaled[fluorouracil_kidney_lines, ]
fluorouracil_test_kidney <- fluorouracil_test[which(fluorouracil_test$Tissue == 'kidney'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
fluorouracil_kidney_most_auc <- auc(fluorouracil_test_kidney$most_sensitive, new_ic50)

fluorouracil_test_scaled_large_intestine <- fluorouracil_rna_seq_test_scaled[fluorouracil_large_intestine_lines, ]
fluorouracil_test_large_intestine <- fluorouracil_test[which(fluorouracil_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
fluorouracil_large_intestine_most_auc <- auc(fluorouracil_test_large_intestine$most_sensitive, new_ic50)

fluorouracil_test_scaled_liver <- fluorouracil_rna_seq_test_scaled[fluorouracil_liver_lines, ]
fluorouracil_test_liver <- fluorouracil_test[which(fluorouracil_test$Tissue == 'liver'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
fluorouracil_liver_most_auc <- auc(fluorouracil_test_liver$most_sensitive, new_ic50)

fluorouracil_test_scaled_lung <- fluorouracil_rna_seq_test_scaled[fluorouracil_lung_lines, ]
fluorouracil_test_lung <- fluorouracil_test[which(fluorouracil_test$Tissue == 'lung'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
fluorouracil_lung_most_auc <- auc(fluorouracil_test_lung$most_sensitive, new_ic50)

fluorouracil_test_scaled_oesaphagus <- fluorouracil_rna_seq_test_scaled[fluorouracil_oesaphagus_lines, ]
fluorouracil_test_oesaphagus <- fluorouracil_test[which(fluorouracil_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
fluorouracil_oesaphagus_most_auc <- auc(fluorouracil_test_oesaphagus$most_sensitive, new_ic50)

fluorouracil_test_scaled_ovary <- fluorouracil_rna_seq_test_scaled[fluorouracil_ovary_lines, ]
fluorouracil_test_ovary <- fluorouracil_test[which(fluorouracil_test$Tissue == 'ovary'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
fluorouracil_ovary_most_auc <- auc(fluorouracil_test_ovary$most_sensitive, new_ic50)

fluorouracil_test_scaled_pancreas <- fluorouracil_rna_seq_test_scaled[fluorouracil_pancreas_lines, ]
fluorouracil_test_pancreas <- fluorouracil_test[which(fluorouracil_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
fluorouracil_pancreas_most_auc <- auc(fluorouracil_test_pancreas$most_sensitive, new_ic50)

fluorouracil_test_scaled_pleura <- fluorouracil_rna_seq_test_scaled[fluorouracil_pleura_lines, ]
fluorouracil_test_pleura <- fluorouracil_test[which(fluorouracil_test$Tissue == 'pleura'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
fluorouracil_pleura_most_auc <- auc(fluorouracil_test_pleura$most_sensitive, new_ic50)

fluorouracil_test_scaled_prostate <- fluorouracil_rna_seq_test_scaled[fluorouracil_prostate_lines, ]
fluorouracil_test_prostate <- fluorouracil_test[which(fluorouracil_test$Tissue == 'prostate'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
fluorouracil_prostate_most_auc <- auc(fluorouracil_test_prostate$most_sensitive, new_ic50)

fluorouracil_test_scaled_salivary_gland <- fluorouracil_rna_seq_test_scaled[fluorouracil_salivary_gland_lines, ]
fluorouracil_test_salivary_gland <- fluorouracil_test[which(fluorouracil_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
fluorouracil_salivary_gland_most_auc <- auc(fluorouracil_test_salivary_gland$most_sensitive, new_ic50)

fluorouracil_test_scaled_skin <- fluorouracil_rna_seq_test_scaled[fluorouracil_skin_lines, ]
fluorouracil_test_skin <- fluorouracil_test[which(fluorouracil_test$Tissue == 'skin'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
fluorouracil_skin_most_auc <- auc(fluorouracil_test_skin$most_sensitive, new_ic50)

fluorouracil_test_scaled_soft_tissue <- fluorouracil_rna_seq_test_scaled[fluorouracil_soft_tissue_lines, ]
fluorouracil_test_soft_tissue <- fluorouracil_test[which(fluorouracil_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
fluorouracil_soft_tissue_most_auc <- auc(fluorouracil_test_soft_tissue$most_sensitive, new_ic50)

fluorouracil_test_scaled_stomach <- fluorouracil_rna_seq_test_scaled[fluorouracil_stomach_lines, ]
fluorouracil_test_stomach <- fluorouracil_test[which(fluorouracil_test$Tissue == 'stomach'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
fluorouracil_stomach_most_auc <- auc(fluorouracil_test_stomach$most_sensitive, new_ic50)

fluorouracil_test_scaled_thyroid <- fluorouracil_rna_seq_test_scaled[fluorouracil_thyroid_lines, ]
fluorouracil_test_thyroid <- fluorouracil_test[which(fluorouracil_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
fluorouracil_thyroid_most_auc <- auc(fluorouracil_test_thyroid$most_sensitive, new_ic50)

fluorouracil_test_scaled_upper_aerodigestive_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_upper_aerodigestive_tract_lines, ]
fluorouracil_test_upper_aerodigestive_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_upper_aerodigestive_tract_most_auc <- auc(fluorouracil_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

fluorouracil_test_scaled_urinary_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_urinary_tract_lines, ]
fluorouracil_test_urinary_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_urinary_tract_most_auc <- auc(fluorouracil_test_urinary_tract$most_sensitive, new_ic50)

fluorouracil_most_auc <- c(fluorouracil_autonomic_ganglia_most_auc, fluorouracil_biliary_tract_most_auc, 
                          fluorouracil_bone_most_auc, fluorouracil_breast_most_auc, fluorouracil_central_nervous_system_most_auc, 
                          fluorouracil_endometrium_most_auc, fluorouracil_kidney_most_auc, 
                          fluorouracil_large_intestine_most_auc, fluorouracil_liver_most_auc, 
                          fluorouracil_lung_most_auc, fluorouracil_oesaphagus_most_auc, 
                          fluorouracil_ovary_most_auc, fluorouracil_pancreas_most_auc, 
                          fluorouracil_pleura_most_auc, fluorouracil_prostate_most_auc, 
                          fluorouracil_salivary_gland_most_auc, fluorouracil_skin_most_auc, 
                          fluorouracil_soft_tissue_most_auc, fluorouracil_stomach_most_auc, 
                          fluorouracil_thyroid_most_auc, fluorouracil_upper_aerodigestive_tract_most_auc, 
                          fluorouracil_urinary_tract_most_auc)

fluorouracil_test_scaled_autonomic_ganglia <- fluorouracil_rna_seq_test_scaled[fluorouracil_autonomic_ganglia_lines, ]
fluorouracil_test_autonomic_ganglia <- fluorouracil_test[which(fluorouracil_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
fluorouracil_autonomic_ganglia_least_auc <- auc(fluorouracil_test_autonomic_ganglia$least_sensitive, new_ic50)

fluorouracil_test_scaled_biliary_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_biliary_tract_lines, ]
fluorouracil_test_biliary_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_biliary_tract_least_auc <- auc(fluorouracil_test_biliary_tract$least_sensitive, new_ic50)

fluorouracil_test_scaled_bone <- fluorouracil_rna_seq_test_scaled[fluorouracil_bone_lines, ]
fluorouracil_test_bone <- fluorouracil_test[which(fluorouracil_test$Tissue == 'bone'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
fluorouracil_bone_least_auc <- auc(fluorouracil_test_bone$least_sensitive, new_ic50)

fluorouracil_test_scaled_breast <- fluorouracil_rna_seq_test_scaled[fluorouracil_breast_lines, ]
fluorouracil_test_breast <- fluorouracil_test[which(fluorouracil_test$Tissue == 'breast'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
fluorouracil_breast_least_auc <- auc(fluorouracil_test_breast$least_sensitive, new_ic50)

fluorouracil_test_scaled_central_nervous_system <- fluorouracil_rna_seq_test_scaled[fluorouracil_central_nervous_system_lines, ]
fluorouracil_test_central_nervous_system <- fluorouracil_test[which(fluorouracil_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
fluorouracil_central_nervous_system_least_auc <- auc(fluorouracil_test_central_nervous_system$least_sensitive, new_ic50)

fluorouracil_test_scaled_endometrium <- fluorouracil_rna_seq_test_scaled[fluorouracil_endometrium_lines, ]
fluorouracil_test_endometrium <- fluorouracil_test[which(fluorouracil_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
fluorouracil_endometrium_least_auc <- auc(fluorouracil_test_endometrium$least_sensitive, new_ic50)

fluorouracil_test_scaled_kidney <- fluorouracil_rna_seq_test_scaled[fluorouracil_kidney_lines, ]
fluorouracil_test_kidney <- fluorouracil_test[which(fluorouracil_test$Tissue == 'kidney'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
fluorouracil_kidney_least_auc <- auc(fluorouracil_test_kidney$least_sensitive, new_ic50)

fluorouracil_test_scaled_large_intestine <- fluorouracil_rna_seq_test_scaled[fluorouracil_large_intestine_lines, ]
fluorouracil_test_large_intestine <- fluorouracil_test[which(fluorouracil_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
fluorouracil_large_intestine_least_auc <- auc(fluorouracil_test_large_intestine$least_sensitive, new_ic50)

fluorouracil_test_scaled_liver <- fluorouracil_rna_seq_test_scaled[fluorouracil_liver_lines, ]
fluorouracil_test_liver <- fluorouracil_test[which(fluorouracil_test$Tissue == 'liver'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
fluorouracil_liver_least_auc <- auc(fluorouracil_test_liver$least_sensitive, new_ic50)

fluorouracil_test_scaled_lung <- fluorouracil_rna_seq_test_scaled[fluorouracil_lung_lines, ]
fluorouracil_test_lung <- fluorouracil_test[which(fluorouracil_test$Tissue == 'lung'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
fluorouracil_lung_least_auc <- auc(fluorouracil_test_lung$least_sensitive, new_ic50)

fluorouracil_test_scaled_oesaphagus <- fluorouracil_rna_seq_test_scaled[fluorouracil_oesaphagus_lines, ]
fluorouracil_test_oesaphagus <- fluorouracil_test[which(fluorouracil_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
fluorouracil_oesaphagus_least_auc <- auc(fluorouracil_test_oesaphagus$least_sensitive, new_ic50)

fluorouracil_test_scaled_ovary <- fluorouracil_rna_seq_test_scaled[fluorouracil_ovary_lines, ]
fluorouracil_test_ovary <- fluorouracil_test[which(fluorouracil_test$Tissue == 'ovary'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
fluorouracil_ovary_least_auc <- auc(fluorouracil_test_ovary$least_sensitive, new_ic50)

fluorouracil_test_scaled_pancreas <- fluorouracil_rna_seq_test_scaled[fluorouracil_pancreas_lines, ]
fluorouracil_test_pancreas <- fluorouracil_test[which(fluorouracil_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
fluorouracil_pancreas_least_auc <- auc(fluorouracil_test_pancreas$least_sensitive, new_ic50)

fluorouracil_test_scaled_pleura <- fluorouracil_rna_seq_test_scaled[fluorouracil_pleura_lines, ]
fluorouracil_test_pleura <- fluorouracil_test[which(fluorouracil_test$Tissue == 'pleura'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
fluorouracil_pleura_least_auc <- auc(fluorouracil_test_pleura$least_sensitive, new_ic50)

fluorouracil_test_scaled_prostate <- fluorouracil_rna_seq_test_scaled[fluorouracil_prostate_lines, ]
fluorouracil_test_prostate <- fluorouracil_test[which(fluorouracil_test$Tissue == 'prostate'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
fluorouracil_prostate_least_auc <- auc(fluorouracil_test_prostate$least_sensitive, new_ic50)

fluorouracil_test_scaled_salivary_gland <- fluorouracil_rna_seq_test_scaled[fluorouracil_salivary_gland_lines, ]
fluorouracil_test_salivary_gland <- fluorouracil_test[which(fluorouracil_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
fluorouracil_salivary_gland_least_auc <- auc(fluorouracil_test_salivary_gland$least_sensitive, new_ic50)

fluorouracil_test_scaled_skin <- fluorouracil_rna_seq_test_scaled[fluorouracil_skin_lines, ]
fluorouracil_test_skin <- fluorouracil_test[which(fluorouracil_test$Tissue == 'skin'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
fluorouracil_skin_least_auc <- auc(fluorouracil_test_skin$least_sensitive, new_ic50)

fluorouracil_test_scaled_soft_tissue <- fluorouracil_rna_seq_test_scaled[fluorouracil_soft_tissue_lines, ]
fluorouracil_test_soft_tissue <- fluorouracil_test[which(fluorouracil_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
fluorouracil_soft_tissue_least_auc <- auc(fluorouracil_test_soft_tissue$least_sensitive, new_ic50)

fluorouracil_test_scaled_stomach <- fluorouracil_rna_seq_test_scaled[fluorouracil_stomach_lines, ]
fluorouracil_test_stomach <- fluorouracil_test[which(fluorouracil_test$Tissue == 'stomach'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
fluorouracil_stomach_least_auc <- auc(fluorouracil_test_stomach$least_sensitive, new_ic50)

fluorouracil_test_scaled_thyroid <- fluorouracil_rna_seq_test_scaled[fluorouracil_thyroid_lines, ]
fluorouracil_test_thyroid <- fluorouracil_test[which(fluorouracil_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
fluorouracil_thyroid_least_auc <- auc(fluorouracil_test_thyroid$least_sensitive, new_ic50)

fluorouracil_test_scaled_upper_aerodigestive_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_upper_aerodigestive_tract_lines, ]
fluorouracil_test_upper_aerodigestive_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_upper_aerodigestive_tract_least_auc <- auc(fluorouracil_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

fluorouracil_test_scaled_urinary_tract <- fluorouracil_rna_seq_test_scaled[fluorouracil_urinary_tract_lines, ]
fluorouracil_test_urinary_tract <- fluorouracil_test[which(fluorouracil_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
fluorouracil_urinary_tract_least_auc <- auc(fluorouracil_test_urinary_tract$least_sensitive, new_ic50)

fluorouracil_least_auc <- c(fluorouracil_autonomic_ganglia_least_auc, fluorouracil_biliary_tract_least_auc, 
                           fluorouracil_bone_least_auc, fluorouracil_breast_least_auc, fluorouracil_central_nervous_system_least_auc, 
                           fluorouracil_endometrium_least_auc, fluorouracil_kidney_least_auc, 
                           fluorouracil_large_intestine_least_auc, fluorouracil_liver_least_auc, 
                           fluorouracil_lung_least_auc, fluorouracil_oesaphagus_least_auc, 
                           fluorouracil_ovary_least_auc, fluorouracil_pancreas_least_auc, 
                           fluorouracil_pleura_least_auc, fluorouracil_prostate_least_auc, 
                           fluorouracil_salivary_gland_least_auc, fluorouracil_skin_least_auc, 
                           fluorouracil_soft_tissue_least_auc, fluorouracil_stomach_least_auc, 
                           fluorouracil_thyroid_least_auc, fluorouracil_upper_aerodigestive_tract_least_auc, 
                           fluorouracil_urinary_tract_least_auc)

gemcitabine_autonomic_ganglia_lines              <- gemcitabine_test$Tissue == 'autonomic_ganglia'
gemcitabine_biliary_tract_lines                  <- gemcitabine_test$Tissue == 'biliary_tract'
gemcitabine_bone_lines                           <- gemcitabine_test$Tissue == 'bone'
gemcitabine_breast_lines                         <- gemcitabine_test$Tissue == 'breast'
gemcitabine_central_nervous_system_lines         <- gemcitabine_test$Tissue == 'central_nervous_system'
gemcitabine_endometrium_lines                    <- gemcitabine_test$Tissue == 'endometrium'
gemcitabine_kidney_lines                         <- gemcitabine_test$Tissue == 'kidney'
gemcitabine_large_intestine_lines                <- gemcitabine_test$Tissue == 'large_intestine'
gemcitabine_liver_lines                          <- gemcitabine_test$Tissue == 'liver'
gemcitabine_lung_lines                           <- gemcitabine_test$Tissue == 'lung'
gemcitabine_oesaphagus_lines                     <- gemcitabine_test$Tissue == 'oesphagus'
gemcitabine_ovary_lines                          <- gemcitabine_test$Tissue == 'ovary'
gemcitabine_pancreas_lines                       <- gemcitabine_test$Tissue == 'pancreas'
gemcitabine_pleura_lines                         <- gemcitabine_test$Tissue == 'pleura'
gemcitabine_prostate_lines                       <- gemcitabine_test$Tissue == 'prostate'
gemcitabine_salivary_gland_lines                 <- gemcitabine_test$Tissue == 'salivary_gland'
gemcitabine_skin_lines                           <- gemcitabine_test$Tissue == 'skin'
gemcitabine_soft_tissue_lines                    <- gemcitabine_test$Tissue == 'soft_tissue'
gemcitabine_stomach_lines                        <- gemcitabine_test$Tissue == 'stomach'
gemcitabine_thyroid_lines                        <- gemcitabine_test$Tissue == 'thyroid'
gemcitabine_upper_aerodigestive_tract_lines      <- gemcitabine_test$Tissue == 'upper_aerodigestive_tract'
gemcitabine_urinary_tract_lines                  <- gemcitabine_test$Tissue == 'urinary_tract'

gemcitabine_test_scaled_autonomic_ganglia <- gemcitabine_rna_seq_test_scaled[gemcitabine_autonomic_ganglia_lines, ]
gemcitabine_test_autonomic_ganglia <- gemcitabine_test[which(gemcitabine_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
gemcitabine_autonomic_ganglia_most_auc <- auc(gemcitabine_test_autonomic_ganglia$most_sensitive, new_ic50)

gemcitabine_test_scaled_biliary_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_biliary_tract_lines, ]
gemcitabine_test_biliary_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_biliary_tract_most_auc <- auc(gemcitabine_test_biliary_tract$most_sensitive, new_ic50)

gemcitabine_test_scaled_bone <- gemcitabine_rna_seq_test_scaled[gemcitabine_bone_lines, ]
gemcitabine_test_bone <- gemcitabine_test[which(gemcitabine_test$Tissue == 'bone'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
gemcitabine_bone_most_auc <- auc(gemcitabine_test_bone$most_sensitive, new_ic50)

gemcitabine_test_scaled_breast <- gemcitabine_rna_seq_test_scaled[gemcitabine_breast_lines, ]
gemcitabine_test_breast <- gemcitabine_test[which(gemcitabine_test$Tissue == 'breast'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
gemcitabine_breast_most_auc <- auc(gemcitabine_test_breast$most_sensitive, new_ic50)

gemcitabine_test_scaled_central_nervous_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_central_nervous_system_lines, ]
gemcitabine_test_central_nervous_system <- gemcitabine_test[which(gemcitabine_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_central_nervous_system_most_auc <- auc(gemcitabine_test_central_nervous_system$most_sensitive, new_ic50)

gemcitabine_test_scaled_endometrium <- gemcitabine_rna_seq_test_scaled[gemcitabine_endometrium_lines, ]
gemcitabine_test_endometrium <- gemcitabine_test[which(gemcitabine_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
gemcitabine_endometrium_most_auc <- auc(gemcitabine_test_endometrium$most_sensitive, new_ic50)

gemcitabine_test_scaled_kidney <- gemcitabine_rna_seq_test_scaled[gemcitabine_kidney_lines, ]
gemcitabine_test_kidney <- gemcitabine_test[which(gemcitabine_test$Tissue == 'kidney'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
gemcitabine_kidney_most_auc <- auc(gemcitabine_test_kidney$most_sensitive, new_ic50)

gemcitabine_test_scaled_large_intestine <- gemcitabine_rna_seq_test_scaled[gemcitabine_large_intestine_lines, ]
gemcitabine_test_large_intestine <- gemcitabine_test[which(gemcitabine_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
gemcitabine_large_intestine_most_auc <- auc(gemcitabine_test_large_intestine$most_sensitive, new_ic50)

gemcitabine_test_scaled_liver <- gemcitabine_rna_seq_test_scaled[gemcitabine_liver_lines, ]
gemcitabine_test_liver <- gemcitabine_test[which(gemcitabine_test$Tissue == 'liver'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
gemcitabine_liver_most_auc <- auc(gemcitabine_test_liver$most_sensitive, new_ic50)

gemcitabine_test_scaled_lung <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_lines, ]
gemcitabine_test_lung <- gemcitabine_test[which(gemcitabine_test$Tissue == 'lung'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_most_auc <- auc(gemcitabine_test_lung$most_sensitive, new_ic50)

gemcitabine_test_scaled_oesaphagus <- gemcitabine_rna_seq_test_scaled[gemcitabine_oesaphagus_lines, ]
gemcitabine_test_oesaphagus <- gemcitabine_test[which(gemcitabine_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
gemcitabine_oesaphagus_most_auc <- auc(gemcitabine_test_oesaphagus$most_sensitive, new_ic50)

gemcitabine_test_scaled_ovary <- gemcitabine_rna_seq_test_scaled[gemcitabine_ovary_lines, ]
gemcitabine_test_ovary <- gemcitabine_test[which(gemcitabine_test$Tissue == 'ovary'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
gemcitabine_ovary_most_auc <- auc(gemcitabine_test_ovary$most_sensitive, new_ic50)

gemcitabine_test_scaled_pancreas <- gemcitabine_rna_seq_test_scaled[gemcitabine_pancreas_lines, ]
gemcitabine_test_pancreas <- gemcitabine_test[which(gemcitabine_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
gemcitabine_pancreas_most_auc <- auc(gemcitabine_test_pancreas$most_sensitive, new_ic50)

gemcitabine_test_scaled_pleura <- gemcitabine_rna_seq_test_scaled[gemcitabine_pleura_lines, ]
gemcitabine_test_pleura <- gemcitabine_test[which(gemcitabine_test$Tissue == 'pleura'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
gemcitabine_pleura_most_auc <- auc(gemcitabine_test_pleura$most_sensitive, new_ic50)

gemcitabine_test_scaled_prostate <- gemcitabine_rna_seq_test_scaled[gemcitabine_prostate_lines, ]
gemcitabine_test_prostate <- gemcitabine_test[which(gemcitabine_test$Tissue == 'prostate'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
gemcitabine_prostate_most_auc <- auc(gemcitabine_test_prostate$most_sensitive, new_ic50)

gemcitabine_test_scaled_salivary_gland <- gemcitabine_rna_seq_test_scaled[gemcitabine_salivary_gland_lines, ]
gemcitabine_test_salivary_gland <- gemcitabine_test[which(gemcitabine_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
gemcitabine_salivary_gland_most_auc <- auc(gemcitabine_test_salivary_gland$most_sensitive, new_ic50)

gemcitabine_test_scaled_skin <- gemcitabine_rna_seq_test_scaled[gemcitabine_skin_lines, ]
gemcitabine_test_skin <- gemcitabine_test[which(gemcitabine_test$Tissue == 'skin'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
gemcitabine_skin_most_auc <- auc(gemcitabine_test_skin$most_sensitive, new_ic50)

gemcitabine_test_scaled_soft_tissue <- gemcitabine_rna_seq_test_scaled[gemcitabine_soft_tissue_lines, ]
gemcitabine_test_soft_tissue <- gemcitabine_test[which(gemcitabine_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
gemcitabine_soft_tissue_most_auc <- auc(gemcitabine_test_soft_tissue$most_sensitive, new_ic50)

gemcitabine_test_scaled_stomach <- gemcitabine_rna_seq_test_scaled[gemcitabine_stomach_lines, ]
gemcitabine_test_stomach <- gemcitabine_test[which(gemcitabine_test$Tissue == 'stomach'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
gemcitabine_stomach_most_auc <- auc(gemcitabine_test_stomach$most_sensitive, new_ic50)

gemcitabine_test_scaled_thyroid <- gemcitabine_rna_seq_test_scaled[gemcitabine_thyroid_lines, ]
gemcitabine_test_thyroid <- gemcitabine_test[which(gemcitabine_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
gemcitabine_thyroid_most_auc <- auc(gemcitabine_test_thyroid$most_sensitive, new_ic50)

gemcitabine_test_scaled_upper_aerodigestive_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_upper_aerodigestive_tract_lines, ]
gemcitabine_test_upper_aerodigestive_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_upper_aerodigestive_tract_most_auc <- auc(gemcitabine_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

gemcitabine_test_scaled_urinary_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_urinary_tract_lines, ]
gemcitabine_test_urinary_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_urinary_tract_most_auc <- auc(gemcitabine_test_urinary_tract$most_sensitive, new_ic50)

gemcitabine_most_auc <- c(gemcitabine_autonomic_ganglia_most_auc, gemcitabine_biliary_tract_most_auc, 
                          gemcitabine_bone_most_auc, gemcitabine_breast_most_auc, gemcitabine_central_nervous_system_most_auc, 
                          gemcitabine_endometrium_most_auc, gemcitabine_kidney_most_auc, 
                          gemcitabine_large_intestine_most_auc, gemcitabine_liver_most_auc, 
                          gemcitabine_lung_most_auc, gemcitabine_oesaphagus_most_auc, 
                          gemcitabine_ovary_most_auc, gemcitabine_pancreas_most_auc, 
                          gemcitabine_pleura_most_auc, gemcitabine_prostate_most_auc, 
                          gemcitabine_salivary_gland_most_auc, gemcitabine_skin_most_auc, 
                          gemcitabine_soft_tissue_most_auc, gemcitabine_stomach_most_auc, 
                          gemcitabine_thyroid_most_auc, gemcitabine_upper_aerodigestive_tract_most_auc, 
                          gemcitabine_urinary_tract_most_auc)

gemcitabine_test_scaled_autonomic_ganglia <- gemcitabine_rna_seq_test_scaled[gemcitabine_autonomic_ganglia_lines, ]
gemcitabine_test_autonomic_ganglia <- gemcitabine_test[which(gemcitabine_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
gemcitabine_autonomic_ganglia_least_auc <- auc(gemcitabine_test_autonomic_ganglia$least_sensitive, new_ic50)

gemcitabine_test_scaled_biliary_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_biliary_tract_lines, ]
gemcitabine_test_biliary_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_biliary_tract_least_auc <- auc(gemcitabine_test_biliary_tract$least_sensitive, new_ic50)

gemcitabine_test_scaled_bone <- gemcitabine_rna_seq_test_scaled[gemcitabine_bone_lines, ]
gemcitabine_test_bone <- gemcitabine_test[which(gemcitabine_test$Tissue == 'bone'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
gemcitabine_bone_least_auc <- auc(gemcitabine_test_bone$least_sensitive, new_ic50)

gemcitabine_test_scaled_breast <- gemcitabine_rna_seq_test_scaled[gemcitabine_breast_lines, ]
gemcitabine_test_breast <- gemcitabine_test[which(gemcitabine_test$Tissue == 'breast'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
gemcitabine_breast_least_auc <- auc(gemcitabine_test_breast$least_sensitive, new_ic50)

gemcitabine_test_scaled_central_nervous_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_central_nervous_system_lines, ]
gemcitabine_test_central_nervous_system <- gemcitabine_test[which(gemcitabine_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_central_nervous_system_least_auc <- auc(gemcitabine_test_central_nervous_system$least_sensitive, new_ic50)

gemcitabine_test_scaled_endometrium <- gemcitabine_rna_seq_test_scaled[gemcitabine_endometrium_lines, ]
gemcitabine_test_endometrium <- gemcitabine_test[which(gemcitabine_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
gemcitabine_endometrium_least_auc <- auc(gemcitabine_test_endometrium$least_sensitive, new_ic50)

gemcitabine_test_scaled_kidney <- gemcitabine_rna_seq_test_scaled[gemcitabine_kidney_lines, ]
gemcitabine_test_kidney <- gemcitabine_test[which(gemcitabine_test$Tissue == 'kidney'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
gemcitabine_kidney_least_auc <- auc(gemcitabine_test_kidney$least_sensitive, new_ic50)

gemcitabine_test_scaled_large_intestine <- gemcitabine_rna_seq_test_scaled[gemcitabine_large_intestine_lines, ]
gemcitabine_test_large_intestine <- gemcitabine_test[which(gemcitabine_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
gemcitabine_large_intestine_least_auc <- auc(gemcitabine_test_large_intestine$least_sensitive, new_ic50)

gemcitabine_test_scaled_liver <- gemcitabine_rna_seq_test_scaled[gemcitabine_liver_lines, ]
gemcitabine_test_liver <- gemcitabine_test[which(gemcitabine_test$Tissue == 'liver'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
gemcitabine_liver_least_auc <- auc(gemcitabine_test_liver$least_sensitive, new_ic50)

gemcitabine_test_scaled_lung <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_lines, ]
gemcitabine_test_lung <- gemcitabine_test[which(gemcitabine_test$Tissue == 'lung'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_least_auc <- auc(gemcitabine_test_lung$least_sensitive, new_ic50)

gemcitabine_test_scaled_oesaphagus <- gemcitabine_rna_seq_test_scaled[gemcitabine_oesaphagus_lines, ]
gemcitabine_test_oesaphagus <- gemcitabine_test[which(gemcitabine_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
gemcitabine_oesaphagus_least_auc <- auc(gemcitabine_test_oesaphagus$least_sensitive, new_ic50)

gemcitabine_test_scaled_ovary <- gemcitabine_rna_seq_test_scaled[gemcitabine_ovary_lines, ]
gemcitabine_test_ovary <- gemcitabine_test[which(gemcitabine_test$Tissue == 'ovary'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
gemcitabine_ovary_least_auc <- auc(gemcitabine_test_ovary$least_sensitive, new_ic50)

gemcitabine_test_scaled_pancreas <- gemcitabine_rna_seq_test_scaled[gemcitabine_pancreas_lines, ]
gemcitabine_test_pancreas <- gemcitabine_test[which(gemcitabine_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
gemcitabine_pancreas_least_auc <- auc(gemcitabine_test_pancreas$least_sensitive, new_ic50)

gemcitabine_test_scaled_pleura <- gemcitabine_rna_seq_test_scaled[gemcitabine_pleura_lines, ]
gemcitabine_test_pleura <- gemcitabine_test[which(gemcitabine_test$Tissue == 'pleura'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
gemcitabine_pleura_least_auc <- auc(gemcitabine_test_pleura$least_sensitive, new_ic50)

gemcitabine_test_scaled_prostate <- gemcitabine_rna_seq_test_scaled[gemcitabine_prostate_lines, ]
gemcitabine_test_prostate <- gemcitabine_test[which(gemcitabine_test$Tissue == 'prostate'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
gemcitabine_prostate_least_auc <- auc(gemcitabine_test_prostate$least_sensitive, new_ic50)

gemcitabine_test_scaled_salivary_gland <- gemcitabine_rna_seq_test_scaled[gemcitabine_salivary_gland_lines, ]
gemcitabine_test_salivary_gland <- gemcitabine_test[which(gemcitabine_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
gemcitabine_salivary_gland_least_auc <- auc(gemcitabine_test_salivary_gland$least_sensitive, new_ic50)

gemcitabine_test_scaled_skin <- gemcitabine_rna_seq_test_scaled[gemcitabine_skin_lines, ]
gemcitabine_test_skin <- gemcitabine_test[which(gemcitabine_test$Tissue == 'skin'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
gemcitabine_skin_least_auc <- auc(gemcitabine_test_skin$least_sensitive, new_ic50)

gemcitabine_test_scaled_soft_tissue <- gemcitabine_rna_seq_test_scaled[gemcitabine_soft_tissue_lines, ]
gemcitabine_test_soft_tissue <- gemcitabine_test[which(gemcitabine_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
gemcitabine_soft_tissue_least_auc <- auc(gemcitabine_test_soft_tissue$least_sensitive, new_ic50)

gemcitabine_test_scaled_stomach <- gemcitabine_rna_seq_test_scaled[gemcitabine_stomach_lines, ]
gemcitabine_test_stomach <- gemcitabine_test[which(gemcitabine_test$Tissue == 'stomach'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
gemcitabine_stomach_least_auc <- auc(gemcitabine_test_stomach$least_sensitive, new_ic50)

gemcitabine_test_scaled_thyroid <- gemcitabine_rna_seq_test_scaled[gemcitabine_thyroid_lines, ]
gemcitabine_test_thyroid <- gemcitabine_test[which(gemcitabine_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
gemcitabine_thyroid_least_auc <- auc(gemcitabine_test_thyroid$least_sensitive, new_ic50)

gemcitabine_test_scaled_upper_aerodigestive_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_upper_aerodigestive_tract_lines, ]
gemcitabine_test_upper_aerodigestive_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_upper_aerodigestive_tract_least_auc <- auc(gemcitabine_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

gemcitabine_test_scaled_urinary_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_urinary_tract_lines, ]
gemcitabine_test_urinary_tract <- gemcitabine_test[which(gemcitabine_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_urinary_tract_least_auc <- auc(gemcitabine_test_urinary_tract$least_sensitive, new_ic50)

gemcitabine_least_auc <- c(gemcitabine_autonomic_ganglia_least_auc, gemcitabine_biliary_tract_least_auc, 
                           gemcitabine_bone_least_auc, gemcitabine_breast_least_auc, gemcitabine_central_nervous_system_least_auc, 
                           gemcitabine_endometrium_least_auc, gemcitabine_kidney_least_auc, 
                           gemcitabine_large_intestine_least_auc, gemcitabine_liver_least_auc, 
                           gemcitabine_lung_least_auc, gemcitabine_oesaphagus_least_auc, 
                           gemcitabine_ovary_least_auc, gemcitabine_pancreas_least_auc, 
                           gemcitabine_pleura_least_auc, gemcitabine_prostate_least_auc, 
                           gemcitabine_salivary_gland_least_auc, gemcitabine_skin_least_auc, 
                           gemcitabine_soft_tissue_least_auc, gemcitabine_stomach_least_auc, 
                           gemcitabine_thyroid_least_auc, gemcitabine_upper_aerodigestive_tract_least_auc, 
                           gemcitabine_urinary_tract_least_auc)


all_ccle_auc_by_type <- data.frame(carboplatin_most_auc, carboplatin_least_auc, 
                                   cyclophosphamide_most_auc, cyclophosphamide_least_auc, 
                                   docetaxel_most_auc, docetaxel_least_auc, 
                                   fluorouracil_most_auc, fluorouracil_least_auc, 
                                   gemcitabine_most_auc, gemcitabine_least_auc)

all_ccle_auc_by_type <- data.frame(t(all_ccle_auc_by_type))
all_ccle_auc_by_type <- cbind(all_ccle_auc_by_type, overall_auc)
all_ccle_auc_by_type <- all_ccle_auc_by_type[, c(23, 1:22)]
rownames(all_ccle_auc_by_type) <- c('carboplatin_most', 'carboplatin_least', 
                                    'cyclophosphamide_most', 'cyclophosphamide_least', 'docetaxel_most', 
                                    'docetaxel_least', 
                                    'fluorouracil_most', 'fluorouracil_least', 
                                    'gemcitabine_most', 'gemcitabine_least')
colnames(all_ccle_auc_by_type) <- c('overall', 'autonomic_ganglia', 'biliary_tract', 'bone', 'breast', 'central_nervous_system', 'endometrium', 'kidney', 
                                    'large_intestine', 'liver', 'lung', 'oesaphagus', 'ovary', 
                                    'pancreas', 'pleura', 'prostate', 'salivary_gland', 'skin', 
                                    'soft_tissue', 'stomach', 'thyroid', 'upper_aerodigestive_tract', 'urinary_tract')

all_ccle_auc_by_type <- round(all_ccle_auc_by_type, digits = 2)
#all_ccle_auc_by_type <- data.frame(t(all_ccle_auc_by_type))
#all_ccle_auc_by_type <- all_ccle_auc_by_type[, c(23, 1:22)]
colors <- colorRampPalette(c("dodgerblue", "white", "red"))(100)

png(filename = 'Images/CCLE_AUC_heatmap.png', width = 1100)
heatmap.2(as.matrix(all_ccle_auc_by_type), trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, density.info = 'none', key.xlab = 'AUC', key.title = '', cexRow = 0.9, cexCol = 0.8, cellnote = all_ccle_auc_by_type, notecol = 'black', colsep = 1, sepwidth = c(0.1,0.1), srtCol = 45, margins = c(8,9))
dev.off()


### testing on TCGA classes treated with drug ----
# BLCA W CARBOPLATIN (8)
blca_clinical <- read.csv('Processed_Clinical_Data/blca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(blca_clinical$most_sensitive)
blca_clinical <- blca_clinical[!na_idx, ]
table(blca_clinical$drug_name)
blca_clinical_carboplatin <- blca_clinical[which(blca_clinical$drug_name == 'Carboplatin'), ]

blca_clinical_carboplatin$most_sensitive  <- ifelse(blca_clinical_carboplatin$PFS < quantile(blca_clinical_carboplatin$PFS, probs = .20), 1, 0)
blca_clinical_carboplatin$least_sensitive <- ifelse(blca_clinical_carboplatin$PFS > quantile(blca_clinical_carboplatin$PFS, probs = .80), 1, 0)

blca_gene <- read.csv('Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(blca_gene) <- gsub('\\.', '-', colnames(blca_gene))
blca_matching_idx <- blca_clinical_carboplatin$submitter_id.samples %in% colnames(blca_gene)
blca_clinical_carboplatin_short <- blca_clinical_carboplatin[blca_matching_idx, ]
blca_matching_idx <- colnames(blca_gene) %in% blca_clinical_carboplatin_short$submitter_id.samples
blca_gene_short <- blca_gene[, blca_matching_idx]
blca_gene_short <- t(blca_gene_short)
blca_gene_short_scaled <- apply(blca_gene_short, 2, scale)

rm(blca_tcga_most_min_auc)
rm(blca_tcga_most_1se_auc)
rm(blca_tcga_least_min_auc)
rm(blca_tcga_least_1se_auc)

new_blca_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_most_min_auc <- auc(blca_clinical_carboplatin_short$most_sensitive, new_blca_tcga_carboplatin_most_sensitive_min)
blca_tcga_most_min_auc <- round(blca_tcga_most_min_auc, digits = 2)

new_blca_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_most_1se_auc <- auc(blca_clinical_carboplatin_short$most_sensitive, new_blca_tcga_carboplatin_most_sensitive_1se)
blca_tcga_most_1se_auc <- round(blca_tcga_most_1se_auc, digits = 2)

new_blca_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_least_min_auc <- auc(blca_clinical_carboplatin_short$least_sensitive, new_blca_tcga_carboplatin_least_sensitive_min)
blca_tcga_least_min_auc <- round(blca_tcga_least_min_auc, digits = 2)

new_blca_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_least_1se_auc <- auc(blca_clinical_carboplatin_short$least_sensitive, new_blca_tcga_carboplatin_least_sensitive_1se)
blca_tcga_least_1se_auc <- round(blca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_blca_tcga_carboplatin_most_sensitive_min, blca_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_blca_tcga_carboplatin_most_sensitive_1se, blca_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_blca_tcga_carboplatin_least_sensitive_min, blca_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_blca_tcga_carboplatin_least_sensitive_1se, blca_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/blca_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BLCA W gemcitabine (8)
blca_clinical <- read.csv('Processed_Clinical_Data/blca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(blca_clinical$most_sensitive)
blca_clinical <- blca_clinical[!na_idx, ]
table(blca_clinical$drug_name)
blca_clinical_gemcitabine <- blca_clinical[which(blca_clinical$drug_name == 'gemcitabine' | blca_clinical$drug_name == 'Gemcitabine' | 
                                                   blca_clinical$drug_name == 'Gemzar'), ]

blca_clinical_gemcitabine$most_sensitive  <- ifelse(blca_clinical_gemcitabine$PFS < quantile(blca_clinical_gemcitabine$PFS, probs = .20), 1, 0)
blca_clinical_gemcitabine$least_sensitive <- ifelse(blca_clinical_gemcitabine$PFS > quantile(blca_clinical_gemcitabine$PFS, probs = .80), 1, 0)

blca_gene <- read.csv('Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(blca_gene) <- gsub('\\.', '-', colnames(blca_gene))
blca_matching_idx <- blca_clinical_gemcitabine$submitter_id.samples %in% colnames(blca_gene)
blca_clinical_gemcitabine_short <- blca_clinical_gemcitabine[blca_matching_idx, ]
blca_matching_idx <- colnames(blca_gene) %in% blca_clinical_gemcitabine_short$submitter_id.samples
blca_gene_short <- blca_gene[, blca_matching_idx]
blca_gene_short <- t(blca_gene_short)
blca_gene_short_scaled <- apply(blca_gene_short, 2, scale)

rm(blca_tcga_most_min_auc)
rm(blca_tcga_most_1se_auc)
rm(blca_tcga_least_min_auc)
rm(blca_tcga_least_1se_auc)

new_blca_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_most_min_auc <- auc(blca_clinical_gemcitabine_short$most_sensitive, new_blca_tcga_gemcitabine_most_sensitive_min)
blca_tcga_most_min_auc <- round(blca_tcga_most_min_auc, digits = 2)

new_blca_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_most_1se_auc <- auc(blca_clinical_gemcitabine_short$most_sensitive, new_blca_tcga_gemcitabine_most_sensitive_1se)
blca_tcga_most_1se_auc <- round(blca_tcga_most_1se_auc, digits = 2)

new_blca_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_least_min_auc <- auc(blca_clinical_gemcitabine_short$least_sensitive, new_blca_tcga_gemcitabine_least_sensitive_min)
blca_tcga_least_min_auc <- round(blca_tcga_least_min_auc, digits = 2)

new_blca_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_least_1se_auc <- auc(blca_clinical_gemcitabine_short$least_sensitive, new_blca_tcga_gemcitabine_least_sensitive_1se)
blca_tcga_least_1se_auc <- round(blca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_blca_tcga_gemcitabine_most_sensitive_min, blca_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_blca_tcga_gemcitabine_most_sensitive_1se, blca_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_blca_tcga_gemcitabine_least_sensitive_min, blca_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_blca_tcga_gemcitabine_least_sensitive_1se, blca_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/blca_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[7], lwd = 2, lty = 2, main = 'BLCA treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[7], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.62)', 'resistant model (AUC = 0.66)'), lwd = 2, lty = c(1,2), col = colors_i_need[7], cex = 0.8, bty = 'n')
dev.off()

# BRCA W CYCLOPHOSPHAMIDE (48)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_cyclophosphamide <- brca_clinical[which(brca_clinical$drug_name == 'cyclophosphamid' | brca_clinical$drug_name == 'cyclophosphamide' | 
                                                        brca_clinical$drug_name == 'Cyclophosphamide' | brca_clinical$drug_name == 'CYCLOPHOSPHAMIDE'), ]

brca_clinical_cyclophosphamide$most_sensitive  <- ifelse(brca_clinical_cyclophosphamide$PFS < quantile(brca_clinical_cyclophosphamide$PFS, probs = .20), 1, 0)
brca_clinical_cyclophosphamide$least_sensitive <- ifelse(brca_clinical_cyclophosphamide$PFS > quantile(brca_clinical_cyclophosphamide$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_cyclophosphamide$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_cyclophosphamide_short <- brca_clinical_cyclophosphamide[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_cyclophosphamide_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_cyclophosphamide_most_sensitive_min <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_cyclophosphamide_short$most_sensitive, new_brca_tcga_cyclophosphamide_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_cyclophosphamide_most_sensitive_1se <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_cyclophosphamide_short$most_sensitive, new_brca_tcga_cyclophosphamide_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_cyclophosphamide_least_sensitive_min <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_cyclophosphamide_short$least_sensitive, new_brca_tcga_cyclophosphamide_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_cyclophosphamide_least_sensitive_1se <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_cyclophosphamide_short$least_sensitive, new_brca_tcga_cyclophosphamide_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_cyclophosphamide_most_sensitive_min, brca_clinical_cyclophosphamide_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_cyclophosphamide_most_sensitive_1se, brca_clinical_cyclophosphamide_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_cyclophosphamide_least_sensitive_min, brca_clinical_cyclophosphamide_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_cyclophosphamide_least_sensitive_1se, brca_clinical_cyclophosphamide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_cyclophosphamide_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
#plot(perf3, col = 'blue', add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.48)', 'most_1se (0.61)', 'least_min (0.50)', 'least_1se (0.50)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BRCA W DOCETAXEL (24)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_docetaxel <- brca_clinical[which(brca_clinical$drug_name == 'docetaxel' | brca_clinical$drug_name == 'Docetaxel' | 
                                                 brca_clinical$drug_name == 'DOCETAXEL'), ]

brca_clinical_docetaxel$most_sensitive  <- ifelse(brca_clinical_docetaxel$PFS < quantile(brca_clinical_docetaxel$PFS, probs = .20), 1, 0)
brca_clinical_docetaxel$least_sensitive <- ifelse(brca_clinical_docetaxel$PFS > quantile(brca_clinical_docetaxel$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_docetaxel$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_docetaxel_short <- brca_clinical_docetaxel[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_docetaxel_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_docetaxel_short$most_sensitive, new_brca_tcga_docetaxel_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_docetaxel_short$most_sensitive, new_brca_tcga_docetaxel_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_docetaxel_short$least_sensitive, new_brca_tcga_docetaxel_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_docetaxel_short$least_sensitive, new_brca_tcga_docetaxel_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_docetaxel_most_sensitive_min, brca_clinical_docetaxel_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_docetaxel_most_sensitive_1se, brca_clinical_docetaxel_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_docetaxel_least_sensitive_min, brca_clinical_docetaxel_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_docetaxel_least_sensitive_1se, brca_clinical_docetaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_docetaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.67)', 'most_1se (0.66)', 'least_min (0.34)', 'least_1se (0.47)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BRCA W DOXORUBICIN (47)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_doxorubicin <- brca_clinical[which(brca_clinical$drug_name == 'doxorubicin' | brca_clinical$drug_name == 'Doxorubicin' | 
                                                   brca_clinical$drug_name == 'DOXORUBICIN' | brca_clinical$drug_name == 'Doxorubicinum'), ]

brca_clinical_doxorubicin$most_sensitive  <- ifelse(brca_clinical_doxorubicin$PFS < quantile(brca_clinical_doxorubicin$PFS, probs = .20), 1, 0)
brca_clinical_doxorubicin$least_sensitive <- ifelse(brca_clinical_doxorubicin$PFS > quantile(brca_clinical_doxorubicin$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_doxorubicin$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_doxorubicin_short <- brca_clinical_doxorubicin[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_doxorubicin_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_doxorubicin_most_sensitive_min <- predict(doxorubicin_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_doxorubicin_short$most_sensitive, new_brca_tcga_doxorubicin_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_doxorubicin_most_sensitive_1se <- predict(doxorubicin_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_doxorubicin_short$most_sensitive, new_brca_tcga_doxorubicin_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_doxorubicin_least_sensitive_min <- predict(doxorubicin_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_doxorubicin_short$least_sensitive, new_brca_tcga_doxorubicin_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_doxorubicin_least_sensitive_1se <- predict(doxorubicin_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_doxorubicin_short$least_sensitive, new_brca_tcga_doxorubicin_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_doxorubicin_most_sensitive_min, brca_clinical_doxorubicin_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_doxorubicin_most_sensitive_1se, brca_clinical_doxorubicin_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_doxorubicin_least_sensitive_min, brca_clinical_doxorubicin_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_doxorubicin_least_sensitive_1se, brca_clinical_doxorubicin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_doxorubicin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.47)', 'most_1se (0.53)', 'least_min (0.52)', 'least_1se (0.52)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BRCA W FLUOROURACIL (14)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_fluorouracil <- brca_clinical[which(brca_clinical$drug_name == '5-Fluorouracil' | brca_clinical$drug_name == '5-FU' | 
                                                    brca_clinical$drug_name == 'fluorouracil' | brca_clinical$drug_name == 'Fluorouracil'), ]

brca_clinical_fluorouracil$most_sensitive  <- ifelse(brca_clinical_fluorouracil$PFS < quantile(brca_clinical_fluorouracil$PFS, probs = .20), 1, 0)
brca_clinical_fluorouracil$least_sensitive <- ifelse(brca_clinical_fluorouracil$PFS > quantile(brca_clinical_fluorouracil$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_fluorouracil$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_fluorouracil_short <- brca_clinical_fluorouracil[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_fluorouracil_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_fluorouracil_short$most_sensitive, new_brca_tcga_fluorouracil_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_fluorouracil_short$most_sensitive, new_brca_tcga_fluorouracil_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_fluorouracil_short$least_sensitive, new_brca_tcga_fluorouracil_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_fluorouracil_short$least_sensitive, new_brca_tcga_fluorouracil_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_fluorouracil_most_sensitive_min, brca_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_fluorouracil_most_sensitive_1se, brca_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_fluorouracil_least_sensitive_min, brca_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_fluorouracil_least_sensitive_1se, brca_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.58)', 'most_1se (0.33)', 'least_min (0.39)', 'least_1se (0.36)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BRCA W GEMCITABINE (4)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_gemcitabine <- brca_clinical[which(brca_clinical$drug_name == 'Gemcitabine' | brca_clinical$drug_name == 'GEMZAR'), ]

brca_clinical_gemcitabine$most_sensitive  <- ifelse(brca_clinical_gemcitabine$PFS < quantile(brca_clinical_gemcitabine$PFS, probs = .20), 1, 0)
brca_clinical_gemcitabine$least_sensitive <- ifelse(brca_clinical_gemcitabine$PFS > quantile(brca_clinical_gemcitabine$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_gemcitabine$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_gemcitabine_short <- brca_clinical_gemcitabine[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_gemcitabine_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_gemcitabine_short$most_sensitive, new_brca_tcga_gemcitabine_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_gemcitabine_short$most_sensitive, new_brca_tcga_gemcitabine_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_gemcitabine_short$least_sensitive, new_brca_tcga_gemcitabine_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_gemcitabine_short$least_sensitive, new_brca_tcga_gemcitabine_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_gemcitabine_most_sensitive_min, brca_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_gemcitabine_most_sensitive_1se, brca_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_gemcitabine_least_sensitive_min, brca_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_gemcitabine_least_sensitive_1se, brca_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[7], lwd = 2, lty = 2, main = 'BRCA treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[7], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 1.0)', 'resistant model (AUC = 1.0)'), cex = 0.8, lty = c(1,2), lwd = 2, col = colors_i_need[7], bty = 'n')
dev.off()

# BRCA W METHOTREXATE (9)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_methotrexate <- brca_clinical[which(brca_clinical$drug_name == 'Methotrexate' | brca_clinical$drug_name == 'METHOTREXATE'), ]

brca_clinical_methotrexate$most_sensitive  <- ifelse(brca_clinical_methotrexate$PFS < quantile(brca_clinical_methotrexate$PFS, probs = .20), 1, 0)
brca_clinical_methotrexate$least_sensitive <- ifelse(brca_clinical_methotrexate$PFS > quantile(brca_clinical_methotrexate$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_methotrexate$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_methotrexate_short <- brca_clinical_methotrexate[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_methotrexate_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_methotrexate_most_sensitive_min <- predict(methotrexate_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_methotrexate_short$most_sensitive, new_brca_tcga_methotrexate_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_methotrexate_most_sensitive_1se <- predict(methotrexate_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_methotrexate_short$most_sensitive, new_brca_tcga_methotrexate_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_methotrexate_least_sensitive_min <- predict(methotrexate_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_methotrexate_short$least_sensitive, new_brca_tcga_methotrexate_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_methotrexate_least_sensitive_1se <- predict(methotrexate_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_methotrexate_short$least_sensitive, new_brca_tcga_methotrexate_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_methotrexate_most_sensitive_min, brca_clinical_methotrexate_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_methotrexate_most_sensitive_1se, brca_clinical_methotrexate_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_methotrexate_least_sensitive_min, brca_clinical_methotrexate_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_methotrexate_least_sensitive_1se, brca_clinical_methotrexate_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_methotrexate_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# BRCA W PACLITAXEL (36)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_paclitaxel <- brca_clinical[which(brca_clinical$drug_name == 'paclitaxel' | brca_clinical$drug_name == 'Paclitaxel' | 
                                                  brca_clinical$drug_name == 'PACLITAXEL' | brca_clinical$drug_name == 'Paclitaxel (Protein-Bound)'), ]

brca_clinical_paclitaxel$most_sensitive  <- ifelse(brca_clinical_paclitaxel$PFS < quantile(brca_clinical_paclitaxel$PFS, probs = .20), 1, 0)
brca_clinical_paclitaxel$least_sensitive <- ifelse(brca_clinical_paclitaxel$PFS > quantile(brca_clinical_paclitaxel$PFS, probs = .80), 1, 0)

brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_paclitaxel$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_paclitaxel_short <- brca_clinical_paclitaxel[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_paclitaxel_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)

rm(brca_tcga_most_min_auc)
rm(brca_tcga_most_1se_auc)
rm(brca_tcga_least_min_auc)
rm(brca_tcga_least_1se_auc)

new_brca_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_paclitaxel_short$most_sensitive, new_brca_tcga_paclitaxel_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_paclitaxel_short$most_sensitive, new_brca_tcga_paclitaxel_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_paclitaxel_short$least_sensitive, new_brca_tcga_paclitaxel_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_least_1se_auc <- auc(brca_clinical_paclitaxel_short$least_sensitive, new_brca_tcga_paclitaxel_least_sensitive_1se)
brca_tcga_least_1se_auc <- round(brca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_brca_tcga_paclitaxel_most_sensitive_min, brca_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_brca_tcga_paclitaxel_most_sensitive_1se, brca_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_brca_tcga_paclitaxel_least_sensitive_min, brca_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_brca_tcga_paclitaxel_least_sensitive_1se, brca_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/brca_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# CESC W CARBOPLATIN (7)
cesc_clinical <- read.csv('Processed_Clinical_Data/cesc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(cesc_clinical$most_sensitive)
cesc_clinical <- cesc_clinical[!na_idx, ]
table(cesc_clinical$drug_name)
cesc_clinical_carboplatin <- cesc_clinical[which(cesc_clinical$drug_name == 'Carboplatin'), ]

cesc_clinical_carboplatin$most_sensitive  <- ifelse(cesc_clinical_carboplatin$PFS < quantile(cesc_clinical_carboplatin$PFS, probs = .20), 1, 0)
cesc_clinical_carboplatin$least_sensitive <- ifelse(cesc_clinical_carboplatin$PFS > quantile(cesc_clinical_carboplatin$PFS, probs = .80), 1, 0)

cesc_gene <- read.csv('Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(cesc_gene) <- gsub('\\.', '-', colnames(cesc_gene))
cesc_matching_idx <- cesc_clinical_carboplatin$submitter_id.samples %in% colnames(cesc_gene)
cesc_clinical_carboplatin_short <- cesc_clinical_carboplatin[cesc_matching_idx, ]
cesc_matching_idx <- colnames(cesc_gene) %in% cesc_clinical_carboplatin_short$submitter_id.samples
cesc_gene_short <- cesc_gene[, cesc_matching_idx]
cesc_gene_short <- t(cesc_gene_short)
cesc_gene_short_scaled <- apply(cesc_gene_short, 2, scale)

rm(cesc_tcga_most_min_auc)
rm(cesc_tcga_most_1se_auc)
rm(cesc_tcga_least_min_auc)
rm(cesc_tcga_least_1se_auc)

new_cesc_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_most_min_auc <- auc(cesc_clinical_carboplatin_short$most_sensitive, new_cesc_tcga_carboplatin_most_sensitive_min)
cesc_tcga_most_min_auc <- round(cesc_tcga_most_min_auc, digits = 2)

new_cesc_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_most_1se_auc <- auc(cesc_clinical_carboplatin_short$most_sensitive, new_cesc_tcga_carboplatin_most_sensitive_1se)
cesc_tcga_most_1se_auc <- round(cesc_tcga_most_1se_auc, digits = 2)

new_cesc_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_least_min_auc <- auc(cesc_clinical_carboplatin_short$least_sensitive, new_cesc_tcga_carboplatin_least_sensitive_min)
cesc_tcga_least_min_auc <- round(cesc_tcga_least_min_auc, digits = 2)

new_cesc_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_least_1se_auc <- auc(cesc_clinical_carboplatin_short$least_sensitive, new_cesc_tcga_carboplatin_least_sensitive_1se)
cesc_tcga_least_1se_auc <- round(cesc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_cesc_tcga_carboplatin_most_sensitive_min, cesc_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_cesc_tcga_carboplatin_most_sensitive_1se, cesc_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_cesc_tcga_carboplatin_least_sensitive_min, cesc_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_cesc_tcga_carboplatin_least_sensitive_1se, cesc_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/cesc_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# CESC W FLUOROURACIL (5)
cesc_clinical <- read.csv('Processed_Clinical_Data/cesc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(cesc_clinical$most_sensitive)
cesc_clinical <- cesc_clinical[!na_idx, ]
table(cesc_clinical$drug_name)
cesc_clinical_fluorouracil <- cesc_clinical[which(cesc_clinical$drug_name == '5-Fluorouracil' | cesc_clinical$drug_name == '5FU' | 
                                                    cesc_clinical$drug_name == 'Fluorouracil' | cesc_clinical$drug_name == 'Fluorouracil (5-FU)' | 
                                                    cesc_clinical$drug_name == 'Fluoruracil'), ]

cesc_clinical_fluorouracil$most_sensitive  <- ifelse(cesc_clinical_fluorouracil$PFS < quantile(cesc_clinical_fluorouracil$PFS, probs = .20), 1, 0)
cesc_clinical_fluorouracil$least_sensitive <- ifelse(cesc_clinical_fluorouracil$PFS > quantile(cesc_clinical_fluorouracil$PFS, probs = .80), 1, 0)

cesc_gene <- read.csv('Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(cesc_gene) <- gsub('\\.', '-', colnames(cesc_gene))
cesc_matching_idx <- cesc_clinical_fluorouracil$submitter_id.samples %in% colnames(cesc_gene)
cesc_clinical_fluorouracil_short <- cesc_clinical_fluorouracil[cesc_matching_idx, ]
cesc_matching_idx <- colnames(cesc_gene) %in% cesc_clinical_fluorouracil_short$submitter_id.samples
cesc_gene_short <- cesc_gene[, cesc_matching_idx]
cesc_gene_short <- t(cesc_gene_short)
cesc_gene_short_scaled <- apply(cesc_gene_short, 2, scale)

rm(cesc_tcga_most_min_auc)
rm(cesc_tcga_most_1se_auc)
rm(cesc_tcga_least_min_auc)
rm(cesc_tcga_least_1se_auc)

new_cesc_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_most_min_auc <- auc(cesc_clinical_fluorouracil_short$most_sensitive, new_cesc_tcga_fluorouracil_most_sensitive_min)
cesc_tcga_most_min_auc <- round(cesc_tcga_most_min_auc, digits = 2)

new_cesc_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_most_1se_auc <- auc(cesc_clinical_fluorouracil_short$most_sensitive, new_cesc_tcga_fluorouracil_most_sensitive_1se)
cesc_tcga_most_1se_auc <- round(cesc_tcga_most_1se_auc, digits = 2)

new_cesc_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_least_min_auc <- auc(cesc_clinical_fluorouracil_short$least_sensitive, new_cesc_tcga_fluorouracil_least_sensitive_min)
cesc_tcga_least_min_auc <- round(cesc_tcga_least_min_auc, digits = 2)

new_cesc_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_least_1se_auc <- auc(cesc_clinical_fluorouracil_short$least_sensitive, new_cesc_tcga_fluorouracil_least_sensitive_1se)
cesc_tcga_least_1se_auc <- round(cesc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_cesc_tcga_fluorouracil_most_sensitive_min, cesc_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_cesc_tcga_fluorouracil_most_sensitive_1se, cesc_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_cesc_tcga_fluorouracil_least_sensitive_min, cesc_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_cesc_tcga_fluorouracil_least_sensitive_1se, cesc_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/cesc_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# CESC W PACLITAXEL (8)
cesc_clinical <- read.csv('Processed_Clinical_Data/cesc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(cesc_clinical$most_sensitive)
cesc_clinical <- cesc_clinical[!na_idx, ]
table(cesc_clinical$drug_name)
cesc_clinical_paclitaxel <- cesc_clinical[which(cesc_clinical$drug_name == 'Paclitaxel'), ]

cesc_clinical_paclitaxel$most_sensitive  <- ifelse(cesc_clinical_paclitaxel$PFS < quantile(cesc_clinical_paclitaxel$PFS, probs = .20), 1, 0)
cesc_clinical_paclitaxel$least_sensitive <- ifelse(cesc_clinical_paclitaxel$PFS > quantile(cesc_clinical_paclitaxel$PFS, probs = .80), 1, 0)

cesc_gene <- read.csv('Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(cesc_gene) <- gsub('\\.', '-', colnames(cesc_gene))
cesc_matching_idx <- cesc_clinical_paclitaxel$submitter_id.samples %in% colnames(cesc_gene)
cesc_clinical_paclitaxel_short <- cesc_clinical_paclitaxel[cesc_matching_idx, ]
cesc_matching_idx <- colnames(cesc_gene) %in% cesc_clinical_paclitaxel_short$submitter_id.samples
cesc_gene_short <- cesc_gene[, cesc_matching_idx]
cesc_gene_short <- t(cesc_gene_short)
cesc_gene_short_scaled <- apply(cesc_gene_short, 2, scale)

rm(cesc_tcga_most_min_auc)
rm(cesc_tcga_most_1se_auc)
rm(cesc_tcga_least_min_auc)
rm(cesc_tcga_least_1se_auc)

new_cesc_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_most_min_auc <- auc(cesc_clinical_paclitaxel_short$most_sensitive, new_cesc_tcga_paclitaxel_most_sensitive_min)
cesc_tcga_most_min_auc <- round(cesc_tcga_most_min_auc, digits = 2)

new_cesc_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_most_1se_auc <- auc(cesc_clinical_paclitaxel_short$most_sensitive, new_cesc_tcga_paclitaxel_most_sensitive_1se)
cesc_tcga_most_1se_auc <- round(cesc_tcga_most_1se_auc, digits = 2)

new_cesc_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_least_min_auc <- auc(cesc_clinical_paclitaxel_short$least_sensitive, new_cesc_tcga_paclitaxel_least_sensitive_min)
cesc_tcga_least_min_auc <- round(cesc_tcga_least_min_auc, digits = 2)

new_cesc_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_least_1se_auc <- auc(cesc_clinical_paclitaxel_short$least_sensitive, new_cesc_tcga_paclitaxel_least_sensitive_1se)
cesc_tcga_least_1se_auc <- round(cesc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_cesc_tcga_paclitaxel_most_sensitive_min, cesc_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_cesc_tcga_paclitaxel_most_sensitive_1se, cesc_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_cesc_tcga_paclitaxel_least_sensitive_min, cesc_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_cesc_tcga_paclitaxel_least_sensitive_1se, cesc_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/cesc_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# CHOL W GEMCITABINE (8)
chol_clinical <- read.csv('Processed_Clinical_Data/chol_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(chol_clinical$most_sensitive)
chol_clinical <- chol_clinical[!na_idx, ]
table(chol_clinical$drug_name)
chol_clinical_gemcitabine <- chol_clinical[which(chol_clinical$drug_name == 'gemcitabine' | chol_clinical$drug_name == 'Gemcitabine' | 
                                                   chol_clinical$drug_name == 'GEMCITABINE' | chol_clinical$drug_name == 'Gemzar'), ]

chol_clinical_gemcitabine$most_sensitive  <- ifelse(chol_clinical_gemcitabine$PFS < quantile(chol_clinical_gemcitabine$PFS, probs = .20), 1, 0)
chol_clinical_gemcitabine$least_sensitive <- ifelse(chol_clinical_gemcitabine$PFS > quantile(chol_clinical_gemcitabine$PFS, probs = .80), 1, 0)

chol_gene <- read.csv('Processed_Gene_Expression/chol_tcga_rna_seq_processed.csv', row.names = 1)
colnames(chol_gene) <- gsub('\\.', '-', colnames(chol_gene))
chol_matching_idx <- chol_clinical_gemcitabine$submitter_id.samples %in% colnames(chol_gene)
chol_clinical_gemcitabine_short <- chol_clinical_gemcitabine[chol_matching_idx, ]
chol_matching_idx <- colnames(chol_gene) %in% chol_clinical_gemcitabine_short$submitter_id.samples
chol_gene_short <- chol_gene[, chol_matching_idx]
chol_gene_short <- t(chol_gene_short)
chol_gene_short_scaled <- apply(chol_gene_short, 2, scale)

rm(chol_tcga_most_min_auc)
rm(chol_tcga_most_1se_auc)
rm(chol_tcga_least_min_auc)
rm(chol_tcga_least_1se_auc)

new_chol_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
chol_tcga_most_min_auc <- auc(chol_clinical_gemcitabine_short$most_sensitive, new_chol_tcga_gemcitabine_most_sensitive_min)
chol_tcga_most_min_auc <- round(chol_tcga_most_min_auc, digits = 2)

new_chol_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
chol_tcga_most_1se_auc <- auc(chol_clinical_gemcitabine_short$most_sensitive, new_chol_tcga_gemcitabine_most_sensitive_1se)
chol_tcga_most_1se_auc <- round(chol_tcga_most_1se_auc, digits = 2)

new_chol_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
chol_tcga_least_min_auc <- auc(chol_clinical_gemcitabine_short$least_sensitive, new_chol_tcga_gemcitabine_least_sensitive_min)
chol_tcga_least_min_auc <- round(chol_tcga_least_min_auc, digits = 2)

new_chol_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = chol_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
chol_tcga_least_1se_auc <- auc(chol_clinical_gemcitabine_short$least_sensitive, new_chol_tcga_gemcitabine_least_sensitive_1se)
chol_tcga_least_1se_auc <- round(chol_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_chol_tcga_gemcitabine_most_sensitive_min, chol_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_chol_tcga_gemcitabine_most_sensitive_1se, chol_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_chol_tcga_gemcitabine_least_sensitive_min, chol_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_chol_tcga_gemcitabine_least_sensitive_1se, chol_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/chol_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# COAD W FLUOROURACIL (56)
coad_clinical <- read.csv('Processed_Clinical_Data/coad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(coad_clinical$most_sensitive)
coad_clinical <- coad_clinical[!na_idx, ]
table(coad_clinical$drug_name)
coad_clinical_fluorouracil <- coad_clinical[which(coad_clinical$drug_name == '5- FU' | coad_clinical$drug_name == '5-Fluorouracil' | 
                                                    coad_clinical$drug_name == '5-FU' | coad_clinical$drug_name == '5 FU' | coad_clinical$drug_name == '5FU' | 
                                                    coad_clinical$drug_name == 'fluorouracil' | coad_clinical$drug_name == 'Fluorouracil'), ]

coad_clinical_fluorouracil$most_sensitive  <- ifelse(coad_clinical_fluorouracil$PFS < quantile(coad_clinical_fluorouracil$PFS, probs = .20), 1, 0)
coad_clinical_fluorouracil$least_sensitive <- ifelse(coad_clinical_fluorouracil$PFS > quantile(coad_clinical_fluorouracil$PFS, probs = .80), 1, 0)

coad_gene <- read.csv('Processed_Gene_Expression/coad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(coad_gene) <- gsub('\\.', '-', colnames(coad_gene))
coad_matching_idx <- coad_clinical_fluorouracil$submitter_id.samples %in% colnames(coad_gene)
coad_clinical_fluorouracil_short <- coad_clinical_fluorouracil[coad_matching_idx, ]
coad_matching_idx <- colnames(coad_gene) %in% coad_clinical_fluorouracil_short$submitter_id.samples
coad_gene_short <- coad_gene[, coad_matching_idx]
coad_gene_short <- t(coad_gene_short)
coad_gene_short_scaled <- apply(coad_gene_short, 2, scale)

rm(coad_tcga_most_min_auc)
rm(coad_tcga_most_1se_auc)
rm(coad_tcga_least_min_auc)
rm(coad_tcga_least_1se_auc)

new_coad_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
coad_tcga_most_min_auc <- auc(coad_clinical_fluorouracil_short$most_sensitive, new_coad_tcga_fluorouracil_most_sensitive_min)
coad_tcga_most_min_auc <- round(coad_tcga_most_min_auc, digits = 2)

new_coad_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
coad_tcga_most_1se_auc <- auc(coad_clinical_fluorouracil_short$most_sensitive, new_coad_tcga_fluorouracil_most_sensitive_1se)
coad_tcga_most_1se_auc <- round(coad_tcga_most_1se_auc, digits = 2)

new_coad_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
coad_tcga_least_min_auc <- auc(coad_clinical_fluorouracil_short$least_sensitive, new_coad_tcga_fluorouracil_least_sensitive_min)
coad_tcga_least_min_auc <- round(coad_tcga_least_min_auc, digits = 2)

new_coad_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
coad_tcga_least_1se_auc <- auc(coad_clinical_fluorouracil_short$least_sensitive, new_coad_tcga_fluorouracil_least_sensitive_1se)
coad_tcga_least_1se_auc <- round(coad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_coad_tcga_fluorouracil_most_sensitive_min, coad_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_coad_tcga_fluorouracil_most_sensitive_1se, coad_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_coad_tcga_fluorouracil_least_sensitive_min, coad_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_coad_tcga_fluorouracil_least_sensitive_1se, coad_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/coad_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[6], lwd = 2, lty = 2, main = 'COAD treated with fluorouracil')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[6], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.62)', 'resistant model (AUC = 0.62)'), cex = 0.8, lty = c(1,2), col = colors_i_need[6], lwd = 2, bty = 'n')
dev.off()

# COAD W OXALIPLATIN (34)
coad_clinical <- read.csv('Processed_Clinical_Data/coad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(coad_clinical$most_sensitive)
coad_clinical <- coad_clinical[!na_idx, ]
table(coad_clinical$drug_name)
coad_clinical_oxalaplatin <- coad_clinical[which(coad_clinical$drug_name == 'oxaliplatin' | coad_clinical$drug_name == 'Oxaliplatin'), ]

coad_clinical_oxalaplatin$most_sensitive  <- ifelse(coad_clinical_oxalaplatin$PFS < quantile(coad_clinical_oxalaplatin$PFS, probs = .20), 1, 0)
coad_clinical_oxalaplatin$least_sensitive <- ifelse(coad_clinical_oxalaplatin$PFS > quantile(coad_clinical_oxalaplatin$PFS, probs = .80), 1, 0)

coad_gene <- read.csv('Processed_Gene_Expression/coad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(coad_gene) <- gsub('\\.', '-', colnames(coad_gene))
coad_matching_idx <- coad_clinical_oxalaplatin$submitter_id.samples %in% colnames(coad_gene)
coad_clinical_oxalaplatin_short <- coad_clinical_oxalaplatin[coad_matching_idx, ]
coad_matching_idx <- colnames(coad_gene) %in% coad_clinical_oxalaplatin_short$submitter_id.samples
coad_gene_short <- coad_gene[, coad_matching_idx]
coad_gene_short <- t(coad_gene_short)
coad_gene_short_scaled <- apply(coad_gene_short, 2, scale)

rm(coad_tcga_most_min_auc)
rm(coad_tcga_most_1se_auc)
rm(coad_tcga_least_min_auc)
rm(coad_tcga_least_1se_auc)

new_coad_tcga_oxalaplatin_most_sensitive_min <- predict(oxalaplatin_ccle_most_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
coad_tcga_most_min_auc <- auc(coad_clinical_oxalaplatin_short$most_sensitive, new_coad_tcga_oxalaplatin_most_sensitive_min)
coad_tcga_most_min_auc <- round(coad_tcga_most_min_auc, digits = 2)

new_coad_tcga_oxalaplatin_most_sensitive_1se <- predict(oxalaplatin_ccle_most_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
coad_tcga_most_1se_auc <- auc(coad_clinical_oxalaplatin_short$most_sensitive, new_coad_tcga_oxalaplatin_most_sensitive_1se)
coad_tcga_most_1se_auc <- round(coad_tcga_most_1se_auc, digits = 2)

new_coad_tcga_oxalaplatin_least_sensitive_min <- predict(oxalaplatin_ccle_least_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
coad_tcga_least_min_auc <- auc(coad_clinical_oxalaplatin_short$least_sensitive, new_coad_tcga_oxalaplatin_least_sensitive_min)
coad_tcga_least_min_auc <- round(coad_tcga_least_min_auc, digits = 2)

new_coad_tcga_oxalaplatin_least_sensitive_1se <- predict(oxalaplatin_ccle_least_fit_elnet, newx = coad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
coad_tcga_least_1se_auc <- auc(coad_clinical_oxalaplatin_short$least_sensitive, new_coad_tcga_oxalaplatin_least_sensitive_1se)
coad_tcga_least_1se_auc <- round(coad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_coad_tcga_oxalaplatin_most_sensitive_min, coad_clinical_oxalaplatin_short$most_sensitive)
pred2 <- prediction(new_coad_tcga_oxalaplatin_most_sensitive_1se, coad_clinical_oxalaplatin_short$most_sensitive)
pred3 <- prediction(new_coad_tcga_oxalaplatin_least_sensitive_min, coad_clinical_oxalaplatin_short$least_sensitive)
pred4 <- prediction(new_coad_tcga_oxalaplatin_least_sensitive_1se, coad_clinical_oxalaplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/coad_oxalaplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# DLBC W CYCLOPHOSPHAMIDE (4)
dlbc_clinical <- read.csv('Processed_Clinical_Data/dlbc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(dlbc_clinical$most_sensitive)
dlbc_clinical <- dlbc_clinical[!na_idx, ]
table(dlbc_clinical$drug_name)
dlbc_clinical_cyclophosphamide <- dlbc_clinical[which(dlbc_clinical$drug_name == 'Cyclophosphamide'), ]

dlbc_clinical_cyclophosphamide$most_sensitive  <- ifelse(dlbc_clinical_cyclophosphamide$PFS < quantile(dlbc_clinical_cyclophosphamide$PFS, probs = .20), 1, 0)
dlbc_clinical_cyclophosphamide$least_sensitive <- ifelse(dlbc_clinical_cyclophosphamide$PFS > quantile(dlbc_clinical_cyclophosphamide$PFS, probs = .80), 1, 0)

dlbc_gene <- read.csv('Processed_Gene_Expression/dlbc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(dlbc_gene) <- gsub('\\.', '-', colnames(dlbc_gene))
dlbc_matching_idx <- dlbc_clinical_cyclophosphamide$submitter_id.samples %in% colnames(dlbc_gene)
dlbc_clinical_cyclophosphamide_short <- dlbc_clinical_cyclophosphamide[dlbc_matching_idx, ]
dlbc_matching_idx <- colnames(dlbc_gene) %in% dlbc_clinical_cyclophosphamide_short$submitter_id.samples
dlbc_gene_short <- dlbc_gene[, dlbc_matching_idx]
dlbc_gene_short <- t(dlbc_gene_short)
dlbc_gene_short_scaled <- apply(dlbc_gene_short, 2, scale)

rm(dlbc_tcga_most_min_auc)
rm(dlbc_tcga_most_1se_auc)
rm(dlbc_tcga_least_min_auc)
rm(dlbc_tcga_least_1se_auc)

new_dlbc_tcga_cyclophosphamide_most_sensitive_min <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
dlbc_tcga_most_min_auc <- auc(dlbc_clinical_cyclophosphamide_short$most_sensitive, new_dlbc_tcga_cyclophosphamide_most_sensitive_min)
dlbc_tcga_most_min_auc <- round(dlbc_tcga_most_min_auc, digits = 2)

new_dlbc_tcga_cyclophosphamide_most_sensitive_1se <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
dlbc_tcga_most_1se_auc <- auc(dlbc_clinical_cyclophosphamide_short$most_sensitive, new_dlbc_tcga_cyclophosphamide_most_sensitive_1se)
dlbc_tcga_most_1se_auc <- round(dlbc_tcga_most_1se_auc, digits = 2)

new_dlbc_tcga_cyclophosphamide_least_sensitive_min <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
dlbc_tcga_least_min_auc <- auc(dlbc_clinical_cyclophosphamide_short$least_sensitive, new_dlbc_tcga_cyclophosphamide_least_sensitive_min)
dlbc_tcga_least_min_auc <- round(dlbc_tcga_least_min_auc, digits = 2)

new_dlbc_tcga_cyclophosphamide_least_sensitive_1se <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
dlbc_tcga_least_1se_auc <- auc(dlbc_clinical_cyclophosphamide_short$least_sensitive, new_dlbc_tcga_cyclophosphamide_least_sensitive_1se)
dlbc_tcga_least_1se_auc <- round(dlbc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_dlbc_tcga_cyclophosphamide_most_sensitive_min, dlbc_clinical_cyclophosphamide_short$most_sensitive)
pred2 <- prediction(new_dlbc_tcga_cyclophosphamide_most_sensitive_1se, dlbc_clinical_cyclophosphamide_short$most_sensitive)
pred3 <- prediction(new_dlbc_tcga_cyclophosphamide_least_sensitive_min, dlbc_clinical_cyclophosphamide_short$least_sensitive)
pred4 <- prediction(new_dlbc_tcga_cyclophosphamide_least_sensitive_1se, dlbc_clinical_cyclophosphamide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/dlbc_cyclophosphamide_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# DLBC W VINCRISTINE (4)
dlbc_clinical <- read.csv('Processed_Clinical_Data/dlbc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(dlbc_clinical$most_sensitive)
dlbc_clinical <- dlbc_clinical[!na_idx, ]
table(dlbc_clinical$drug_name)
dlbc_clinical_vincristine <- dlbc_clinical[which(dlbc_clinical$drug_name == 'vincristina' | dlbc_clinical$drug_name == 'Vincristine'), ]

dlbc_clinical_vincristine$most_sensitive  <- ifelse(dlbc_clinical_vincristine$PFS < quantile(dlbc_clinical_vincristine$PFS, probs = .20), 1, 0)
dlbc_clinical_vincristine$least_sensitive <- ifelse(dlbc_clinical_vincristine$PFS > quantile(dlbc_clinical_vincristine$PFS, probs = .80), 1, 0)

dlbc_gene <- read.csv('Processed_Gene_Expression/dlbc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(dlbc_gene) <- gsub('\\.', '-', colnames(dlbc_gene))
dlbc_matching_idx <- dlbc_clinical_vincristine$submitter_id.samples %in% colnames(dlbc_gene)
dlbc_clinical_vincristine_short <- dlbc_clinical_vincristine[dlbc_matching_idx, ]
dlbc_matching_idx <- colnames(dlbc_gene) %in% dlbc_clinical_vincristine_short$submitter_id.samples
dlbc_gene_short <- dlbc_gene[, dlbc_matching_idx]
dlbc_gene_short <- t(dlbc_gene_short)
dlbc_gene_short_scaled <- apply(dlbc_gene_short, 2, scale)

rm(dlbc_tcga_most_min_auc)
rm(dlbc_tcga_most_1se_auc)
rm(dlbc_tcga_least_min_auc)
rm(dlbc_tcga_least_1se_auc)

new_dlbc_tcga_vincristine_most_sensitive_min <- predict(vincristine_ccle_most_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
dlbc_tcga_most_min_auc <- auc(dlbc_clinical_vincristine_short$most_sensitive, new_dlbc_tcga_vincristine_most_sensitive_min)
dlbc_tcga_most_min_auc <- round(dlbc_tcga_most_min_auc, digits = 2)

new_dlbc_tcga_vincristine_most_sensitive_1se <- predict(vincristine_ccle_most_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
dlbc_tcga_most_1se_auc <- auc(dlbc_clinical_vincristine_short$most_sensitive, new_dlbc_tcga_vincristine_most_sensitive_1se)
dlbc_tcga_most_1se_auc <- round(dlbc_tcga_most_1se_auc, digits = 2)

new_dlbc_tcga_vincristine_least_sensitive_min <- predict(vincristine_ccle_least_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
dlbc_tcga_least_min_auc <- auc(dlbc_clinical_vincristine_short$least_sensitive, new_dlbc_tcga_vincristine_least_sensitive_min)
dlbc_tcga_least_min_auc <- round(dlbc_tcga_least_min_auc, digits = 2)

new_dlbc_tcga_vincristine_least_sensitive_1se <- predict(vincristine_ccle_least_fit_elnet, newx = dlbc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
dlbc_tcga_least_1se_auc <- auc(dlbc_clinical_vincristine_short$least_sensitive, new_dlbc_tcga_vincristine_least_sensitive_1se)
dlbc_tcga_least_1se_auc <- round(dlbc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_dlbc_tcga_vincristine_most_sensitive_min, dlbc_clinical_vincristine_short$most_sensitive)
pred2 <- prediction(new_dlbc_tcga_vincristine_most_sensitive_1se, dlbc_clinical_vincristine_short$most_sensitive)
pred3 <- prediction(new_dlbc_tcga_vincristine_least_sensitive_min, dlbc_clinical_vincristine_short$least_sensitive)
pred4 <- prediction(new_dlbc_tcga_vincristine_least_sensitive_1se, dlbc_clinical_vincristine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/dlbc_vincristine_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# ESCA W FLUOROURACIL (8)
esca_clinical <- read.csv('Processed_Clinical_Data/esca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(esca_clinical$most_sensitive)
esca_clinical <- esca_clinical[!na_idx, ]
table(esca_clinical$drug_name)
esca_clinical_fluorouracil <- esca_clinical[which(esca_clinical$drug_name == '5 fluorouracil' | esca_clinical$drug_name == '5 FU' | 
                                                    esca_clinical$drug_name == '5FU' | esca_clinical$drug_name == 'fluorouracil' | 
                                                    esca_clinical$drug_name == 'Fluorouracil'), ]

esca_clinical_fluorouracil$most_sensitive  <- ifelse(esca_clinical_fluorouracil$PFS < quantile(esca_clinical_fluorouracil$PFS, probs = .20), 1, 0)
esca_clinical_fluorouracil$least_sensitive <- ifelse(esca_clinical_fluorouracil$PFS > quantile(esca_clinical_fluorouracil$PFS, probs = .80), 1, 0)

esca_gene <- read.csv('Processed_Gene_Expression/esca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(esca_gene) <- gsub('\\.', '-', colnames(esca_gene))
esca_matching_idx <- esca_clinical_fluorouracil$submitter_id.samples %in% colnames(esca_gene)
esca_clinical_fluorouracil_short <- esca_clinical_fluorouracil[esca_matching_idx, ]
esca_matching_idx <- colnames(esca_gene) %in% esca_clinical_fluorouracil_short$submitter_id.samples
esca_gene_short <- esca_gene[, esca_matching_idx]
esca_gene_short <- t(esca_gene_short)
esca_gene_short_scaled <- apply(esca_gene_short, 2, scale)

rm(esca_tcga_most_min_auc)
rm(esca_tcga_most_1se_auc)
rm(esca_tcga_least_min_auc)
rm(esca_tcga_least_1se_auc)

new_esca_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
esca_tcga_most_min_auc <- auc(esca_clinical_fluorouracil_short$most_sensitive, new_esca_tcga_fluorouracil_most_sensitive_min)
esca_tcga_most_min_auc <- round(esca_tcga_most_min_auc, digits = 2)

new_esca_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
esca_tcga_most_1se_auc <- auc(esca_clinical_fluorouracil_short$most_sensitive, new_esca_tcga_fluorouracil_most_sensitive_1se)
esca_tcga_most_1se_auc <- round(esca_tcga_most_1se_auc, digits = 2)

new_esca_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
esca_tcga_least_min_auc <- auc(esca_clinical_fluorouracil_short$least_sensitive, new_esca_tcga_fluorouracil_least_sensitive_min)
esca_tcga_least_min_auc <- round(esca_tcga_least_min_auc, digits = 2)

new_esca_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = esca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
esca_tcga_least_1se_auc <- auc(esca_clinical_fluorouracil_short$least_sensitive, new_esca_tcga_fluorouracil_least_sensitive_1se)
esca_tcga_least_1se_auc <- round(esca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_esca_tcga_fluorouracil_most_sensitive_min, esca_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_esca_tcga_fluorouracil_most_sensitive_1se, esca_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_esca_tcga_fluorouracil_least_sensitive_min, esca_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_esca_tcga_fluorouracil_least_sensitive_1se, esca_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/esca_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[6], lwd = 2, lty = 2, main = 'ESCA treated with fluorouracil')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[6], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.92)', 'resistant model (AUC = 0.58)'), cex = 0.8, lwd = 2, lty = c(1,2), col = colors_i_need[6], bty = 'n')
dev.off()

# HNSC W CARBOPLATIN (30)
hnsc_clinical <- read.csv('Processed_Clinical_Data/hnsc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(hnsc_clinical$most_sensitive)
hnsc_clinical <- hnsc_clinical[!na_idx, ]
table(hnsc_clinical$drug_name)
hnsc_clinical_carboplatin <- hnsc_clinical[which(hnsc_clinical$drug_name == 'carboplatin' | hnsc_clinical$drug_name == 'Carboplatin' | 
                                                   hnsc_clinical$drug_name == 'CARBOPLATIN' | hnsc_clinical$drug_name == 'Carboplatinum'), ]

hnsc_clinical_carboplatin$most_sensitive  <- ifelse(hnsc_clinical_carboplatin$PFS < quantile(hnsc_clinical_carboplatin$PFS, probs = .20), 1, 0)
hnsc_clinical_carboplatin$least_sensitive <- ifelse(hnsc_clinical_carboplatin$PFS > quantile(hnsc_clinical_carboplatin$PFS, probs = .80), 1, 0)

hnsc_gene <- read.csv('Processed_Gene_Expression/hnsc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(hnsc_gene) <- gsub('\\.', '-', colnames(hnsc_gene))
hnsc_matching_idx <- hnsc_clinical_carboplatin$submitter_id.samples %in% colnames(hnsc_gene)
hnsc_clinical_carboplatin_short <- hnsc_clinical_carboplatin[hnsc_matching_idx, ]
hnsc_matching_idx <- colnames(hnsc_gene) %in% hnsc_clinical_carboplatin_short$submitter_id.samples
hnsc_gene_short <- hnsc_gene[, hnsc_matching_idx]
hnsc_gene_short <- t(hnsc_gene_short)
hnsc_gene_short_scaled <- apply(hnsc_gene_short, 2, scale)

rm(hnsc_tcga_most_min_auc)
rm(hnsc_tcga_most_1se_auc)
rm(hnsc_tcga_least_min_auc)
rm(hnsc_tcga_least_1se_auc)

new_hnsc_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_min_auc <- auc(hnsc_clinical_carboplatin_short$most_sensitive, new_hnsc_tcga_carboplatin_most_sensitive_min)
hnsc_tcga_most_min_auc <- round(hnsc_tcga_most_min_auc, digits = 2)

new_hnsc_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_1se_auc <- auc(hnsc_clinical_carboplatin_short$most_sensitive, new_hnsc_tcga_carboplatin_most_sensitive_1se)
hnsc_tcga_most_1se_auc <- round(hnsc_tcga_most_1se_auc, digits = 2)

new_hnsc_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_min_auc <- auc(hnsc_clinical_carboplatin_short$least_sensitive, new_hnsc_tcga_carboplatin_least_sensitive_min)
hnsc_tcga_least_min_auc <- round(hnsc_tcga_least_min_auc, digits = 2)

new_hnsc_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_1se_auc <- auc(hnsc_clinical_carboplatin_short$least_sensitive, new_hnsc_tcga_carboplatin_least_sensitive_1se)
hnsc_tcga_least_1se_auc <- round(hnsc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_hnsc_tcga_carboplatin_most_sensitive_min, hnsc_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_hnsc_tcga_carboplatin_most_sensitive_1se, hnsc_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_hnsc_tcga_carboplatin_least_sensitive_min, hnsc_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_hnsc_tcga_carboplatin_least_sensitive_1se, hnsc_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/hnsc_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# HNSC W FLUOROURACIL (5)
hnsc_clinical <- read.csv('Processed_Clinical_Data/hnsc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(hnsc_clinical$most_sensitive)
hnsc_clinical <- hnsc_clinical[!na_idx, ]
table(hnsc_clinical$drug_name)
hnsc_clinical_fluorouracil <- hnsc_clinical[which(hnsc_clinical$drug_name == '5-Fluorouracil' | hnsc_clinical$drug_name == '5-FU' | 
                                                    hnsc_clinical$drug_name == '5FU' | hnsc_clinical$drug_name == 'Fluorouracil'), ]

hnsc_clinical_fluorouracil$most_sensitive  <- ifelse(hnsc_clinical_fluorouracil$PFS < quantile(hnsc_clinical_fluorouracil$PFS, probs = .20), 1, 0)
hnsc_clinical_fluorouracil$least_sensitive <- ifelse(hnsc_clinical_fluorouracil$PFS > quantile(hnsc_clinical_fluorouracil$PFS, probs = .80), 1, 0)

hnsc_gene <- read.csv('Processed_Gene_Expression/hnsc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(hnsc_gene) <- gsub('\\.', '-', colnames(hnsc_gene))
hnsc_matching_idx <- hnsc_clinical_fluorouracil$submitter_id.samples %in% colnames(hnsc_gene)
hnsc_clinical_fluorouracil_short <- hnsc_clinical_fluorouracil[hnsc_matching_idx, ]
hnsc_matching_idx <- colnames(hnsc_gene) %in% hnsc_clinical_fluorouracil_short$submitter_id.samples
hnsc_gene_short <- hnsc_gene[, hnsc_matching_idx]
hnsc_gene_short <- t(hnsc_gene_short)
hnsc_gene_short_scaled <- apply(hnsc_gene_short, 2, scale)

rm(hnsc_tcga_most_min_auc)
rm(hnsc_tcga_most_1se_auc)
rm(hnsc_tcga_least_min_auc)
rm(hnsc_tcga_least_1se_auc)

new_hnsc_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_min_auc <- auc(hnsc_clinical_fluorouracil_short$most_sensitive, new_hnsc_tcga_fluorouracil_most_sensitive_min)
hnsc_tcga_most_min_auc <- round(hnsc_tcga_most_min_auc, digits = 2)

new_hnsc_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_1se_auc <- auc(hnsc_clinical_fluorouracil_short$most_sensitive, new_hnsc_tcga_fluorouracil_most_sensitive_1se)
hnsc_tcga_most_1se_auc <- round(hnsc_tcga_most_1se_auc, digits = 2)

new_hnsc_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_min_auc <- auc(hnsc_clinical_fluorouracil_short$least_sensitive, new_hnsc_tcga_fluorouracil_least_sensitive_min)
hnsc_tcga_least_min_auc <- round(hnsc_tcga_least_min_auc, digits = 2)

new_hnsc_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_1se_auc <- auc(hnsc_clinical_fluorouracil_short$least_sensitive, new_hnsc_tcga_fluorouracil_least_sensitive_1se)
hnsc_tcga_least_1se_auc <- round(hnsc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_hnsc_tcga_fluorouracil_most_sensitive_min, hnsc_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_hnsc_tcga_fluorouracil_most_sensitive_1se, hnsc_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_hnsc_tcga_fluorouracil_least_sensitive_min, hnsc_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_hnsc_tcga_fluorouracil_least_sensitive_1se, hnsc_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/hnsc_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# HNSC W PACLITAXEL (18)
hnsc_clinical <- read.csv('Processed_Clinical_Data/hnsc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(hnsc_clinical$most_sensitive)
hnsc_clinical <- hnsc_clinical[!na_idx, ]
table(hnsc_clinical$drug_name)
hnsc_clinical_paclitaxel <- hnsc_clinical[which(hnsc_clinical$drug_name == 'paclitaxel' | hnsc_clinical$drug_name == 'Paclitaxel' | 
                                                  hnsc_clinical$drug_name == 'Palixtaxel'), ]

hnsc_clinical_paclitaxel$most_sensitive  <- ifelse(hnsc_clinical_paclitaxel$PFS < quantile(hnsc_clinical_paclitaxel$PFS, probs = .20), 1, 0)
hnsc_clinical_paclitaxel$least_sensitive <- ifelse(hnsc_clinical_paclitaxel$PFS > quantile(hnsc_clinical_paclitaxel$PFS, probs = .80), 1, 0)

hnsc_gene <- read.csv('Processed_Gene_Expression/hnsc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(hnsc_gene) <- gsub('\\.', '-', colnames(hnsc_gene))
hnsc_matching_idx <- hnsc_clinical_paclitaxel$submitter_id.samples %in% colnames(hnsc_gene)
hnsc_clinical_paclitaxel_short <- hnsc_clinical_paclitaxel[hnsc_matching_idx, ]
hnsc_matching_idx <- colnames(hnsc_gene) %in% hnsc_clinical_paclitaxel_short$submitter_id.samples
hnsc_gene_short <- hnsc_gene[, hnsc_matching_idx]
hnsc_gene_short <- t(hnsc_gene_short)
hnsc_gene_short_scaled <- apply(hnsc_gene_short, 2, scale)

rm(hnsc_tcga_most_min_auc)
rm(hnsc_tcga_most_1se_auc)
rm(hnsc_tcga_least_min_auc)
rm(hnsc_tcga_least_1se_auc)

new_hnsc_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_min_auc <- auc(hnsc_clinical_paclitaxel_short$most_sensitive, new_hnsc_tcga_paclitaxel_most_sensitive_min)
hnsc_tcga_most_min_auc <- round(hnsc_tcga_most_min_auc, digits = 2)

new_hnsc_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_most_1se_auc <- auc(hnsc_clinical_paclitaxel_short$most_sensitive, new_hnsc_tcga_paclitaxel_most_sensitive_1se)
hnsc_tcga_most_1se_auc <- round(hnsc_tcga_most_1se_auc, digits = 2)

new_hnsc_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_min_auc <- auc(hnsc_clinical_paclitaxel_short$least_sensitive, new_hnsc_tcga_paclitaxel_least_sensitive_min)
hnsc_tcga_least_min_auc <- round(hnsc_tcga_least_min_auc, digits = 2)

new_hnsc_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = hnsc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
hnsc_tcga_least_1se_auc <- auc(hnsc_clinical_paclitaxel_short$least_sensitive, new_hnsc_tcga_paclitaxel_least_sensitive_1se)
hnsc_tcga_least_1se_auc <- round(hnsc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_hnsc_tcga_paclitaxel_most_sensitive_min, hnsc_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_hnsc_tcga_paclitaxel_most_sensitive_1se, hnsc_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_hnsc_tcga_paclitaxel_least_sensitive_min, hnsc_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_hnsc_tcga_paclitaxel_least_sensitive_1se, hnsc_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/hnsc_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LIHC W GEMCITABINE (8)
lihc_clinical <- read.csv('Processed_Clinical_Data/lihc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lihc_clinical$most_sensitive)
lihc_clinical <- lihc_clinical[!na_idx, ]
table(lihc_clinical$drug_name)
lihc_clinical_gemcitabine <- lihc_clinical[which(lihc_clinical$drug_name == 'Gemcitabine'), ]

lihc_clinical_gemcitabine$most_sensitive  <- ifelse(lihc_clinical_gemcitabine$PFS < quantile(lihc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
lihc_clinical_gemcitabine$least_sensitive <- ifelse(lihc_clinical_gemcitabine$PFS > quantile(lihc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

lihc_gene <- read.csv('Processed_Gene_Expression/lihc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lihc_gene) <- gsub('\\.', '-', colnames(lihc_gene))
lihc_matching_idx <- lihc_clinical_gemcitabine$submitter_id.samples %in% colnames(lihc_gene)
lihc_clinical_gemcitabine_short <- lihc_clinical_gemcitabine[lihc_matching_idx, ]
lihc_matching_idx <- colnames(lihc_gene) %in% lihc_clinical_gemcitabine_short$submitter_id.samples
lihc_gene_short <- lihc_gene[, lihc_matching_idx]
lihc_gene_short <- t(lihc_gene_short)
lihc_gene_short_scaled <- apply(lihc_gene_short, 2, scale)

rm(lihc_tcga_most_min_auc)
rm(lihc_tcga_most_1se_auc)
rm(lihc_tcga_least_min_auc)
rm(lihc_tcga_least_1se_auc)

new_lihc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_most_min_auc <- auc(lihc_clinical_gemcitabine_short$most_sensitive, new_lihc_tcga_gemcitabine_most_sensitive_min)
lihc_tcga_most_min_auc <- round(lihc_tcga_most_min_auc, digits = 2)

new_lihc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_most_1se_auc <- auc(lihc_clinical_gemcitabine_short$most_sensitive, new_lihc_tcga_gemcitabine_most_sensitive_1se)
lihc_tcga_most_1se_auc <- round(lihc_tcga_most_1se_auc, digits = 2)

new_lihc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lihc_tcga_least_min_auc <- auc(lihc_clinical_gemcitabine_short$least_sensitive, new_lihc_tcga_gemcitabine_least_sensitive_min)
lihc_tcga_least_min_auc <- round(lihc_tcga_least_min_auc, digits = 2)

new_lihc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = lihc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lihc_tcga_least_1se_auc <- auc(lihc_clinical_gemcitabine_short$least_sensitive, new_lihc_tcga_gemcitabine_least_sensitive_1se)
lihc_tcga_least_1se_auc <- round(lihc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lihc_tcga_gemcitabine_most_sensitive_min, lihc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_lihc_tcga_gemcitabine_most_sensitive_1se, lihc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_lihc_tcga_gemcitabine_least_sensitive_min, lihc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_lihc_tcga_gemcitabine_least_sensitive_1se, lihc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lihc_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[7], lwd = 2, lty = 2, main = 'LIHC treated with gemcitabine')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[7], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.35, legend = c('sensitive model (AUC = 0.58)', 'resistant model (AUC = 0.58)'), cex = 0.8, lty = c(1,2), lwd = 2, col = colors_i_need[7], bty = 'n')
dev.off()

# LUAD W CARBOPLATIN (26)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_carboplatin <- luad_clinical[which(luad_clinical$drug_name == 'carboplatin' | luad_clinical$drug_name == 'Carboplatin' | 
                                                   luad_clinical$drug_name == 'CARBOPLATIN'), ]

luad_clinical_carboplatin$most_sensitive  <- ifelse(luad_clinical_carboplatin$PFS < quantile(luad_clinical_carboplatin$PFS, probs = .20), 1, 0)
luad_clinical_carboplatin$least_sensitive <- ifelse(luad_clinical_carboplatin$PFS > quantile(luad_clinical_carboplatin$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_carboplatin$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_carboplatin_short <- luad_clinical_carboplatin[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_carboplatin_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_carboplatin_short$most_sensitive, new_luad_tcga_carboplatin_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_carboplatin_short$most_sensitive, new_luad_tcga_carboplatin_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_carboplatin_short$least_sensitive, new_luad_tcga_carboplatin_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_carboplatin_short$least_sensitive, new_luad_tcga_carboplatin_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_carboplatin_most_sensitive_min, luad_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_carboplatin_most_sensitive_1se, luad_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_carboplatin_least_sensitive_min, luad_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_carboplatin_least_sensitive_1se, luad_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUAD W DOCETAXEL (9)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_docetaxel <- luad_clinical[which(luad_clinical$drug_name == 'Docetaxel' | luad_clinical$drug_name == 'DOCETAXEL'), ]

luad_clinical_docetaxel$most_sensitive  <- ifelse(luad_clinical_docetaxel$PFS < quantile(luad_clinical_docetaxel$PFS, probs = .20), 1, 0)
luad_clinical_docetaxel$least_sensitive <- ifelse(luad_clinical_docetaxel$PFS > quantile(luad_clinical_docetaxel$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_docetaxel$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_docetaxel_short <- luad_clinical_docetaxel[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_docetaxel_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_docetaxel_short$most_sensitive, new_luad_tcga_docetaxel_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_docetaxel_short$most_sensitive, new_luad_tcga_docetaxel_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_docetaxel_short$least_sensitive, new_luad_tcga_docetaxel_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_docetaxel_short$least_sensitive, new_luad_tcga_docetaxel_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_docetaxel_most_sensitive_min, luad_clinical_docetaxel_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_docetaxel_most_sensitive_1se, luad_clinical_docetaxel_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_docetaxel_least_sensitive_min, luad_clinical_docetaxel_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_docetaxel_least_sensitive_1se, luad_clinical_docetaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_docetaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUAD W ETOPOSIDE (7)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_etoposide <- luad_clinical[which(luad_clinical$drug_name == 'Etoposide'), ]

luad_clinical_etoposide$most_sensitive  <- ifelse(luad_clinical_etoposide$PFS < quantile(luad_clinical_etoposide$PFS, probs = .20), 1, 0)
luad_clinical_etoposide$least_sensitive <- ifelse(luad_clinical_etoposide$PFS > quantile(luad_clinical_etoposide$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_etoposide$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_etoposide_short <- luad_clinical_etoposide[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_etoposide_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_etoposide_most_sensitive_min <- predict(etoposide_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_etoposide_short$most_sensitive, new_luad_tcga_etoposide_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_etoposide_most_sensitive_1se <- predict(etoposide_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_etoposide_short$most_sensitive, new_luad_tcga_etoposide_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_etoposide_least_sensitive_min <- predict(etoposide_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_etoposide_short$least_sensitive, new_luad_tcga_etoposide_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_etoposide_least_sensitive_1se <- predict(etoposide_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_etoposide_short$least_sensitive, new_luad_tcga_etoposide_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_etoposide_most_sensitive_min, luad_clinical_etoposide_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_etoposide_most_sensitive_1se, luad_clinical_etoposide_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_etoposide_least_sensitive_min, luad_clinical_etoposide_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_etoposide_least_sensitive_1se, luad_clinical_etoposide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_etoposide_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUAD W GEMCITABINE (18)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_gemcitabine <- luad_clinical[which(luad_clinical$drug_name == 'gemcitabine' | luad_clinical$drug_name == 'Gemcitabine' | 
                                                   luad_clinical$drug_name == 'Gemzar'), ]

luad_clinical_gemcitabine$most_sensitive  <- ifelse(luad_clinical_gemcitabine$PFS < quantile(luad_clinical_gemcitabine$PFS, probs = .20), 1, 0)
luad_clinical_gemcitabine$least_sensitive <- ifelse(luad_clinical_gemcitabine$PFS > quantile(luad_clinical_gemcitabine$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_gemcitabine$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_gemcitabine_short <- luad_clinical_gemcitabine[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_gemcitabine_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_gemcitabine_short$most_sensitive, new_luad_tcga_gemcitabine_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_gemcitabine_short$most_sensitive, new_luad_tcga_gemcitabine_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_gemcitabine_short$least_sensitive, new_luad_tcga_gemcitabine_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_gemcitabine_short$least_sensitive, new_luad_tcga_gemcitabine_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_gemcitabine_most_sensitive_min, luad_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_gemcitabine_most_sensitive_1se, luad_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_gemcitabine_least_sensitive_min, luad_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_gemcitabine_least_sensitive_1se, luad_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUAD W PACLITAXEL (8)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_paclitaxel <- luad_clinical[which(luad_clinical$drug_name == 'paclitaxel' | luad_clinical$drug_name == 'Paclitaxel'), ]

luad_clinical_paclitaxel$most_sensitive  <- ifelse(luad_clinical_paclitaxel$PFS < quantile(luad_clinical_paclitaxel$PFS, probs = .20), 1, 0)
luad_clinical_paclitaxel$least_sensitive <- ifelse(luad_clinical_paclitaxel$PFS > quantile(luad_clinical_paclitaxel$PFS, probs = .80), 1, 0)

luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_paclitaxel$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_paclitaxel_short <- luad_clinical_paclitaxel[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_paclitaxel_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)

rm(luad_tcga_most_min_auc)
rm(luad_tcga_most_1se_auc)
rm(luad_tcga_least_min_auc)
rm(luad_tcga_least_1se_auc)

new_luad_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_most_min_auc <- auc(luad_clinical_paclitaxel_short$most_sensitive, new_luad_tcga_paclitaxel_most_sensitive_min)
luad_tcga_most_min_auc <- round(luad_tcga_most_min_auc, digits = 2)

new_luad_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_most_1se_auc <- auc(luad_clinical_paclitaxel_short$most_sensitive, new_luad_tcga_paclitaxel_most_sensitive_1se)
luad_tcga_most_1se_auc <- round(luad_tcga_most_1se_auc, digits = 2)

new_luad_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
luad_tcga_least_min_auc <- auc(luad_clinical_paclitaxel_short$least_sensitive, new_luad_tcga_paclitaxel_least_sensitive_min)
luad_tcga_least_min_auc <- round(luad_tcga_least_min_auc, digits = 2)

new_luad_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = luad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
luad_tcga_least_1se_auc <- auc(luad_clinical_paclitaxel_short$least_sensitive, new_luad_tcga_paclitaxel_least_sensitive_1se)
luad_tcga_least_1se_auc <- round(luad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_luad_tcga_paclitaxel_most_sensitive_min, luad_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_luad_tcga_paclitaxel_most_sensitive_1se, luad_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_luad_tcga_paclitaxel_least_sensitive_min, luad_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_luad_tcga_paclitaxel_least_sensitive_1se, luad_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/luad_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUSC W CARBOPLATIN (28)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_carboplatin <- lusc_clinical[which(lusc_clinical$drug_name == 'Carboplatin' | lusc_clinical$drug_name == 'carboplatin'), ]

lusc_clinical_carboplatin$most_sensitive  <- ifelse(lusc_clinical_carboplatin$PFS < quantile(lusc_clinical_carboplatin$PFS, probs = .20), 1, 0)
lusc_clinical_carboplatin$least_sensitive <- ifelse(lusc_clinical_carboplatin$PFS > quantile(lusc_clinical_carboplatin$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_carboplatin$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_carboplatin_short <- lusc_clinical_carboplatin[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_carboplatin_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_carboplatin_short$most_sensitive, new_lusc_tcga_carboplatin_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_carboplatin_short$most_sensitive, new_lusc_tcga_carboplatin_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_carboplatin_short$least_sensitive, new_lusc_tcga_carboplatin_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_carboplatin_short$least_sensitive, new_lusc_tcga_carboplatin_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_carboplatin_most_sensitive_min, lusc_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_carboplatin_most_sensitive_1se, lusc_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_carboplatin_least_sensitive_min, lusc_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_carboplatin_least_sensitive_1se, lusc_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[1], lwd = 2, lty = 2, main = 'LUSC treated with carboplatin')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[1], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.50)', 'resistant model (AUC = 0.85)'), cex = 0.8, lwd = 2, lty = c(1,2), col = colors_i_need[1], bty = 'n')
dev.off()

# LUSC W DOCETAXEL (6)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_docetaxel <- lusc_clinical[which(lusc_clinical$drug_name == 'Docetaxel'), ]

lusc_clinical_docetaxel$most_sensitive  <- ifelse(lusc_clinical_docetaxel$PFS < quantile(lusc_clinical_docetaxel$PFS, probs = .20), 1, 0)
lusc_clinical_docetaxel$least_sensitive <- ifelse(lusc_clinical_docetaxel$PFS > quantile(lusc_clinical_docetaxel$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_docetaxel$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_docetaxel_short <- lusc_clinical_docetaxel[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_docetaxel_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_docetaxel_short$most_sensitive, new_lusc_tcga_docetaxel_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_docetaxel_short$most_sensitive, new_lusc_tcga_docetaxel_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_docetaxel_short$least_sensitive, new_lusc_tcga_docetaxel_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_docetaxel_short$least_sensitive, new_lusc_tcga_docetaxel_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_docetaxel_most_sensitive_min, lusc_clinical_docetaxel_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_docetaxel_most_sensitive_1se, lusc_clinical_docetaxel_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_docetaxel_least_sensitive_min, lusc_clinical_docetaxel_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_docetaxel_least_sensitive_1se, lusc_clinical_docetaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_docetaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# LUSC W GEMCITABINE (30)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_gemcitabine <- lusc_clinical[which(lusc_clinical$drug_name == 'gemcitabine' | lusc_clinical$drug_name == 'Gemcitabine' | 
                                                   lusc_clinical$drug_name == 'Gemzar' | lusc_clinical$drug_name == 'GEMZAR' | 
                                                   lusc_clinical$drug_name == 'Gemzar (Gemcitabine)'), ]

lusc_clinical_gemcitabine$most_sensitive  <- ifelse(lusc_clinical_gemcitabine$PFS < quantile(lusc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
lusc_clinical_gemcitabine$least_sensitive <- ifelse(lusc_clinical_gemcitabine$PFS > quantile(lusc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_gemcitabine$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_gemcitabine_short <- lusc_clinical_gemcitabine[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_gemcitabine_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)

rm(lusc_tcga_most_min_auc)
rm(lusc_tcga_most_1se_auc)
rm(lusc_tcga_least_min_auc)
rm(lusc_tcga_least_1se_auc)

new_lusc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_most_min_auc <- auc(lusc_clinical_gemcitabine_short$most_sensitive, new_lusc_tcga_gemcitabine_most_sensitive_min)
lusc_tcga_most_min_auc <- round(lusc_tcga_most_min_auc, digits = 2)

new_lusc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_most_1se_auc <- auc(lusc_clinical_gemcitabine_short$most_sensitive, new_lusc_tcga_gemcitabine_most_sensitive_1se)
lusc_tcga_most_1se_auc <- round(lusc_tcga_most_1se_auc, digits = 2)

new_lusc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
lusc_tcga_least_min_auc <- auc(lusc_clinical_gemcitabine_short$least_sensitive, new_lusc_tcga_gemcitabine_least_sensitive_min)
lusc_tcga_least_min_auc <- round(lusc_tcga_least_min_auc, digits = 2)

new_lusc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = lusc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
lusc_tcga_least_1se_auc <- auc(lusc_clinical_gemcitabine_short$least_sensitive, new_lusc_tcga_gemcitabine_least_sensitive_1se)
lusc_tcga_least_1se_auc <- round(lusc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_lusc_tcga_gemcitabine_most_sensitive_min, lusc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_lusc_tcga_gemcitabine_most_sensitive_1se, lusc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_lusc_tcga_gemcitabine_least_sensitive_min, lusc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_lusc_tcga_gemcitabine_least_sensitive_1se, lusc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/lusc_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W CARBOPLATIN (231)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_carboplatin <- ov_clinical[which(ov_clinical$drug_name == 'Carboplatin' | ov_clinical$drug_name == 'Carbo' | 
                                               ov_clinical$drug_name == 'carboplatin' | ov_clinical$drug_name == 'Carbplatin'), ]

ov_clinical_carboplatin$most_sensitive  <- ifelse(ov_clinical_carboplatin$PFS < quantile(ov_clinical_carboplatin$PFS, probs = .20), 1, 0)
ov_clinical_carboplatin$least_sensitive <- ifelse(ov_clinical_carboplatin$PFS > quantile(ov_clinical_carboplatin$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_carboplatin$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_carboplatin_short <- ov_clinical_carboplatin[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_carboplatin_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_carboplatin_short$most_sensitive, new_ov_tcga_carboplatin_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_carboplatin_short$most_sensitive, new_ov_tcga_carboplatin_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_carboplatin_short$least_sensitive, new_ov_tcga_carboplatin_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_carboplatin_short$least_sensitive, new_ov_tcga_carboplatin_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_carboplatin_most_sensitive_min, ov_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_carboplatin_most_sensitive_1se, ov_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_carboplatin_least_sensitive_min, ov_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_carboplatin_least_sensitive_1se, ov_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W DOCETAXEL (15)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_docetaxel <- ov_clinical[which(ov_clinical$drug_name == 'Docetaxel' | ov_clinical$drug_name == 'Doxetaxel'), ]

ov_clinical_docetaxel$most_sensitive  <- ifelse(ov_clinical_docetaxel$PFS < quantile(ov_clinical_docetaxel$PFS, probs = .20), 1, 0)
ov_clinical_docetaxel$least_sensitive <- ifelse(ov_clinical_docetaxel$PFS > quantile(ov_clinical_docetaxel$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_docetaxel$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_docetaxel_short <- ov_clinical_docetaxel[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_docetaxel_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_docetaxel_short$most_sensitive, new_ov_tcga_docetaxel_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_docetaxel_short$most_sensitive, new_ov_tcga_docetaxel_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_docetaxel_short$least_sensitive, new_ov_tcga_docetaxel_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_docetaxel_short$least_sensitive, new_ov_tcga_docetaxel_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_docetaxel_most_sensitive_min, ov_clinical_docetaxel_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_docetaxel_most_sensitive_1se, ov_clinical_docetaxel_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_docetaxel_least_sensitive_min, ov_clinical_docetaxel_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_docetaxel_least_sensitive_1se, ov_clinical_docetaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_docetaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W DOXORUBICIN (5)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_doxorubicin <- ov_clinical[which(ov_clinical$drug_name == 'Doxorubicin' | ov_clinical$drug_name == 'Doxorubicin HCL'), ]

ov_clinical_doxorubicin$most_sensitive  <- ifelse(ov_clinical_doxorubicin$PFS < quantile(ov_clinical_doxorubicin$PFS, probs = .20), 1, 0)
ov_clinical_doxorubicin$least_sensitive <- ifelse(ov_clinical_doxorubicin$PFS > quantile(ov_clinical_doxorubicin$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_doxorubicin$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_doxorubicin_short <- ov_clinical_doxorubicin[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_doxorubicin_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_doxorubicin_most_sensitive_min <- predict(doxorubicin_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_doxorubicin_short$most_sensitive, new_ov_tcga_doxorubicin_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_doxorubicin_most_sensitive_1se <- predict(doxorubicin_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_doxorubicin_short$most_sensitive, new_ov_tcga_doxorubicin_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_doxorubicin_least_sensitive_min <- predict(doxorubicin_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_doxorubicin_short$least_sensitive, new_ov_tcga_doxorubicin_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_doxorubicin_least_sensitive_1se <- predict(doxorubicin_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_doxorubicin_short$least_sensitive, new_ov_tcga_doxorubicin_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_doxorubicin_most_sensitive_min, ov_clinical_doxorubicin_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_doxorubicin_most_sensitive_1se, ov_clinical_doxorubicin_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_doxorubicin_least_sensitive_min, ov_clinical_doxorubicin_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_doxorubicin_least_sensitive_1se, ov_clinical_doxorubicin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_doxorubicin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W GEMCITABINE (41)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_gemcitabine <- ov_clinical[which(ov_clinical$drug_name == 'gemcitabine' | ov_clinical$drug_name == 'Gamzar' | 
                                               ov_clinical$drug_name == 'Gemcitabine' | ov_clinical$drug_name == 'Gemcitabine HCL' | 
                                               ov_clinical$drug_name == 'Gemcitibine' | ov_clinical$drug_name = 'Gemzar'), ]

ov_clinical_gemcitabine$most_sensitive  <- ifelse(ov_clinical_gemcitabine$PFS < quantile(ov_clinical_gemcitabine$PFS, probs = .20), 1, 0)
ov_clinical_gemcitabine$least_sensitive <- ifelse(ov_clinical_gemcitabine$PFS > quantile(ov_clinical_gemcitabine$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_gemcitabine$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_gemcitabine_short <- ov_clinical_gemcitabine[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_gemcitabine_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_gemcitabine_short$most_sensitive, new_ov_tcga_gemcitabine_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_gemcitabine_short$most_sensitive, new_ov_tcga_gemcitabine_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_gemcitabine_short$least_sensitive, new_ov_tcga_gemcitabine_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_gemcitabine_short$least_sensitive, new_ov_tcga_gemcitabine_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_gemcitabine_most_sensitive_min, ov_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_gemcitabine_most_sensitive_1se, ov_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_gemcitabine_least_sensitive_min, ov_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_gemcitabine_least_sensitive_1se, ov_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W PACLITAXEL (75)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_paclitaxel <- ov_clinical[which(ov_clinical$drug_name == 'Pacilatxel' | ov_clinical$drug_name == 'Paciltaxel' | 
                                              ov_clinical$drug_name == 'Paciltaxle' | ov_clinical$drug_name == 'Pacliltaxel' | 
                                              ov_clinical$drug_name == 'Paclitaxel' | ov_clinical$drug_name == 'Paclitaxel; Albumin-Bount'), ]

ov_clinical_paclitaxel$most_sensitive  <- ifelse(ov_clinical_paclitaxel$PFS < quantile(ov_clinical_paclitaxel$PFS, probs = .20), 1, 0)
ov_clinical_paclitaxel$least_sensitive <- ifelse(ov_clinical_paclitaxel$PFS > quantile(ov_clinical_paclitaxel$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_paclitaxel$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_paclitaxel_short <- ov_clinical_paclitaxel[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_paclitaxel_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_paclitaxel_short$most_sensitive, new_ov_tcga_paclitaxel_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_paclitaxel_short$most_sensitive, new_ov_tcga_paclitaxel_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_paclitaxel_short$least_sensitive, new_ov_tcga_paclitaxel_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_paclitaxel_short$least_sensitive, new_ov_tcga_paclitaxel_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_paclitaxel_most_sensitive_min, ov_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_paclitaxel_most_sensitive_1se, ov_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_paclitaxel_least_sensitive_min, ov_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_paclitaxel_least_sensitive_1se, ov_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# OV W TOPOTECAN (33)
ov_clinical <- read.csv('Processed_Clinical_Data/ov_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ov_clinical$most_sensitive)
ov_clinical <- ov_clinical[!na_idx, ]
table(ov_clinical$drug_name)
ov_clinical_topotecan <- ov_clinical[which(ov_clinical$drug_name == 'Topotecan' | ov_clinical$drug_name == 'Topotecan HCL' | 
                                             ov_clinical$drug_name == 'Toptecan'), ]

ov_clinical_topotecan$most_sensitive  <- ifelse(ov_clinical_topotecan$PFS < quantile(ov_clinical_topotecan$PFS, probs = .20), 1, 0)
ov_clinical_topotecan$least_sensitive <- ifelse(ov_clinical_topotecan$PFS > quantile(ov_clinical_topotecan$PFS, probs = .80), 1, 0)

ov_gene <- read.csv('Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ov_gene) <- gsub('\\.', '-', colnames(ov_gene))
ov_matching_idx <- ov_clinical_topotecan$submitter_id.samples %in% colnames(ov_gene)
ov_clinical_topotecan_short <- ov_clinical_topotecan[ov_matching_idx, ]
ov_matching_idx <- colnames(ov_gene) %in% ov_clinical_topotecan_short$submitter_id.samples
ov_gene_short <- ov_gene[, ov_matching_idx]
ov_gene_short <- t(ov_gene_short)
ov_gene_short_scaled <- apply(ov_gene_short, 2, scale)

rm(ov_tcga_most_min_auc)
rm(ov_tcga_most_1se_auc)
rm(ov_tcga_least_min_auc)
rm(ov_tcga_least_1se_auc)

new_ov_tcga_topotecan_most_sensitive_min <- predict(topotecan_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_most_min_auc <- auc(ov_clinical_topotecan_short$most_sensitive, new_ov_tcga_topotecan_most_sensitive_min)
ov_tcga_most_min_auc <- round(ov_tcga_most_min_auc, digits = 2)

new_ov_tcga_topotecan_most_sensitive_1se <- predict(topotecan_ccle_most_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_most_1se_auc <- auc(ov_clinical_topotecan_short$most_sensitive, new_ov_tcga_topotecan_most_sensitive_1se)
ov_tcga_most_1se_auc <- round(ov_tcga_most_1se_auc, digits = 2)

new_ov_tcga_topotecan_least_sensitive_min <- predict(topotecan_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ov_tcga_least_min_auc <- auc(ov_clinical_topotecan_short$least_sensitive, new_ov_tcga_topotecan_least_sensitive_min)
ov_tcga_least_min_auc <- round(ov_tcga_least_min_auc, digits = 2)

new_ov_tcga_topotecan_least_sensitive_1se <- predict(topotecan_ccle_least_fit_elnet, newx = ov_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ov_tcga_least_1se_auc <- auc(ov_clinical_topotecan_short$least_sensitive, new_ov_tcga_topotecan_least_sensitive_1se)
ov_tcga_least_1se_auc <- round(ov_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ov_tcga_topotecan_most_sensitive_min, ov_clinical_topotecan_short$most_sensitive)
pred2 <- prediction(new_ov_tcga_topotecan_most_sensitive_1se, ov_clinical_topotecan_short$most_sensitive)
pred3 <- prediction(new_ov_tcga_topotecan_least_sensitive_min, ov_clinical_topotecan_short$least_sensitive)
pred4 <- prediction(new_ov_tcga_topotecan_least_sensitive_1se, ov_clinical_topotecan_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ov_topotecan_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# PAAD W FLUOROURACIL (14)
paad_clinical <- read.csv('Processed_Clinical_Data/paad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(paad_clinical$most_sensitive)
paad_clinical <- paad_clinical[!na_idx, ]
table(paad_clinical$drug_name)
paad_clinical_fluorouracil <- paad_clinical[which(paad_clinical$drug_name == '5-fluorouracil' | paad_clinical$drug_name == '5-Fluorouracil' | 
                                                    paad_clinical$drug_name == '5-fu' | paad_clinical$drug_name == '5-FU' | paad_clinical$drug_name == '5 FU' | 
                                                    paad_clinical$drug_name == '5FU' | paad_clinical$drug_name == 'Fluorouracil'), ]

paad_clinical_fluorouracil$most_sensitive  <- ifelse(paad_clinical_fluorouracil$PFS < quantile(paad_clinical_fluorouracil$PFS, probs = .20), 1, 0)
paad_clinical_fluorouracil$least_sensitive <- ifelse(paad_clinical_fluorouracil$PFS > quantile(paad_clinical_fluorouracil$PFS, probs = .80), 1, 0)

paad_gene <- read.csv('Processed_Gene_Expression/paad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(paad_gene) <- gsub('\\.', '-', colnames(paad_gene))
paad_matching_idx <- paad_clinical_fluorouracil$submitter_id.samples %in% colnames(paad_gene)
paad_clinical_fluorouracil_short <- paad_clinical_fluorouracil[paad_matching_idx, ]
paad_matching_idx <- colnames(paad_gene) %in% paad_clinical_fluorouracil_short$submitter_id.samples
paad_gene_short <- paad_gene[, paad_matching_idx]
paad_gene_short <- t(paad_gene_short)
paad_gene_short_scaled <- apply(paad_gene_short, 2, scale)

rm(paad_tcga_most_min_auc)
rm(paad_tcga_most_1se_auc)
rm(paad_tcga_least_min_auc)
rm(paad_tcga_least_1se_auc)

new_paad_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_most_min_auc <- auc(paad_clinical_fluorouracil_short$most_sensitive, new_paad_tcga_fluorouracil_most_sensitive_min)
paad_tcga_most_min_auc <- round(paad_tcga_most_min_auc, digits = 2)

new_paad_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_most_1se_auc <- auc(paad_clinical_fluorouracil_short$most_sensitive, new_paad_tcga_fluorouracil_most_sensitive_1se)
paad_tcga_most_1se_auc <- round(paad_tcga_most_1se_auc, digits = 2)

new_paad_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_least_min_auc <- auc(paad_clinical_fluorouracil_short$least_sensitive, new_paad_tcga_fluorouracil_least_sensitive_min)
paad_tcga_least_min_auc <- round(paad_tcga_least_min_auc, digits = 2)

new_paad_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_least_1se_auc <- auc(paad_clinical_fluorouracil_short$least_sensitive, new_paad_tcga_fluorouracil_least_sensitive_1se)
paad_tcga_least_1se_auc <- round(paad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_paad_tcga_fluorouracil_most_sensitive_min, paad_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_paad_tcga_fluorouracil_most_sensitive_1se, paad_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_paad_tcga_fluorouracil_least_sensitive_min, paad_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_paad_tcga_fluorouracil_least_sensitive_1se, paad_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/paad_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# PAAD W GEMCITABINE (39)
paad_clinical <- read.csv('Processed_Clinical_Data/paad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(paad_clinical$most_sensitive)
paad_clinical <- paad_clinical[!na_idx, ]
table(paad_clinical$drug_name)
paad_clinical_gemcitabine <- paad_clinical[which(paad_clinical$drug_name == 'Gemcitabine' | paad_clinical$drug_name == 'Gemcitabine Injection' | 
                                                   paad_clinical$drug_name == 'gemzar' | paad_clinical$drug_name == 'Gemzar'), ]

paad_clinical_gemcitabine$most_sensitive  <- ifelse(paad_clinical_gemcitabine$PFS < quantile(paad_clinical_gemcitabine$PFS, probs = .20), 1, 0)
paad_clinical_gemcitabine$least_sensitive <- ifelse(paad_clinical_gemcitabine$PFS > quantile(paad_clinical_gemcitabine$PFS, probs = .80), 1, 0)

paad_gene <- read.csv('Processed_Gene_Expression/paad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(paad_gene) <- gsub('\\.', '-', colnames(paad_gene))
paad_matching_idx <- paad_clinical_gemcitabine$submitter_id.samples %in% colnames(paad_gene)
paad_clinical_gemcitabine_short <- paad_clinical_gemcitabine[paad_matching_idx, ]
paad_matching_idx <- colnames(paad_gene) %in% paad_clinical_gemcitabine_short$submitter_id.samples
paad_gene_short <- paad_gene[, paad_matching_idx]
paad_gene_short <- t(paad_gene_short)
paad_gene_short_scaled <- apply(paad_gene_short, 2, scale)

rm(paad_tcga_most_min_auc)
rm(paad_tcga_most_1se_auc)
rm(paad_tcga_least_min_auc)
rm(paad_tcga_least_1se_auc)

new_paad_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_most_min_auc <- auc(paad_clinical_gemcitabine_short$most_sensitive, new_paad_tcga_gemcitabine_most_sensitive_min)
paad_tcga_most_min_auc <- round(paad_tcga_most_min_auc, digits = 2)

new_paad_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_most_1se_auc <- auc(paad_clinical_gemcitabine_short$most_sensitive, new_paad_tcga_gemcitabine_most_sensitive_1se)
paad_tcga_most_1se_auc <- round(paad_tcga_most_1se_auc, digits = 2)

new_paad_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_least_min_auc <- auc(paad_clinical_gemcitabine_short$least_sensitive, new_paad_tcga_gemcitabine_least_sensitive_min)
paad_tcga_least_min_auc <- round(paad_tcga_least_min_auc, digits = 2)

new_paad_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_least_1se_auc <- auc(paad_clinical_gemcitabine_short$least_sensitive, new_paad_tcga_gemcitabine_least_sensitive_1se)
paad_tcga_least_1se_auc <- round(paad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_paad_tcga_gemcitabine_most_sensitive_min, paad_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_paad_tcga_gemcitabine_most_sensitive_1se, paad_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_paad_tcga_gemcitabine_least_sensitive_min, paad_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_paad_tcga_gemcitabine_least_sensitive_1se, paad_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/paad_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[7], lwd = 2, lty = 2, main = 'PAAD treated with gemcitabine')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[7], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.56)', 'resistant model (AUC = 0.61)'), cex = 0.8, lwd = 2, lty = c(1,2), col = colors_i_need[7], bty = 'n')
dev.off()

# PAAD W OXALAPLATIN (11)
paad_clinical <- read.csv('Processed_Clinical_Data/paad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(paad_clinical$most_sensitive)
paad_clinical <- paad_clinical[!na_idx, ]
table(paad_clinical$drug_name)
paad_clinical_oxalaplatin <- paad_clinical[which(paad_clinical$drug_name == 'oxaliplatin' | paad_clinical$drug_name == 'Oxaliplatin'), ]

paad_clinical_oxalaplatin$most_sensitive  <- ifelse(paad_clinical_oxalaplatin$PFS < quantile(paad_clinical_oxalaplatin$PFS, probs = .20), 1, 0)
paad_clinical_oxalaplatin$least_sensitive <- ifelse(paad_clinical_oxalaplatin$PFS > quantile(paad_clinical_oxalaplatin$PFS, probs = .80), 1, 0)

paad_gene <- read.csv('Processed_Gene_Expression/paad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(paad_gene) <- gsub('\\.', '-', colnames(paad_gene))
paad_matching_idx <- paad_clinical_oxalaplatin$submitter_id.samples %in% colnames(paad_gene)
paad_clinical_oxalaplatin_short <- paad_clinical_oxalaplatin[paad_matching_idx, ]
paad_matching_idx <- colnames(paad_gene) %in% paad_clinical_oxalaplatin_short$submitter_id.samples
paad_gene_short <- paad_gene[, paad_matching_idx]
paad_gene_short <- t(paad_gene_short)
paad_gene_short_scaled <- apply(paad_gene_short, 2, scale)

rm(paad_tcga_most_min_auc)
rm(paad_tcga_most_1se_auc)
rm(paad_tcga_least_min_auc)
rm(paad_tcga_least_1se_auc)

new_paad_tcga_oxalaplatin_most_sensitive_min <- predict(oxalaplatin_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_most_min_auc <- auc(paad_clinical_oxalaplatin_short$most_sensitive, new_paad_tcga_oxalaplatin_most_sensitive_min)
paad_tcga_most_min_auc <- round(paad_tcga_most_min_auc, digits = 2)

new_paad_tcga_oxalaplatin_most_sensitive_1se <- predict(oxalaplatin_ccle_most_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_most_1se_auc <- auc(paad_clinical_oxalaplatin_short$most_sensitive, new_paad_tcga_oxalaplatin_most_sensitive_1se)
paad_tcga_most_1se_auc <- round(paad_tcga_most_1se_auc, digits = 2)

new_paad_tcga_oxalaplatin_least_sensitive_min <- predict(oxalaplatin_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
paad_tcga_least_min_auc <- auc(paad_clinical_oxalaplatin_short$least_sensitive, new_paad_tcga_oxalaplatin_least_sensitive_min)
paad_tcga_least_min_auc <- round(paad_tcga_least_min_auc, digits = 2)

new_paad_tcga_oxalaplatin_least_sensitive_1se <- predict(oxalaplatin_ccle_least_fit_elnet, newx = paad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
paad_tcga_least_1se_auc <- auc(paad_clinical_oxalaplatin_short$least_sensitive, new_paad_tcga_oxalaplatin_least_sensitive_1se)
paad_tcga_least_1se_auc <- round(paad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_paad_tcga_oxalaplatin_most_sensitive_min, paad_clinical_oxalaplatin_short$most_sensitive)
pred2 <- prediction(new_paad_tcga_oxalaplatin_most_sensitive_1se, paad_clinical_oxalaplatin_short$most_sensitive)
pred3 <- prediction(new_paad_tcga_oxalaplatin_least_sensitive_min, paad_clinical_oxalaplatin_short$least_sensitive)
pred4 <- prediction(new_paad_tcga_oxalaplatin_least_sensitive_1se, paad_clinical_oxalaplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/paad_oxalaplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# READ W FLUOROURACIL (24)
read_clinical <- read.csv('Processed_Clinical_Data/read_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(read_clinical$most_sensitive)
read_clinical <- read_clinical[!na_idx, ]
table(read_clinical$drug_name)
read_clinical_fluorouracil <- read_clinical[which(read_clinical$drug_name == '5-Fluorouracil' | read_clinical$drug_name == '5 FU' | 
                                                    read_clinical$drug_name == '5FU' | read_clinical$drug_name == 'Fluorouracil'), ]

read_clinical_fluorouracil$most_sensitive  <- ifelse(read_clinical_fluorouracil$PFS < quantile(read_clinical_fluorouracil$PFS, probs = .20), 1, 0)
read_clinical_fluorouracil$least_sensitive <- ifelse(read_clinical_fluorouracil$PFS > quantile(read_clinical_fluorouracil$PFS, probs = .80), 1, 0)

read_gene <- read.csv('Processed_Gene_Expression/read_tcga_rna_seq_processed.csv', row.names = 1)
colnames(read_gene) <- gsub('\\.', '-', colnames(read_gene))
read_matching_idx <- read_clinical_fluorouracil$submitter_id.samples %in% colnames(read_gene)
read_clinical_fluorouracil_short <- read_clinical_fluorouracil[read_matching_idx, ]
read_matching_idx <- colnames(read_gene) %in% read_clinical_fluorouracil_short$submitter_id.samples
read_gene_short <- read_gene[, read_matching_idx]
read_gene_short <- t(read_gene_short)
read_gene_short_scaled <- apply(read_gene_short, 2, scale)

rm(read_tcga_most_min_auc)
rm(read_tcga_most_1se_auc)
rm(read_tcga_least_min_auc)
rm(read_tcga_least_1se_auc)

new_read_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
read_tcga_most_min_auc <- auc(read_clinical_fluorouracil_short$most_sensitive, new_read_tcga_fluorouracil_most_sensitive_min)
read_tcga_most_min_auc <- round(read_tcga_most_min_auc, digits = 2)

new_read_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
read_tcga_most_1se_auc <- auc(read_clinical_fluorouracil_short$most_sensitive, new_read_tcga_fluorouracil_most_sensitive_1se)
read_tcga_most_1se_auc <- round(read_tcga_most_1se_auc, digits = 2)

new_read_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
read_tcga_least_min_auc <- auc(read_clinical_fluorouracil_short$least_sensitive, new_read_tcga_fluorouracil_least_sensitive_min)
read_tcga_least_min_auc <- round(read_tcga_least_min_auc, digits = 2)

new_read_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
read_tcga_least_1se_auc <- auc(read_clinical_fluorouracil_short$least_sensitive, new_read_tcga_fluorouracil_least_sensitive_1se)
read_tcga_least_1se_auc <- round(read_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_read_tcga_fluorouracil_most_sensitive_min, read_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_read_tcga_fluorouracil_most_sensitive_1se, read_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_read_tcga_fluorouracil_least_sensitive_min, read_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_read_tcga_fluorouracil_least_sensitive_1se, read_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/read_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# READ W OXALIPLATIN (18)
read_clinical <- read.csv('Processed_Clinical_Data/read_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(read_clinical$most_sensitive)
read_clinical <- read_clinical[!na_idx, ]
table(read_clinical$drug_name)
read_clinical_oxalaplatin <- read_clinical[which(read_clinical$drug_name == 'Oxaliplatin'), ]

read_clinical_oxalaplatin$most_sensitive  <- ifelse(read_clinical_oxalaplatin$PFS < quantile(read_clinical_oxalaplatin$PFS, probs = .20), 1, 0)
read_clinical_oxalaplatin$least_sensitive <- ifelse(read_clinical_oxalaplatin$PFS > quantile(read_clinical_oxalaplatin$PFS, probs = .80), 1, 0)

read_gene <- read.csv('Processed_Gene_Expression/read_tcga_rna_seq_processed.csv', row.names = 1)
colnames(read_gene) <- gsub('\\.', '-', colnames(read_gene))
read_matching_idx <- read_clinical_oxalaplatin$submitter_id.samples %in% colnames(read_gene)
read_clinical_oxalaplatin_short <- read_clinical_oxalaplatin[read_matching_idx, ]
read_matching_idx <- colnames(read_gene) %in% read_clinical_oxalaplatin_short$submitter_id.samples
read_gene_short <- read_gene[, read_matching_idx]
read_gene_short <- t(read_gene_short)
read_gene_short_scaled <- apply(read_gene_short, 2, scale)

rm(read_tcga_most_min_auc)
rm(read_tcga_most_1se_auc)
rm(read_tcga_least_min_auc)
rm(read_tcga_least_1se_auc)

new_read_tcga_oxalaplatin_most_sensitive_min <- predict(oxalaplatin_ccle_most_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
read_tcga_most_min_auc <- auc(read_clinical_oxalaplatin_short$most_sensitive, new_read_tcga_oxalaplatin_most_sensitive_min)
read_tcga_most_min_auc <- round(read_tcga_most_min_auc, digits = 2)

new_read_tcga_oxalaplatin_most_sensitive_1se <- predict(oxalaplatin_ccle_most_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
read_tcga_most_1se_auc <- auc(read_clinical_oxalaplatin_short$most_sensitive, new_read_tcga_oxalaplatin_most_sensitive_1se)
read_tcga_most_1se_auc <- round(read_tcga_most_1se_auc, digits = 2)

new_read_tcga_oxalaplatin_least_sensitive_min <- predict(oxalaplatin_ccle_least_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
read_tcga_least_min_auc <- auc(read_clinical_oxalaplatin_short$least_sensitive, new_read_tcga_oxalaplatin_least_sensitive_min)
read_tcga_least_min_auc <- round(read_tcga_least_min_auc, digits = 2)

new_read_tcga_oxalaplatin_least_sensitive_1se <- predict(oxalaplatin_ccle_least_fit_elnet, newx = read_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
read_tcga_least_1se_auc <- auc(read_clinical_oxalaplatin_short$least_sensitive, new_read_tcga_oxalaplatin_least_sensitive_1se)
read_tcga_least_1se_auc <- round(read_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_read_tcga_oxalaplatin_most_sensitive_min, read_clinical_oxalaplatin_short$most_sensitive)
pred2 <- prediction(new_read_tcga_oxalaplatin_most_sensitive_1se, read_clinical_oxalaplatin_short$most_sensitive)
pred3 <- prediction(new_read_tcga_oxalaplatin_least_sensitive_min, read_clinical_oxalaplatin_short$least_sensitive)
pred4 <- prediction(new_read_tcga_oxalaplatin_least_sensitive_1se, read_clinical_oxalaplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/read_oxalaplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# SARC W DOCETAXEL (11)
sarc_clinical <- read.csv('Processed_Clinical_Data/sarc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(sarc_clinical$most_sensitive)
sarc_clinical <- sarc_clinical[!na_idx, ]
table(sarc_clinical$drug_name)
sarc_clinical_docetaxel <- sarc_clinical[which(sarc_clinical$drug_name == 'Docetaxel'), ]

sarc_clinical_docetaxel$most_sensitive  <- ifelse(sarc_clinical_docetaxel$PFS < quantile(sarc_clinical_docetaxel$PFS, probs = .20), 1, 0)
sarc_clinical_docetaxel$least_sensitive <- ifelse(sarc_clinical_docetaxel$PFS > quantile(sarc_clinical_docetaxel$PFS, probs = .80), 1, 0)

sarc_gene <- read.csv('Processed_Gene_Expression/sarc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(sarc_gene) <- gsub('\\.', '-', colnames(sarc_gene))
sarc_matching_idx <- sarc_clinical_docetaxel$submitter_id.samples %in% colnames(sarc_gene)
sarc_clinical_docetaxel_short <- sarc_clinical_docetaxel[sarc_matching_idx, ]
sarc_matching_idx <- colnames(sarc_gene) %in% sarc_clinical_docetaxel_short$submitter_id.samples
sarc_gene_short <- sarc_gene[, sarc_matching_idx]
sarc_gene_short <- t(sarc_gene_short)
sarc_gene_short_scaled <- apply(sarc_gene_short, 2, scale)

rm(sarc_tcga_most_min_auc)
rm(sarc_tcga_most_1se_auc)
rm(sarc_tcga_least_min_auc)
rm(sarc_tcga_least_1se_auc)

new_sarc_tcga_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_most_min_auc <- auc(sarc_clinical_docetaxel_short$most_sensitive, new_sarc_tcga_docetaxel_most_sensitive_min)
sarc_tcga_most_min_auc <- round(sarc_tcga_most_min_auc, digits = 2)

new_sarc_tcga_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_most_1se_auc <- auc(sarc_clinical_docetaxel_short$most_sensitive, new_sarc_tcga_docetaxel_most_sensitive_1se)
sarc_tcga_most_1se_auc <- round(sarc_tcga_most_1se_auc, digits = 2)

new_sarc_tcga_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_least_min_auc <- auc(sarc_clinical_docetaxel_short$least_sensitive, new_sarc_tcga_docetaxel_least_sensitive_min)
sarc_tcga_least_min_auc <- round(sarc_tcga_least_min_auc, digits = 2)

new_sarc_tcga_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_least_1se_auc <- auc(sarc_clinical_docetaxel_short$least_sensitive, new_sarc_tcga_docetaxel_least_sensitive_1se)
sarc_tcga_least_1se_auc <- round(sarc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_sarc_tcga_docetaxel_most_sensitive_min, sarc_clinical_docetaxel_short$most_sensitive)
pred2 <- prediction(new_sarc_tcga_docetaxel_most_sensitive_1se, sarc_clinical_docetaxel_short$most_sensitive)
pred3 <- prediction(new_sarc_tcga_docetaxel_least_sensitive_min, sarc_clinical_docetaxel_short$least_sensitive)
pred4 <- prediction(new_sarc_tcga_docetaxel_least_sensitive_1se, sarc_clinical_docetaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/sarc_docetaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# SARC W GEMCITABINE (17)
sarc_clinical <- read.csv('Processed_Clinical_Data/sarc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(sarc_clinical$most_sensitive)
sarc_clinical <- sarc_clinical[!na_idx, ]
table(sarc_clinical$drug_name)
sarc_clinical_gemcitabine <- sarc_clinical[which(sarc_clinical$drug_name == 'gemcitabine' | sarc_clinical$drug_name == 'Gemcitabine' | 
                                                   sarc_clinical$drug_name == 'GEMCITABINE' | sarc_clinical$drug_name == 'Gemzar'), ]

sarc_clinical_gemcitabine$most_sensitive  <- ifelse(sarc_clinical_gemcitabine$PFS < quantile(sarc_clinical_gemcitabine$PFS, probs = .20), 1, 0)
sarc_clinical_gemcitabine$least_sensitive <- ifelse(sarc_clinical_gemcitabine$PFS > quantile(sarc_clinical_gemcitabine$PFS, probs = .80), 1, 0)

sarc_gene <- read.csv('Processed_Gene_Expression/sarc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(sarc_gene) <- gsub('\\.', '-', colnames(sarc_gene))
sarc_matching_idx <- sarc_clinical_gemcitabine$submitter_id.samples %in% colnames(sarc_gene)
sarc_clinical_gemcitabine_short <- sarc_clinical_gemcitabine[sarc_matching_idx, ]
sarc_matching_idx <- colnames(sarc_gene) %in% sarc_clinical_gemcitabine_short$submitter_id.samples
sarc_gene_short <- sarc_gene[, sarc_matching_idx]
sarc_gene_short <- t(sarc_gene_short)
sarc_gene_short_scaled <- apply(sarc_gene_short, 2, scale)

rm(sarc_tcga_most_min_auc)
rm(sarc_tcga_most_1se_auc)
rm(sarc_tcga_least_min_auc)
rm(sarc_tcga_least_1se_auc)

new_sarc_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_most_min_auc <- auc(sarc_clinical_gemcitabine_short$most_sensitive, new_sarc_tcga_gemcitabine_most_sensitive_min)
sarc_tcga_most_min_auc <- round(sarc_tcga_most_min_auc, digits = 2)

new_sarc_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_most_1se_auc <- auc(sarc_clinical_gemcitabine_short$most_sensitive, new_sarc_tcga_gemcitabine_most_sensitive_1se)
sarc_tcga_most_1se_auc <- round(sarc_tcga_most_1se_auc, digits = 2)

new_sarc_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
sarc_tcga_least_min_auc <- auc(sarc_clinical_gemcitabine_short$least_sensitive, new_sarc_tcga_gemcitabine_least_sensitive_min)
sarc_tcga_least_min_auc <- round(sarc_tcga_least_min_auc, digits = 2)

new_sarc_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = sarc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
sarc_tcga_least_1se_auc <- auc(sarc_clinical_gemcitabine_short$least_sensitive, new_sarc_tcga_gemcitabine_least_sensitive_1se)
sarc_tcga_least_1se_auc <- round(sarc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_sarc_tcga_gemcitabine_most_sensitive_min, sarc_clinical_gemcitabine_short$most_sensitive)
pred2 <- prediction(new_sarc_tcga_gemcitabine_most_sensitive_1se, sarc_clinical_gemcitabine_short$most_sensitive)
pred3 <- prediction(new_sarc_tcga_gemcitabine_least_sensitive_min, sarc_clinical_gemcitabine_short$least_sensitive)
pred4 <- prediction(new_sarc_tcga_gemcitabine_least_sensitive_1se, sarc_clinical_gemcitabine_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/sarc_gemcitabine_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[7], lwd = 2, lty = 2, main = 'SARC treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[7], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 1.0)', 'resistant model (AUC = 0.64)'), cex = 0.8, lwd = 2, lty = c(1,2), col = colors_i_need[7], bty = 'n')
dev.off()

# STAD W ETOPOSIDE (7)
stad_clinical <- read.csv('Processed_Clinical_Data/stad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(stad_clinical$most_sensitive)
stad_clinical <- stad_clinical[!na_idx, ]
table(stad_clinical$drug_name)
stad_clinical_etoposide <- stad_clinical[which(stad_clinical$drug_name == 'Etoposide'), ]

stad_clinical_etoposide$most_sensitive  <- ifelse(stad_clinical_etoposide$PFS < quantile(stad_clinical_etoposide$PFS, probs = .20), 1, 0)
stad_clinical_etoposide$least_sensitive <- ifelse(stad_clinical_etoposide$PFS > quantile(stad_clinical_etoposide$PFS, probs = .80), 1, 0)

stad_gene <- read.csv('Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(stad_gene) <- gsub('\\.', '-', colnames(stad_gene))
stad_matching_idx <- stad_clinical_etoposide$submitter_id.samples %in% colnames(stad_gene)
stad_clinical_etoposide_short <- stad_clinical_etoposide[stad_matching_idx, ]
stad_matching_idx <- colnames(stad_gene) %in% stad_clinical_etoposide_short$submitter_id.samples
stad_gene_short <- stad_gene[, stad_matching_idx]
stad_gene_short <- t(stad_gene_short)
stad_gene_short_scaled <- apply(stad_gene_short, 2, scale)

rm(stad_tcga_most_min_auc)
rm(stad_tcga_most_1se_auc)
rm(stad_tcga_least_min_auc)
rm(stad_tcga_least_1se_auc)

new_stad_tcga_etoposide_most_sensitive_min <- predict(etoposide_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_most_min_auc <- auc(stad_clinical_etoposide_short$most_sensitive, new_stad_tcga_etoposide_most_sensitive_min)
stad_tcga_most_min_auc <- round(stad_tcga_most_min_auc, digits = 2)

new_stad_tcga_etoposide_most_sensitive_1se <- predict(etoposide_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_most_1se_auc <- auc(stad_clinical_etoposide_short$most_sensitive, new_stad_tcga_etoposide_most_sensitive_1se)
stad_tcga_most_1se_auc <- round(stad_tcga_most_1se_auc, digits = 2)

new_stad_tcga_etoposide_least_sensitive_min <- predict(etoposide_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_least_min_auc <- auc(stad_clinical_etoposide_short$least_sensitive, new_stad_tcga_etoposide_least_sensitive_min)
stad_tcga_least_min_auc <- round(stad_tcga_least_min_auc, digits = 2)

new_stad_tcga_etoposide_least_sensitive_1se <- predict(etoposide_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_least_1se_auc <- auc(stad_clinical_etoposide_short$least_sensitive, new_stad_tcga_etoposide_least_sensitive_1se)
stad_tcga_least_1se_auc <- round(stad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_stad_tcga_etoposide_most_sensitive_min, stad_clinical_etoposide_short$most_sensitive)
pred2 <- prediction(new_stad_tcga_etoposide_most_sensitive_1se, stad_clinical_etoposide_short$most_sensitive)
pred3 <- prediction(new_stad_tcga_etoposide_least_sensitive_min, stad_clinical_etoposide_short$least_sensitive)
pred4 <- prediction(new_stad_tcga_etoposide_least_sensitive_1se, stad_clinical_etoposide_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/stad_etoposide_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# STAD W FLUOROURACIL (80)
stad_clinical <- read.csv('Processed_Clinical_Data/stad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(stad_clinical$most_sensitive)
stad_clinical <- stad_clinical[!na_idx, ]
table(stad_clinical$drug_name)
stad_clinical_fluorouracil <- stad_clinical[which(stad_clinical$drug_name == '5-Fluorouracil' | stad_clinical$drug_name == '5-fluorouracil' | 
                                                    stad_clinical$drug_name == '5-Flourouracil' | stad_clinical$drug_name == '5-fluorouracilum' | 
                                                    stad_clinical$drug_name == '5-FU' | stad_clinical$drug_name == '5FU' | stad_clinical$drug_name == 'Fluorouracil' | 
                                                    stad_clinical$drug_name == 'fluorouracilum' | stad_clinical$drug_name == 'Fluorouracilum'), ]

stad_clinical_fluorouracil$most_sensitive  <- ifelse(stad_clinical_fluorouracil$PFS < quantile(stad_clinical_fluorouracil$PFS, probs = .20), 1, 0)
stad_clinical_fluorouracil$least_sensitive <- ifelse(stad_clinical_fluorouracil$PFS > quantile(stad_clinical_fluorouracil$PFS, probs = .80), 1, 0)

stad_gene <- read.csv('Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(stad_gene) <- gsub('\\.', '-', colnames(stad_gene))
stad_matching_idx <- stad_clinical_fluorouracil$submitter_id.samples %in% colnames(stad_gene)
stad_clinical_fluorouracil_short <- stad_clinical_fluorouracil[stad_matching_idx, ]
stad_matching_idx <- colnames(stad_gene) %in% stad_clinical_fluorouracil_short$submitter_id.samples
stad_gene_short <- stad_gene[, stad_matching_idx]
stad_gene_short <- t(stad_gene_short)
stad_gene_short_scaled <- apply(stad_gene_short, 2, scale)

rm(stad_tcga_most_min_auc)
rm(stad_tcga_most_1se_auc)
rm(stad_tcga_least_min_auc)
rm(stad_tcga_least_1se_auc)

new_stad_tcga_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_most_min_auc <- auc(stad_clinical_fluorouracil_short$most_sensitive, new_stad_tcga_fluorouracil_most_sensitive_min)
stad_tcga_most_min_auc <- round(stad_tcga_most_min_auc, digits = 2)

new_stad_tcga_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_most_1se_auc <- auc(stad_clinical_fluorouracil_short$most_sensitive, new_stad_tcga_fluorouracil_most_sensitive_1se)
stad_tcga_most_1se_auc <- round(stad_tcga_most_1se_auc, digits = 2)

new_stad_tcga_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_least_min_auc <- auc(stad_clinical_fluorouracil_short$least_sensitive, new_stad_tcga_fluorouracil_least_sensitive_min)
stad_tcga_least_min_auc <- round(stad_tcga_least_min_auc, digits = 2)

new_stad_tcga_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_least_1se_auc <- auc(stad_clinical_fluorouracil_short$least_sensitive, new_stad_tcga_fluorouracil_least_sensitive_1se)
stad_tcga_least_1se_auc <- round(stad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_stad_tcga_fluorouracil_most_sensitive_min, stad_clinical_fluorouracil_short$most_sensitive)
pred2 <- prediction(new_stad_tcga_fluorouracil_most_sensitive_1se, stad_clinical_fluorouracil_short$most_sensitive)
pred3 <- prediction(new_stad_tcga_fluorouracil_least_sensitive_min, stad_clinical_fluorouracil_short$least_sensitive)
pred4 <- prediction(new_stad_tcga_fluorouracil_least_sensitive_1se, stad_clinical_fluorouracil_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/stad_fluorouracil_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# STAD W OXALIPLATIN (7)
stad_clinical <- read.csv('Processed_Clinical_Data/stad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(stad_clinical$most_sensitive)
stad_clinical <- stad_clinical[!na_idx, ]
table(stad_clinical$drug_name)
stad_clinical_oxalaplatin <- stad_clinical[which(stad_clinical$drug_name == 'Oxaliplatin' | stad_clinical$drug_name == 'Oxaloplatin'), ]

stad_clinical_oxalaplatin$most_sensitive  <- ifelse(stad_clinical_oxalaplatin$PFS < quantile(stad_clinical_oxalaplatin$PFS, probs = .20), 1, 0)
stad_clinical_oxalaplatin$least_sensitive <- ifelse(stad_clinical_oxalaplatin$PFS > quantile(stad_clinical_oxalaplatin$PFS, probs = .80), 1, 0)

stad_gene <- read.csv('Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(stad_gene) <- gsub('\\.', '-', colnames(stad_gene))
stad_matching_idx <- stad_clinical_oxalaplatin$submitter_id.samples %in% colnames(stad_gene)
stad_clinical_oxalaplatin_short <- stad_clinical_oxalaplatin[stad_matching_idx, ]
stad_matching_idx <- colnames(stad_gene) %in% stad_clinical_oxalaplatin_short$submitter_id.samples
stad_gene_short <- stad_gene[, stad_matching_idx]
stad_gene_short <- t(stad_gene_short)
stad_gene_short_scaled <- apply(stad_gene_short, 2, scale)

rm(stad_tcga_most_min_auc)
rm(stad_tcga_most_1se_auc)
rm(stad_tcga_least_min_auc)
rm(stad_tcga_least_1se_auc)

new_stad_tcga_oxalaplatin_most_sensitive_min <- predict(oxalaplatin_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_most_min_auc <- auc(stad_clinical_oxalaplatin_short$most_sensitive, new_stad_tcga_oxalaplatin_most_sensitive_min)
stad_tcga_most_min_auc <- round(stad_tcga_most_min_auc, digits = 2)

new_stad_tcga_oxalaplatin_most_sensitive_1se <- predict(oxalaplatin_ccle_most_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_most_1se_auc <- auc(stad_clinical_oxalaplatin_short$most_sensitive, new_stad_tcga_oxalaplatin_most_sensitive_1se)
stad_tcga_most_1se_auc <- round(stad_tcga_most_1se_auc, digits = 2)

new_stad_tcga_oxalaplatin_least_sensitive_min <- predict(oxalaplatin_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
stad_tcga_least_min_auc <- auc(stad_clinical_oxalaplatin_short$least_sensitive, new_stad_tcga_oxalaplatin_least_sensitive_min)
stad_tcga_least_min_auc <- round(stad_tcga_least_min_auc, digits = 2)

new_stad_tcga_oxalaplatin_least_sensitive_1se <- predict(oxalaplatin_ccle_least_fit_elnet, newx = stad_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
stad_tcga_least_1se_auc <- auc(stad_clinical_oxalaplatin_short$least_sensitive, new_stad_tcga_oxalaplatin_least_sensitive_1se)
stad_tcga_least_1se_auc <- round(stad_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_stad_tcga_oxalaplatin_most_sensitive_min, stad_clinical_oxalaplatin_short$most_sensitive)
pred2 <- prediction(new_stad_tcga_oxalaplatin_most_sensitive_1se, stad_clinical_oxalaplatin_short$most_sensitive)
pred3 <- prediction(new_stad_tcga_oxalaplatin_least_sensitive_min, stad_clinical_oxalaplatin_short$least_sensitive)
pred4 <- prediction(new_stad_tcga_oxalaplatin_least_sensitive_1se, stad_clinical_oxalaplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/stad_oxalaplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCEC W CARBOPLATIN (52)
ucec_clinical <- read.csv('Processed_Clinical_Data/ucec_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucec_clinical$most_sensitive)
ucec_clinical <- ucec_clinical[!na_idx, ]
table(ucec_clinical$drug_name)
ucec_clinical_carboplatin <- ucec_clinical[which(ucec_clinical$drug_name == 'Carboplatin' | ucec_clinical$drug_name == 'carboplatin' | 
                                                   ucec_clinical$drug_name == 'CARBOPLATIN'), ]

ucec_clinical_carboplatin$most_sensitive  <- ifelse(ucec_clinical_carboplatin$PFS < quantile(ucec_clinical_carboplatin$PFS, probs = .20), 1, 0)
ucec_clinical_carboplatin$least_sensitive <- ifelse(ucec_clinical_carboplatin$PFS > quantile(ucec_clinical_carboplatin$PFS, probs = .80), 1, 0)

ucec_gene <- read.csv('Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucec_gene) <- gsub('\\.', '-', colnames(ucec_gene))
ucec_matching_idx <- ucec_clinical_carboplatin$submitter_id.samples %in% colnames(ucec_gene)
ucec_clinical_carboplatin_short <- ucec_clinical_carboplatin[ucec_matching_idx, ]
ucec_matching_idx <- colnames(ucec_gene) %in% ucec_clinical_carboplatin_short$submitter_id.samples
ucec_gene_short <- ucec_gene[, ucec_matching_idx]
ucec_gene_short <- t(ucec_gene_short)
ucec_gene_short_scaled <- apply(ucec_gene_short, 2, scale)

rm(ucec_tcga_most_min_auc)
rm(ucec_tcga_most_1se_auc)
rm(ucec_tcga_least_min_auc)
rm(ucec_tcga_least_1se_auc)

new_ucec_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_most_min_auc <- auc(ucec_clinical_carboplatin_short$most_sensitive, new_ucec_tcga_carboplatin_most_sensitive_min)
ucec_tcga_most_min_auc <- round(ucec_tcga_most_min_auc, digits = 2)

new_ucec_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_most_1se_auc <- auc(ucec_clinical_carboplatin_short$most_sensitive, new_ucec_tcga_carboplatin_most_sensitive_1se)
ucec_tcga_most_1se_auc <- round(ucec_tcga_most_1se_auc, digits = 2)

new_ucec_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_least_min_auc <- auc(ucec_clinical_carboplatin_short$least_sensitive, new_ucec_tcga_carboplatin_least_sensitive_min)
ucec_tcga_least_min_auc <- round(ucec_tcga_least_min_auc, digits = 2)

new_ucec_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_least_1se_auc <- auc(ucec_clinical_carboplatin_short$least_sensitive, new_ucec_tcga_carboplatin_least_sensitive_1se)
ucec_tcga_least_1se_auc <- round(ucec_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucec_tcga_carboplatin_most_sensitive_min, ucec_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_ucec_tcga_carboplatin_most_sensitive_1se, ucec_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_ucec_tcga_carboplatin_least_sensitive_min, ucec_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_ucec_tcga_carboplatin_least_sensitive_1se, ucec_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucec_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCEC W DOXORUBICIN (7)
ucec_clinical <- read.csv('Processed_Clinical_Data/ucec_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucec_clinical$most_sensitive)
ucec_clinical <- ucec_clinical[!na_idx, ]
table(ucec_clinical$drug_name)
ucec_clinical_doxorubicin <- ucec_clinical[which(ucec_clinical$drug_name == 'doxorubicin' | ucec_clinical$drug_name == 'Doxorubicin'), ]

ucec_clinical_doxorubicin$most_sensitive  <- ifelse(ucec_clinical_doxorubicin$PFS < quantile(ucec_clinical_doxorubicin$PFS, probs = .20), 1, 0)
ucec_clinical_doxorubicin$least_sensitive <- ifelse(ucec_clinical_doxorubicin$PFS > quantile(ucec_clinical_doxorubicin$PFS, probs = .80), 1, 0)

ucec_gene <- read.csv('Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucec_gene) <- gsub('\\.', '-', colnames(ucec_gene))
ucec_matching_idx <- ucec_clinical_doxorubicin$submitter_id.samples %in% colnames(ucec_gene)
ucec_clinical_doxorubicin_short <- ucec_clinical_doxorubicin[ucec_matching_idx, ]
ucec_matching_idx <- colnames(ucec_gene) %in% ucec_clinical_doxorubicin_short$submitter_id.samples
ucec_gene_short <- ucec_gene[, ucec_matching_idx]
ucec_gene_short <- t(ucec_gene_short)
ucec_gene_short_scaled <- apply(ucec_gene_short, 2, scale)

rm(ucec_tcga_most_min_auc)
rm(ucec_tcga_most_1se_auc)
rm(ucec_tcga_least_min_auc)
rm(ucec_tcga_least_1se_auc)

new_ucec_tcga_doxorubicin_most_sensitive_min <- predict(doxorubicin_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_most_min_auc <- auc(ucec_clinical_doxorubicin_short$most_sensitive, new_ucec_tcga_doxorubicin_most_sensitive_min)
ucec_tcga_most_min_auc <- round(ucec_tcga_most_min_auc, digits = 2)

new_ucec_tcga_doxorubicin_most_sensitive_1se <- predict(doxorubicin_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_most_1se_auc <- auc(ucec_clinical_doxorubicin_short$most_sensitive, new_ucec_tcga_doxorubicin_most_sensitive_1se)
ucec_tcga_most_1se_auc <- round(ucec_tcga_most_1se_auc, digits = 2)

new_ucec_tcga_doxorubicin_least_sensitive_min <- predict(doxorubicin_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_least_min_auc <- auc(ucec_clinical_doxorubicin_short$least_sensitive, new_ucec_tcga_doxorubicin_least_sensitive_min)
ucec_tcga_least_min_auc <- round(ucec_tcga_least_min_auc, digits = 2)

new_ucec_tcga_doxorubicin_least_sensitive_1se <- predict(doxorubicin_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_least_1se_auc <- auc(ucec_clinical_doxorubicin_short$least_sensitive, new_ucec_tcga_doxorubicin_least_sensitive_1se)
ucec_tcga_least_1se_auc <- round(ucec_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucec_tcga_doxorubicin_most_sensitive_min, ucec_clinical_doxorubicin_short$most_sensitive)
pred2 <- prediction(new_ucec_tcga_doxorubicin_most_sensitive_1se, ucec_clinical_doxorubicin_short$most_sensitive)
pred3 <- prediction(new_ucec_tcga_doxorubicin_least_sensitive_min, ucec_clinical_doxorubicin_short$least_sensitive)
pred4 <- prediction(new_ucec_tcga_doxorubicin_least_sensitive_1se, ucec_clinical_doxorubicin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucec_doxorubicin_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCEC W PACLITAXEL (49)
ucec_clinical <- read.csv('Processed_Clinical_Data/ucec_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucec_clinical$most_sensitive)
ucec_clinical <- ucec_clinical[!na_idx, ]
table(ucec_clinical$drug_name)
ucec_clinical_paclitaxel <- ucec_clinical[which(ucec_clinical$drug_name == 'paclitaxel' | ucec_clinical$drug_name == 'Paclitaxel' | 
                                                  ucec_clinical$drug_name == 'PACLITAXEL' | ucec_clinical$drug_name == 'Paclitaxol'), ]

ucec_clinical_paclitaxel$most_sensitive  <- ifelse(ucec_clinical_paclitaxel$PFS < quantile(ucec_clinical_paclitaxel$PFS, probs = .20), 1, 0)
ucec_clinical_paclitaxel$least_sensitive <- ifelse(ucec_clinical_paclitaxel$PFS > quantile(ucec_clinical_paclitaxel$PFS, probs = .80), 1, 0)

ucec_gene <- read.csv('Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucec_gene) <- gsub('\\.', '-', colnames(ucec_gene))
ucec_matching_idx <- ucec_clinical_paclitaxel$submitter_id.samples %in% colnames(ucec_gene)
ucec_clinical_paclitaxel_short <- ucec_clinical_paclitaxel[ucec_matching_idx, ]
ucec_matching_idx <- colnames(ucec_gene) %in% ucec_clinical_paclitaxel_short$submitter_id.samples
ucec_gene_short <- ucec_gene[, ucec_matching_idx]
ucec_gene_short <- t(ucec_gene_short)
ucec_gene_short_scaled <- apply(ucec_gene_short, 2, scale)

rm(ucec_tcga_most_min_auc)
rm(ucec_tcga_most_1se_auc)
rm(ucec_tcga_least_min_auc)
rm(ucec_tcga_least_1se_auc)

new_ucec_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_most_min_auc <- auc(ucec_clinical_paclitaxel_short$most_sensitive, new_ucec_tcga_paclitaxel_most_sensitive_min)
ucec_tcga_most_min_auc <- round(ucec_tcga_most_min_auc, digits = 2)

new_ucec_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_most_1se_auc <- auc(ucec_clinical_paclitaxel_short$most_sensitive, new_ucec_tcga_paclitaxel_most_sensitive_1se)
ucec_tcga_most_1se_auc <- round(ucec_tcga_most_1se_auc, digits = 2)

new_ucec_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucec_tcga_least_min_auc <- auc(ucec_clinical_paclitaxel_short$least_sensitive, new_ucec_tcga_paclitaxel_least_sensitive_min)
ucec_tcga_least_min_auc <- round(ucec_tcga_least_min_auc, digits = 2)

new_ucec_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = ucec_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucec_tcga_least_1se_auc <- auc(ucec_clinical_paclitaxel_short$least_sensitive, new_ucec_tcga_paclitaxel_least_sensitive_1se)
ucec_tcga_least_1se_auc <- round(ucec_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucec_tcga_paclitaxel_most_sensitive_min, ucec_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_ucec_tcga_paclitaxel_most_sensitive_1se, ucec_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_ucec_tcga_paclitaxel_least_sensitive_min, ucec_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_ucec_tcga_paclitaxel_least_sensitive_1se, ucec_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucec_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()

# UCS W CARBOPLATIN (16)
ucs_clinical <- read.csv('Processed_Clinical_Data/ucs_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucs_clinical$most_sensitive)
ucs_clinical <- ucs_clinical[!na_idx, ]
table(ucs_clinical$drug_name)
ucs_clinical_carboplatin <- ucs_clinical[which(ucs_clinical$drug_name == 'Carboplatin'), ]

ucs_clinical_carboplatin$most_sensitive  <- ifelse(ucs_clinical_carboplatin$PFS < quantile(ucs_clinical_carboplatin$PFS, probs = .20), 1, 0)
ucs_clinical_carboplatin$least_sensitive <- ifelse(ucs_clinical_carboplatin$PFS > quantile(ucs_clinical_carboplatin$PFS, probs = .80), 1, 0)

ucs_gene <- read.csv('Processed_Gene_Expression/ucs_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucs_gene) <- gsub('\\.', '-', colnames(ucs_gene))
ucs_matching_idx <- ucs_clinical_carboplatin$submitter_id.samples %in% colnames(ucs_gene)
ucs_clinical_carboplatin_short <- ucs_clinical_carboplatin[ucs_matching_idx, ]
ucs_matching_idx <- colnames(ucs_gene) %in% ucs_clinical_carboplatin_short$submitter_id.samples
ucs_gene_short <- ucs_gene[, ucs_matching_idx]
ucs_gene_short <- t(ucs_gene_short)
ucs_gene_short_scaled <- apply(ucs_gene_short, 2, scale)

rm(ucs_tcga_most_min_auc)
rm(ucs_tcga_most_1se_auc)
rm(ucs_tcga_least_min_auc)
rm(ucs_tcga_least_1se_auc)

new_ucs_tcga_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_most_min_auc <- auc(ucs_clinical_carboplatin_short$most_sensitive, new_ucs_tcga_carboplatin_most_sensitive_min)
ucs_tcga_most_min_auc <- round(ucs_tcga_most_min_auc, digits = 2)

new_ucs_tcga_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_most_1se_auc <- auc(ucs_clinical_carboplatin_short$most_sensitive, new_ucs_tcga_carboplatin_most_sensitive_1se)
ucs_tcga_most_1se_auc <- round(ucs_tcga_most_1se_auc, digits = 2)

new_ucs_tcga_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_least_min_auc <- auc(ucs_clinical_carboplatin_short$least_sensitive, new_ucs_tcga_carboplatin_least_sensitive_min)
ucs_tcga_least_min_auc <- round(ucs_tcga_least_min_auc, digits = 2)

new_ucs_tcga_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_least_1se_auc <- auc(ucs_clinical_carboplatin_short$least_sensitive, new_ucs_tcga_carboplatin_least_sensitive_1se)
ucs_tcga_least_1se_auc <- round(ucs_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucs_tcga_carboplatin_most_sensitive_min, ucs_clinical_carboplatin_short$most_sensitive)
pred2 <- prediction(new_ucs_tcga_carboplatin_most_sensitive_1se, ucs_clinical_carboplatin_short$most_sensitive)
pred3 <- prediction(new_ucs_tcga_carboplatin_least_sensitive_min, ucs_clinical_carboplatin_short$least_sensitive)
pred4 <- prediction(new_ucs_tcga_carboplatin_least_sensitive_1se, ucs_clinical_carboplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucs_carboplatin_tcga_ccle_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = colors_i_need[1], lwd = 2, lty = 2, main = 'UCS treated with carboplatin')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = colors_i_need[1], lwd = 2, lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.6, y = 0.4, legend = c('sensitive model (AUC = 0.59)', 'resistant model (AUC = 0.80)'), cex = 0.8, lwd = 2, lty = c(1,2), col = colors_i_need[1], bty = 'n')
dev.off()

# UCS W PACLITAXEL (10)
ucs_clinical <- read.csv('Processed_Clinical_Data/ucs_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(ucs_clinical$most_sensitive)
ucs_clinical <- ucs_clinical[!na_idx, ]
table(ucs_clinical$drug_name)
ucs_clinical_paclitaxel <- ucs_clinical[which(ucs_clinical$drug_name == 'Paclitaxel'), ]

ucs_clinical_paclitaxel$most_sensitive  <- ifelse(ucs_clinical_paclitaxel$PFS < quantile(ucs_clinical_paclitaxel$PFS, probs = .20), 1, 0)
ucs_clinical_paclitaxel$least_sensitive <- ifelse(ucs_clinical_paclitaxel$PFS > quantile(ucs_clinical_paclitaxel$PFS, probs = .80), 1, 0)

ucs_gene <- read.csv('Processed_Gene_Expression/ucs_tcga_rna_seq_processed.csv', row.names = 1)
colnames(ucs_gene) <- gsub('\\.', '-', colnames(ucs_gene))
ucs_matching_idx <- ucs_clinical_paclitaxel$submitter_id.samples %in% colnames(ucs_gene)
ucs_clinical_paclitaxel_short <- ucs_clinical_paclitaxel[ucs_matching_idx, ]
ucs_matching_idx <- colnames(ucs_gene) %in% ucs_clinical_paclitaxel_short$submitter_id.samples
ucs_gene_short <- ucs_gene[, ucs_matching_idx]
ucs_gene_short <- t(ucs_gene_short)
ucs_gene_short_scaled <- apply(ucs_gene_short, 2, scale)

rm(ucs_tcga_most_min_auc)
rm(ucs_tcga_most_1se_auc)
rm(ucs_tcga_least_min_auc)
rm(ucs_tcga_least_1se_auc)

new_ucs_tcga_paclitaxel_most_sensitive_min <- predict(paclitaxel_ccle_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_most_min_auc <- auc(ucs_clinical_paclitaxel_short$most_sensitive, new_ucs_tcga_paclitaxel_most_sensitive_min)
ucs_tcga_most_min_auc <- round(ucs_tcga_most_min_auc, digits = 2)

new_ucs_tcga_paclitaxel_most_sensitive_1se <- predict(paclitaxel_ccle_most_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_most_1se_auc <- auc(ucs_clinical_paclitaxel_short$most_sensitive, new_ucs_tcga_paclitaxel_most_sensitive_1se)
ucs_tcga_most_1se_auc <- round(ucs_tcga_most_1se_auc, digits = 2)

new_ucs_tcga_paclitaxel_least_sensitive_min <- predict(paclitaxel_ccle_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
ucs_tcga_least_min_auc <- auc(ucs_clinical_paclitaxel_short$least_sensitive, new_ucs_tcga_paclitaxel_least_sensitive_min)
ucs_tcga_least_min_auc <- round(ucs_tcga_least_min_auc, digits = 2)

new_ucs_tcga_paclitaxel_least_sensitive_1se <- predict(paclitaxel_ccle_least_fit_elnet, newx = ucs_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
ucs_tcga_least_1se_auc <- auc(ucs_clinical_paclitaxel_short$least_sensitive, new_ucs_tcga_paclitaxel_least_sensitive_1se)
ucs_tcga_least_1se_auc <- round(ucs_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_ucs_tcga_paclitaxel_most_sensitive_min, ucs_clinical_paclitaxel_short$most_sensitive)
pred2 <- prediction(new_ucs_tcga_paclitaxel_most_sensitive_1se, ucs_clinical_paclitaxel_short$most_sensitive)
pred3 <- prediction(new_ucs_tcga_paclitaxel_least_sensitive_min, ucs_clinical_paclitaxel_short$least_sensitive)
pred4 <- prediction(new_ucs_tcga_paclitaxel_least_sensitive_1se, ucs_clinical_paclitaxel_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/ucs_paclitaxel_tcga_ccle_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
dev.off()