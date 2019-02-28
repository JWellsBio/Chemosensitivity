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

### GDSC ------
## load clinical data ----
cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)
etoposide     <- read.csv('Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv', row.names = 1)
gemcitabine   <- read.csv('Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv', row.names = 1)
methotrexate  <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)

# GDSC cell line info
gdsc_cell_line_info <- read.csv('Clinical_Files/Cell_listTue Dec 18 19_10_36 2018 (version 1).csv', header = TRUE, stringsAsFactors = FALSE)
gdsc_cell_line_info <- gdsc_cell_line_info[, c(1,2,5,6)]
colnames(gdsc_cell_line_info)[c(3,4)] <- c('TCGA_class', 'Tissue_subtype')

## load gene expression data ----
gdsc_rna_seq <- read.csv('Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE)
gdsc_rna_seq <- gdsc_rna_seq[, -1]
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))

### set up data for model building ----------
# get names of GDSC cell lines treated with each drug
cisplatin_lines           <- cisplatin$COSMIC_ID #680
etoposide_lines           <- etoposide$COSMIC_ID #718
gemcitabine_lines         <- gemcitabine$COSMIC_ID #707
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
gemcitabine_rna_seq_train   <- gdsc_train[intersect(gemcitabine_lines, rownames(gdsc_train)), ]
# 347 x 14209
gemcitabine_rna_seq_test    <- gdsc_test[intersect(gemcitabine_lines, rownames(gdsc_test)), ]
# 360 x 14209
methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# 337 x 14209
methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# 342 x 14209

# split clinical data
cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]

gemcitabine_train      <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test       <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_test)), ]

methotrexate_train     <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test      <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]

# scale data
cisplatin_rna_seq_train_scaled          <- apply(cisplatin_rna_seq_train, 2, scale)
cisplatin_rna_seq_test_scaled           <- as.data.frame(apply(cisplatin_rna_seq_test, 2, scale))

etoposide_rna_seq_train_scaled          <- apply(etoposide_rna_seq_train, 2, scale)
etoposide_rna_seq_test_scaled           <- as.data.frame(apply(etoposide_rna_seq_test, 2, scale))

gemcitabine_rna_seq_train_scaled        <- apply(gemcitabine_rna_seq_train, 2, scale)
gemcitabine_rna_seq_test_scaled         <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))

methotrexate_rna_seq_train_scaled       <- apply(methotrexate_rna_seq_train, 2, scale)
methotrexate_rna_seq_test_scaled        <- as.data.frame(apply(methotrexate_rna_seq_test, 2, scale))

### load models ----
cisplatin_most_fit_elnet      <- readRDS('GLM_Models/cisplatin_most_model.rds')
cisplatin_least_fit_elnet     <- readRDS('GLM_Models/cisplatin_least_model.rds')

etoposide_most_fit_elnet      <- readRDS('GLM_Models/etoposide_most_model.rds')
etoposide_least_fit_elnet     <- readRDS('GLM_Models/etoposide_least_model.rds')

gemcitabine_most_fit_elnet    <- readRDS('GLM_Models/gemcitabine_most_model.rds')
gemcitabine_least_fit_elnet   <- readRDS('GLM_Models/gemcitabine_least_model.rds')

methotrexate_most_fit_elnet   <- readRDS('GLM_Models/methotrexate_most_model.rds')
methotrexate_least_fit_elnet  <- readRDS('GLM_Models/methotrexate_least_model.rds')

### get accuracy on testing data --------

## CISPLATIN
#new_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#cisplatin_most_test_gdsc_auc_min <- auc(cisplatin_test$most_sensitive, new_cisplatin_most_sensitive_min)
#cisplatin_most_test_gdsc_auc_min <- round(cisplatin_most_test_gdsc_auc_min, digits = 2)

new_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_most_test_gdsc_auc_1se <- auc(cisplatin_test$most_sensitive, new_cisplatin_most_sensitive_1se)
cisplatin_most_test_gdsc_auc_1se <- round(cisplatin_most_test_gdsc_auc_1se, digits = 2)


new_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_least_test_gdsc_auc_min <- auc(cisplatin_test$least_sensitive, new_cisplatin_least_sensitive_min)
cisplatin_least_test_gdsc_auc_min <- round(cisplatin_least_test_gdsc_auc_min, digits = 2)

#new_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

#cisplatin_least_test_gdsc_auc_1se <- auc(cisplatin_test$least_sensitive, new_cisplatin_least_sensitive_1se)
#cisplatin_least_test_gdsc_auc_1se <- round(cisplatin_least_test_gdsc_auc_1se, digits = 2)


## ETOPOSIDE
new_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_most_test_gdsc_auc_min <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_min)
etoposide_most_test_gdsc_auc_min <- round(etoposide_most_test_gdsc_auc_min, digits = 2)

#new_etoposide_most_sensitive_1se <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

#etoposide_most_test_gdsc_auc_1se <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_1se)
#etoposide_most_test_gdsc_auc_1se <- round(etoposide_most_test_gdsc_auc_1se, digits = 2)


new_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_least_test_gdsc_auc_min <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_min)
etoposide_least_test_gdsc_auc_min <- round(etoposide_least_test_gdsc_auc_min, digits = 2)

#new_etoposide_least_sensitive_1se <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

#etoposide_least_test_gdsc_auc_1se <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_1se)
#etoposide_least_test_gdsc_auc_1se <- round(etoposide_least_test_gdsc_auc_1se, digits = 2)


## GEMCITABINE
#new_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#gemcitabine_most_test_gdsc_auc_min <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_min)
#gemcitabine_most_test_gdsc_auc_min <- round(gemcitabine_most_test_gdsc_auc_min, digits = 2)

new_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_most_test_gdsc_auc_1se <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_1se)
gemcitabine_most_test_gdsc_auc_1se <- round(gemcitabine_most_test_gdsc_auc_1se, digits = 2)


#new_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#gemcitabine_least_test_gdsc_auc_min <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_min)
#gemcitabine_least_test_gdsc_auc_min <- round(gemcitabine_least_test_gdsc_auc_min, digits = 2)

new_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_least_test_gdsc_auc_1se <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_1se)
gemcitabine_least_test_gdsc_auc_1se <- round(gemcitabine_least_test_gdsc_auc_1se, digits = 2)


## METHOTREXATE
new_methotrexate_most_sensitive_min <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_most_test_gdsc_auc_min <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_min)
methotrexate_most_test_gdsc_auc_min <- round(methotrexate_most_test_gdsc_auc_min, digits = 2)

#new_methotrexate_most_sensitive_1se <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

#methotrexate_most_test_gdsc_auc_1se <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_1se)
#methotrexate_most_test_gdsc_auc_1se <- round(methotrexate_most_test_gdsc_auc_1se, digits = 2)


#new_methotrexate_least_sensitive_min <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#methotrexate_least_test_gdsc_auc_min <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_min)
#methotrexate_least_test_gdsc_auc_min <- round(methotrexate_least_test_gdsc_auc_min, digits = 2)

new_methotrexate_least_sensitive_1se <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_least_test_gdsc_auc_1se <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_1se)
methotrexate_least_test_gdsc_auc_1se <- round(methotrexate_least_test_gdsc_auc_1se, digits = 2)

overall_auc <- c(cisplatin_most_test_gdsc_auc_1se, cisplatin_least_test_gdsc_auc_min, 
                 etoposide_most_test_gdsc_auc_min, etoposide_least_test_gdsc_auc_min, 
                 gemcitabine_most_test_gdsc_auc_1se, gemcitabine_least_test_gdsc_auc_1se, 
                 methotrexate_most_test_gdsc_auc_min, methotrexate_least_test_gdsc_auc_1se)

# subsetting lines by cancer type
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

cisplatin_test_scaled_aero_dig_tract <- cisplatin_rna_seq_test_scaled[cisplatin_aero_dig_tract_lines, ]
cisplatin_test_aero_dig_tract <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_aero_dig_tract), s = 'lambda.1se', interval = 'conf')
cisplatin_aero_dig_tract_most_auc <- auc(cisplatin_test_aero_dig_tract$most_sensitive, new_ic50)

cisplatin_test_scaled_bone <- cisplatin_rna_seq_test_scaled[cisplatin_bone_lines, ]
cisplatin_test_bone <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
cisplatin_bone_most_auc <- auc(cisplatin_test_bone$most_sensitive, new_ic50)

cisplatin_test_scaled_breast <- cisplatin_rna_seq_test_scaled[cisplatin_breast_lines, ]
cisplatin_test_breast <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
cisplatin_breast_most_auc <- auc(cisplatin_test_breast$most_sensitive, new_ic50)

cisplatin_test_scaled_digestive_system <- cisplatin_rna_seq_test_scaled[cisplatin_digestive_system_lines, ]
cisplatin_test_digestive_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_digestive_system), s = 'lambda.1se', interval = 'conf')
cisplatin_digestive_system_most_auc <- auc(cisplatin_test_digestive_system$most_sensitive, new_ic50)

cisplatin_test_scaled_kidney <- cisplatin_rna_seq_test_scaled[cisplatin_kidney_lines, ]
cisplatin_test_kidney <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
cisplatin_kidney_most_auc <- auc(cisplatin_test_kidney$most_sensitive, new_ic50)

cisplatin_test_scaled_large_intestine <- cisplatin_rna_seq_test_scaled[cisplatin_large_intestine_lines, ]
cisplatin_test_large_intestine <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
cisplatin_large_intestine_most_auc <- auc(cisplatin_test_large_intestine$most_sensitive, new_ic50)

cisplatin_test_scaled_lung <- cisplatin_rna_seq_test_scaled[cisplatin_lung_lines, ]
cisplatin_test_lung <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
cisplatin_lung_most_auc <- auc(cisplatin_test_lung$most_sensitive, new_ic50)

cisplatin_test_scaled_lung_NSCLC <- cisplatin_rna_seq_test_scaled[cisplatin_lung_NSCLC_lines, ]
cisplatin_test_lung_NSCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung_NSCLC), s = 'lambda.1se', interval = 'conf')
cisplatin_lung_NSCLC_most_auc <- auc(cisplatin_test_lung_NSCLC$most_sensitive, new_ic50)

cisplatin_test_scaled_lung_SCLC <- cisplatin_rna_seq_test_scaled[cisplatin_lung_SCLC_lines, ]
cisplatin_test_lung_SCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung_SCLC), s = 'lambda.1se', interval = 'conf')
cisplatin_lung_SCLC_most_auc <- auc(cisplatin_test_lung_SCLC$most_sensitive, new_ic50)

cisplatin_test_scaled_nervous_system <- cisplatin_rna_seq_test_scaled[cisplatin_nervous_system_lines, ]
cisplatin_test_nervous_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_nervous_system), s = 'lambda.1se', interval = 'conf')
cisplatin_nervous_system_most_auc <- auc(cisplatin_test_nervous_system$most_sensitive, new_ic50)

cisplatin_test_scaled_neuroblastoma <- cisplatin_rna_seq_test_scaled[cisplatin_neuroblastoma_lines, ]
cisplatin_test_neuroblastoma <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_neuroblastoma), s = 'lambda.1se', interval = 'conf')
cisplatin_neuroblastoma_most_auc <- auc(cisplatin_test_neuroblastoma$most_sensitive, new_ic50)

cisplatin_test_scaled_pancreas <- cisplatin_rna_seq_test_scaled[cisplatin_pancreas_lines, ]
cisplatin_test_pancreas <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
cisplatin_pancreas_most_auc <- auc(cisplatin_test_pancreas$most_sensitive, new_ic50)

cisplatin_test_scaled_skin <- cisplatin_rna_seq_test_scaled[cisplatin_skin_lines, ]
cisplatin_test_skin <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
cisplatin_skin_most_auc <- auc(cisplatin_test_skin$most_sensitive, new_ic50)

cisplatin_test_scaled_soft_tissue <- cisplatin_rna_seq_test_scaled[cisplatin_soft_tissue_lines, ]
cisplatin_test_soft_tissue <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
cisplatin_soft_tissue_most_auc <- auc(cisplatin_test_soft_tissue$most_sensitive, new_ic50)

cisplatin_test_scaled_thyroid <- cisplatin_rna_seq_test_scaled[cisplatin_thyroid_lines, ]
cisplatin_test_thyroid <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
cisplatin_thyroid_most_auc <- auc(cisplatin_test_thyroid$most_sensitive, new_ic50)

cisplatin_test_scaled_urogenital_system <- cisplatin_rna_seq_test_scaled[cisplatin_urogenital_system_lines, ]
cisplatin_test_urogenital_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_test_scaled_urogenital_system), s = 'lambda.1se', interval = 'conf')
cisplatin_urogenital_system_most_auc <- auc(cisplatin_test_urogenital_system$most_sensitive, new_ic50)

cisplatin_most_auc <- c(cisplatin_aero_dig_tract_most_auc, cisplatin_bone_most_auc, 
                        cisplatin_breast_most_auc, cisplatin_digestive_system_most_auc, 
                        cisplatin_kidney_most_auc, cisplatin_large_intestine_most_auc, 
                        cisplatin_lung_most_auc, cisplatin_lung_NSCLC_most_auc, cisplatin_lung_SCLC_most_auc, 
                        cisplatin_nervous_system_most_auc, cisplatin_neuroblastoma_most_auc, 
                        cisplatin_pancreas_most_auc, cisplatin_skin_most_auc, 
                        cisplatin_soft_tissue_most_auc, cisplatin_thyroid_most_auc, 
                        cisplatin_urogenital_system_most_auc)

cisplatin_test_scaled_aero_dig_tract <- cisplatin_rna_seq_test_scaled[cisplatin_aero_dig_tract_lines, ]
cisplatin_test_aero_dig_tract <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_aero_dig_tract), s = 'lambda.min', interval = 'conf')
cisplatin_aero_dig_tract_least_auc <- auc(cisplatin_test_aero_dig_tract$least_sensitive, new_ic50)

cisplatin_test_scaled_bone <- cisplatin_rna_seq_test_scaled[cisplatin_bone_lines, ]
cisplatin_test_bone <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_bone), s = 'lambda.min', interval = 'conf')
cisplatin_bone_least_auc <- auc(cisplatin_test_bone$least_sensitive, new_ic50)

cisplatin_test_scaled_breast <- cisplatin_rna_seq_test_scaled[cisplatin_breast_lines, ]
cisplatin_test_breast <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_breast), s = 'lambda.min', interval = 'conf')
cisplatin_breast_least_auc <- auc(cisplatin_test_breast$least_sensitive, new_ic50)

cisplatin_test_scaled_digestive_system <- cisplatin_rna_seq_test_scaled[cisplatin_digestive_system_lines, ]
cisplatin_test_digestive_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_digestive_system), s = 'lambda.min', interval = 'conf')
cisplatin_digestive_system_least_auc <- auc(cisplatin_test_digestive_system$least_sensitive, new_ic50)

cisplatin_test_scaled_kidney <- cisplatin_rna_seq_test_scaled[cisplatin_kidney_lines, ]
cisplatin_test_kidney <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_kidney), s = 'lambda.min', interval = 'conf')
cisplatin_kidney_least_auc <- auc(cisplatin_test_kidney$least_sensitive, new_ic50)

cisplatin_test_scaled_large_intestine <- cisplatin_rna_seq_test_scaled[cisplatin_large_intestine_lines, ]
cisplatin_test_large_intestine <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_large_intestine), s = 'lambda.min', interval = 'conf')
cisplatin_large_intestine_least_auc <- auc(cisplatin_test_large_intestine$least_sensitive, new_ic50)

cisplatin_test_scaled_lung <- cisplatin_rna_seq_test_scaled[cisplatin_lung_lines, ]
cisplatin_test_lung <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung), s = 'lambda.min', interval = 'conf')
cisplatin_lung_least_auc <- auc(cisplatin_test_lung$least_sensitive, new_ic50)

cisplatin_test_scaled_lung_NSCLC <- cisplatin_rna_seq_test_scaled[cisplatin_lung_NSCLC_lines, ]
cisplatin_test_lung_NSCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung_NSCLC), s = 'lambda.min', interval = 'conf')
cisplatin_lung_NSCLC_least_auc <- auc(cisplatin_test_lung_NSCLC$least_sensitive, new_ic50)

cisplatin_test_scaled_lung_SCLC <- cisplatin_rna_seq_test_scaled[cisplatin_lung_SCLC_lines, ]
cisplatin_test_lung_SCLC <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_lung_SCLC), s = 'lambda.min', interval = 'conf')
cisplatin_lung_SCLC_least_auc <- auc(cisplatin_test_lung_SCLC$least_sensitive, new_ic50)

cisplatin_test_scaled_nervous_system <- cisplatin_rna_seq_test_scaled[cisplatin_nervous_system_lines, ]
cisplatin_test_nervous_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_nervous_system), s = 'lambda.min', interval = 'conf')
cisplatin_nervous_system_least_auc <- auc(cisplatin_test_nervous_system$least_sensitive, new_ic50)

cisplatin_test_scaled_neuroblastoma <- cisplatin_rna_seq_test_scaled[cisplatin_neuroblastoma_lines, ]
cisplatin_test_neuroblastoma <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_neuroblastoma), s = 'lambda.min', interval = 'conf')
cisplatin_neuroblastoma_least_auc <- auc(cisplatin_test_neuroblastoma$least_sensitive, new_ic50)

cisplatin_test_scaled_pancreas <- cisplatin_rna_seq_test_scaled[cisplatin_pancreas_lines, ]
cisplatin_test_pancreas <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_pancreas), s = 'lambda.min', interval = 'conf')
cisplatin_pancreas_least_auc <- auc(cisplatin_test_pancreas$least_sensitive, new_ic50)

cisplatin_test_scaled_skin <- cisplatin_rna_seq_test_scaled[cisplatin_skin_lines, ]
cisplatin_test_skin <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_skin), s = 'lambda.min', interval = 'conf')
cisplatin_skin_least_auc <- auc(cisplatin_test_skin$least_sensitive, new_ic50)

cisplatin_test_scaled_soft_tissue <- cisplatin_rna_seq_test_scaled[cisplatin_soft_tissue_lines, ]
cisplatin_test_soft_tissue <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_soft_tissue), s = 'lambda.min', interval = 'conf')
cisplatin_soft_tissue_least_auc <- auc(cisplatin_test_soft_tissue$least_sensitive, new_ic50)

cisplatin_test_scaled_thyroid <- cisplatin_rna_seq_test_scaled[cisplatin_thyroid_lines, ]
cisplatin_test_thyroid <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_thyroid), s = 'lambda.min', interval = 'conf')
cisplatin_thyroid_least_auc <- auc(cisplatin_test_thyroid$least_sensitive, new_ic50)

cisplatin_test_scaled_urogenital_system <- cisplatin_rna_seq_test_scaled[cisplatin_urogenital_system_lines, ]
cisplatin_test_urogenital_system <- cisplatin_test[which(cisplatin_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_test_scaled_urogenital_system), s = 'lambda.min', interval = 'conf')
cisplatin_urogenital_system_least_auc <- auc(cisplatin_test_urogenital_system$least_sensitive, new_ic50)

cisplatin_least_auc <- c(cisplatin_aero_dig_tract_least_auc, cisplatin_bone_least_auc, 
                        cisplatin_breast_least_auc, cisplatin_digestive_system_least_auc, 
                        cisplatin_kidney_least_auc, cisplatin_large_intestine_least_auc, 
                        cisplatin_lung_least_auc, cisplatin_lung_NSCLC_least_auc, cisplatin_lung_SCLC_least_auc, 
                        cisplatin_nervous_system_least_auc, cisplatin_neuroblastoma_least_auc, 
                        cisplatin_pancreas_least_auc, cisplatin_skin_least_auc, 
                        cisplatin_soft_tissue_least_auc, cisplatin_thyroid_least_auc, 
                        cisplatin_urogenital_system_least_auc)

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

etoposide_test_scaled_aero_dig_tract <- etoposide_rna_seq_test_scaled[etoposide_aero_dig_tract_lines, ]
etoposide_test_aero_dig_tract <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_aero_dig_tract), s = 'lambda.min', interval = 'conf')
etoposide_aero_dig_tract_most_auc <- auc(etoposide_test_aero_dig_tract$most_sensitive, new_ic50)

etoposide_test_scaled_bone <- etoposide_rna_seq_test_scaled[etoposide_bone_lines, ]
etoposide_test_bone <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_bone), s = 'lambda.min', interval = 'conf')
etoposide_bone_most_auc <- auc(etoposide_test_bone$most_sensitive, new_ic50)

etoposide_test_scaled_breast <- etoposide_rna_seq_test_scaled[etoposide_breast_lines, ]
etoposide_test_breast <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_breast), s = 'lambda.min', interval = 'conf')
etoposide_breast_most_auc <- auc(etoposide_test_breast$most_sensitive, new_ic50)

etoposide_test_scaled_digestive_system <- etoposide_rna_seq_test_scaled[etoposide_digestive_system_lines, ]
etoposide_test_digestive_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_digestive_system), s = 'lambda.min', interval = 'conf')
etoposide_digestive_system_most_auc <- auc(etoposide_test_digestive_system$most_sensitive, new_ic50)

etoposide_test_scaled_kidney <- etoposide_rna_seq_test_scaled[etoposide_kidney_lines, ]
etoposide_test_kidney <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_kidney), s = 'lambda.min', interval = 'conf')
etoposide_kidney_most_auc <- auc(etoposide_test_kidney$most_sensitive, new_ic50)

etoposide_test_scaled_large_intestine <- etoposide_rna_seq_test_scaled[etoposide_large_intestine_lines, ]
etoposide_test_large_intestine <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_large_intestine), s = 'lambda.min', interval = 'conf')
etoposide_large_intestine_most_auc <- auc(etoposide_test_large_intestine$most_sensitive, new_ic50)

etoposide_test_scaled_lung <- etoposide_rna_seq_test_scaled[etoposide_lung_lines, ]
etoposide_test_lung <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung), s = 'lambda.min', interval = 'conf')
etoposide_lung_most_auc <- auc(etoposide_test_lung$most_sensitive, new_ic50)

etoposide_test_scaled_lung_NSCLC <- etoposide_rna_seq_test_scaled[etoposide_lung_NSCLC_lines, ]
etoposide_test_lung_NSCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung_NSCLC), s = 'lambda.min', interval = 'conf')
etoposide_lung_NSCLC_most_auc <- auc(etoposide_test_lung_NSCLC$most_sensitive, new_ic50)

etoposide_test_scaled_lung_SCLC <- etoposide_rna_seq_test_scaled[etoposide_lung_SCLC_lines, ]
etoposide_test_lung_SCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung_SCLC), s = 'lambda.min', interval = 'conf')
etoposide_lung_SCLC_most_auc <- auc(etoposide_test_lung_SCLC$most_sensitive, new_ic50)

etoposide_test_scaled_nervous_system <- etoposide_rna_seq_test_scaled[etoposide_nervous_system_lines, ]
etoposide_test_nervous_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_nervous_system), s = 'lambda.min', interval = 'conf')
etoposide_nervous_system_most_auc <- auc(etoposide_test_nervous_system$most_sensitive, new_ic50)

etoposide_test_scaled_neuroblastoma <- etoposide_rna_seq_test_scaled[etoposide_neuroblastoma_lines, ]
etoposide_test_neuroblastoma <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_neuroblastoma), s = 'lambda.min', interval = 'conf')
etoposide_neuroblastoma_most_auc <- auc(etoposide_test_neuroblastoma$most_sensitive, new_ic50)

etoposide_test_scaled_pancreas <- etoposide_rna_seq_test_scaled[etoposide_pancreas_lines, ]
etoposide_test_pancreas <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_pancreas), s = 'lambda.min', interval = 'conf')
etoposide_pancreas_most_auc <- auc(etoposide_test_pancreas$most_sensitive, new_ic50)

etoposide_test_scaled_skin <- etoposide_rna_seq_test_scaled[etoposide_skin_lines, ]
etoposide_test_skin <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_skin), s = 'lambda.min', interval = 'conf')
etoposide_skin_most_auc <- auc(etoposide_test_skin$most_sensitive, new_ic50)

etoposide_test_scaled_soft_tissue <- etoposide_rna_seq_test_scaled[etoposide_soft_tissue_lines, ]
etoposide_test_soft_tissue <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_soft_tissue), s = 'lambda.min', interval = 'conf')
etoposide_soft_tissue_most_auc <- auc(etoposide_test_soft_tissue$most_sensitive, new_ic50)

etoposide_test_scaled_thyroid <- etoposide_rna_seq_test_scaled[etoposide_thyroid_lines, ]
etoposide_test_thyroid <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_thyroid), s = 'lambda.min', interval = 'conf')
etoposide_thyroid_most_auc <- auc(etoposide_test_thyroid$most_sensitive, new_ic50)

etoposide_test_scaled_urogenital_system <- etoposide_rna_seq_test_scaled[etoposide_urogenital_system_lines, ]
etoposide_test_urogenital_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_urogenital_system), s = 'lambda.min', interval = 'conf')
etoposide_urogenital_system_most_auc <- auc(etoposide_test_urogenital_system$most_sensitive, new_ic50)

etoposide_most_auc <- c(etoposide_aero_dig_tract_most_auc, etoposide_bone_most_auc, 
                        etoposide_breast_most_auc, etoposide_digestive_system_most_auc, 
                        etoposide_kidney_most_auc, etoposide_large_intestine_most_auc, 
                        etoposide_lung_most_auc, etoposide_lung_NSCLC_most_auc, etoposide_lung_SCLC_most_auc, 
                        etoposide_nervous_system_most_auc, etoposide_neuroblastoma_most_auc, 
                        etoposide_pancreas_most_auc, etoposide_skin_most_auc, 
                        etoposide_soft_tissue_most_auc, etoposide_thyroid_most_auc, 
                        etoposide_urogenital_system_most_auc)

etoposide_test_scaled_aero_dig_tract <- etoposide_rna_seq_test_scaled[etoposide_aero_dig_tract_lines, ]
etoposide_test_aero_dig_tract <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_aero_dig_tract), s = 'lambda.min', interval = 'conf')
etoposide_aero_dig_tract_least_auc <- auc(etoposide_test_aero_dig_tract$least_sensitive, new_ic50)

etoposide_test_scaled_bone <- etoposide_rna_seq_test_scaled[etoposide_bone_lines, ]
etoposide_test_bone <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_bone), s = 'lambda.min', interval = 'conf')
etoposide_bone_least_auc <- auc(etoposide_test_bone$least_sensitive, new_ic50)

etoposide_test_scaled_breast <- etoposide_rna_seq_test_scaled[etoposide_breast_lines, ]
etoposide_test_breast <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_breast), s = 'lambda.min', interval = 'conf')
etoposide_breast_least_auc <- auc(etoposide_test_breast$least_sensitive, new_ic50)

etoposide_test_scaled_digestive_system <- etoposide_rna_seq_test_scaled[etoposide_digestive_system_lines, ]
etoposide_test_digestive_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_digestive_system), s = 'lambda.min', interval = 'conf')
etoposide_digestive_system_least_auc <- auc(etoposide_test_digestive_system$least_sensitive, new_ic50)

etoposide_test_scaled_kidney <- etoposide_rna_seq_test_scaled[etoposide_kidney_lines, ]
etoposide_test_kidney <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_kidney), s = 'lambda.min', interval = 'conf')
etoposide_kidney_least_auc <- auc(etoposide_test_kidney$least_sensitive, new_ic50)

etoposide_test_scaled_large_intestine <- etoposide_rna_seq_test_scaled[etoposide_large_intestine_lines, ]
etoposide_test_large_intestine <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_large_intestine), s = 'lambda.min', interval = 'conf')
etoposide_large_intestine_least_auc <- auc(etoposide_test_large_intestine$least_sensitive, new_ic50)

etoposide_test_scaled_lung <- etoposide_rna_seq_test_scaled[etoposide_lung_lines, ]
etoposide_test_lung <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung), s = 'lambda.min', interval = 'conf')
etoposide_lung_least_auc <- auc(etoposide_test_lung$least_sensitive, new_ic50)

etoposide_test_scaled_lung_NSCLC <- etoposide_rna_seq_test_scaled[etoposide_lung_NSCLC_lines, ]
etoposide_test_lung_NSCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung_NSCLC), s = 'lambda.min', interval = 'conf')
etoposide_lung_NSCLC_least_auc <- auc(etoposide_test_lung_NSCLC$least_sensitive, new_ic50)

etoposide_test_scaled_lung_SCLC <- etoposide_rna_seq_test_scaled[etoposide_lung_SCLC_lines, ]
etoposide_test_lung_SCLC <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung_SCLC), s = 'lambda.min', interval = 'conf')
etoposide_lung_SCLC_least_auc <- auc(etoposide_test_lung_SCLC$least_sensitive, new_ic50)

etoposide_test_scaled_nervous_system <- etoposide_rna_seq_test_scaled[etoposide_nervous_system_lines, ]
etoposide_test_nervous_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_nervous_system), s = 'lambda.min', interval = 'conf')
etoposide_nervous_system_least_auc <- auc(etoposide_test_nervous_system$least_sensitive, new_ic50)

etoposide_test_scaled_neuroblastoma <- etoposide_rna_seq_test_scaled[etoposide_neuroblastoma_lines, ]
etoposide_test_neuroblastoma <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_neuroblastoma), s = 'lambda.min', interval = 'conf')
etoposide_neuroblastoma_least_auc <- auc(etoposide_test_neuroblastoma$least_sensitive, new_ic50)

etoposide_test_scaled_pancreas <- etoposide_rna_seq_test_scaled[etoposide_pancreas_lines, ]
etoposide_test_pancreas <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_pancreas), s = 'lambda.min', interval = 'conf')
etoposide_pancreas_least_auc <- auc(etoposide_test_pancreas$least_sensitive, new_ic50)

etoposide_test_scaled_skin <- etoposide_rna_seq_test_scaled[etoposide_skin_lines, ]
etoposide_test_skin <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_skin), s = 'lambda.min', interval = 'conf')
etoposide_skin_least_auc <- auc(etoposide_test_skin$least_sensitive, new_ic50)

etoposide_test_scaled_soft_tissue <- etoposide_rna_seq_test_scaled[etoposide_soft_tissue_lines, ]
etoposide_test_soft_tissue <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_soft_tissue), s = 'lambda.min', interval = 'conf')
etoposide_soft_tissue_least_auc <- auc(etoposide_test_soft_tissue$least_sensitive, new_ic50)

etoposide_test_scaled_thyroid <- etoposide_rna_seq_test_scaled[etoposide_thyroid_lines, ]
etoposide_test_thyroid <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_thyroid), s = 'lambda.min', interval = 'conf')
etoposide_thyroid_least_auc <- auc(etoposide_test_thyroid$least_sensitive, new_ic50)

etoposide_test_scaled_urogenital_system <- etoposide_rna_seq_test_scaled[etoposide_urogenital_system_lines, ]
etoposide_test_urogenital_system <- etoposide_test[which(etoposide_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_urogenital_system), s = 'lambda.min', interval = 'conf')
etoposide_urogenital_system_least_auc <- auc(etoposide_test_urogenital_system$least_sensitive, new_ic50)

etoposide_least_auc <- c(etoposide_aero_dig_tract_least_auc, etoposide_bone_least_auc, 
                         etoposide_breast_least_auc, etoposide_digestive_system_least_auc, 
                         etoposide_kidney_least_auc, etoposide_large_intestine_least_auc, 
                         etoposide_lung_least_auc, etoposide_lung_NSCLC_least_auc, etoposide_lung_SCLC_least_auc, 
                         etoposide_nervous_system_least_auc, etoposide_neuroblastoma_least_auc, 
                         etoposide_pancreas_least_auc, etoposide_skin_least_auc, 
                         etoposide_soft_tissue_least_auc, etoposide_thyroid_least_auc, 
                         etoposide_urogenital_system_least_auc)

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

gemcitabine_test_scaled_aero_dig_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_aero_dig_tract_lines, ]
gemcitabine_test_aero_dig_tract <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_aero_dig_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_aero_dig_tract_most_auc <- auc(gemcitabine_test_aero_dig_tract$most_sensitive, new_ic50)

gemcitabine_test_scaled_bone <- gemcitabine_rna_seq_test_scaled[gemcitabine_bone_lines, ]
gemcitabine_test_bone <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
gemcitabine_bone_most_auc <- auc(gemcitabine_test_bone$most_sensitive, new_ic50)

gemcitabine_test_scaled_breast <- gemcitabine_rna_seq_test_scaled[gemcitabine_breast_lines, ]
gemcitabine_test_breast <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
gemcitabine_breast_most_auc <- auc(gemcitabine_test_breast$most_sensitive, new_ic50)

gemcitabine_test_scaled_digestive_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_digestive_system_lines, ]
gemcitabine_test_digestive_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_digestive_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_digestive_system_most_auc <- auc(gemcitabine_test_digestive_system$most_sensitive, new_ic50)

gemcitabine_test_scaled_kidney <- gemcitabine_rna_seq_test_scaled[gemcitabine_kidney_lines, ]
gemcitabine_test_kidney <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
gemcitabine_kidney_most_auc <- auc(gemcitabine_test_kidney$most_sensitive, new_ic50)

gemcitabine_test_scaled_large_intestine <- gemcitabine_rna_seq_test_scaled[gemcitabine_large_intestine_lines, ]
gemcitabine_test_large_intestine <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
gemcitabine_large_intestine_most_auc <- auc(gemcitabine_test_large_intestine$most_sensitive, new_ic50)

gemcitabine_test_scaled_lung <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_lines, ]
gemcitabine_test_lung <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_most_auc <- auc(gemcitabine_test_lung$most_sensitive, new_ic50)

gemcitabine_test_scaled_lung_NSCLC <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_NSCLC_lines, ]
gemcitabine_test_lung_NSCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung_NSCLC), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_NSCLC_most_auc <- auc(gemcitabine_test_lung_NSCLC$most_sensitive, new_ic50)

gemcitabine_test_scaled_lung_SCLC <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_SCLC_lines, ]
gemcitabine_test_lung_SCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung_SCLC), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_SCLC_most_auc <- auc(gemcitabine_test_lung_SCLC$most_sensitive, new_ic50)

gemcitabine_test_scaled_nervous_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_nervous_system_lines, ]
gemcitabine_test_nervous_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_nervous_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_nervous_system_most_auc <- auc(gemcitabine_test_nervous_system$most_sensitive, new_ic50)

gemcitabine_test_scaled_neuroblastoma <- gemcitabine_rna_seq_test_scaled[gemcitabine_neuroblastoma_lines, ]
gemcitabine_test_neuroblastoma <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_neuroblastoma), s = 'lambda.1se', interval = 'conf')
gemcitabine_neuroblastoma_most_auc <- auc(gemcitabine_test_neuroblastoma$most_sensitive, new_ic50)

gemcitabine_test_scaled_pancreas <- gemcitabine_rna_seq_test_scaled[gemcitabine_pancreas_lines, ]
gemcitabine_test_pancreas <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
gemcitabine_pancreas_most_auc <- auc(gemcitabine_test_pancreas$most_sensitive, new_ic50)

gemcitabine_test_scaled_skin <- gemcitabine_rna_seq_test_scaled[gemcitabine_skin_lines, ]
gemcitabine_test_skin <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
gemcitabine_skin_most_auc <- auc(gemcitabine_test_skin$most_sensitive, new_ic50)

gemcitabine_test_scaled_soft_tissue <- gemcitabine_rna_seq_test_scaled[gemcitabine_soft_tissue_lines, ]
gemcitabine_test_soft_tissue <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
gemcitabine_soft_tissue_most_auc <- auc(gemcitabine_test_soft_tissue$most_sensitive, new_ic50)

gemcitabine_test_scaled_thyroid <- gemcitabine_rna_seq_test_scaled[gemcitabine_thyroid_lines, ]
gemcitabine_test_thyroid <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
gemcitabine_thyroid_most_auc <- auc(gemcitabine_test_thyroid$most_sensitive, new_ic50)

gemcitabine_test_scaled_urogenital_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_urogenital_system_lines, ]
gemcitabine_test_urogenital_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_urogenital_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_urogenital_system_most_auc <- auc(gemcitabine_test_urogenital_system$most_sensitive, new_ic50)

gemcitabine_most_auc <- c(gemcitabine_aero_dig_tract_most_auc, gemcitabine_bone_most_auc, 
                        gemcitabine_breast_most_auc, gemcitabine_digestive_system_most_auc, 
                        gemcitabine_kidney_most_auc, gemcitabine_large_intestine_most_auc, 
                        gemcitabine_lung_most_auc, gemcitabine_lung_NSCLC_most_auc, gemcitabine_lung_SCLC_most_auc, 
                        gemcitabine_nervous_system_most_auc, gemcitabine_neuroblastoma_most_auc, 
                        gemcitabine_pancreas_most_auc, gemcitabine_skin_most_auc, 
                        gemcitabine_soft_tissue_most_auc, gemcitabine_thyroid_most_auc, 
                        gemcitabine_urogenital_system_most_auc)

gemcitabine_test_scaled_aero_dig_tract <- gemcitabine_rna_seq_test_scaled[gemcitabine_aero_dig_tract_lines, ]
gemcitabine_test_aero_dig_tract <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_aero_dig_tract), s = 'lambda.1se', interval = 'conf')
gemcitabine_aero_dig_tract_least_auc <- auc(gemcitabine_test_aero_dig_tract$least_sensitive, new_ic50)

gemcitabine_test_scaled_bone <- gemcitabine_rna_seq_test_scaled[gemcitabine_bone_lines, ]
gemcitabine_test_bone <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
gemcitabine_bone_least_auc <- auc(gemcitabine_test_bone$least_sensitive, new_ic50)

gemcitabine_test_scaled_breast <- gemcitabine_rna_seq_test_scaled[gemcitabine_breast_lines, ]
gemcitabine_test_breast <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
gemcitabine_breast_least_auc <- auc(gemcitabine_test_breast$least_sensitive, new_ic50)

gemcitabine_test_scaled_digestive_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_digestive_system_lines, ]
gemcitabine_test_digestive_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_digestive_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_digestive_system_least_auc <- auc(gemcitabine_test_digestive_system$least_sensitive, new_ic50)

gemcitabine_test_scaled_kidney <- gemcitabine_rna_seq_test_scaled[gemcitabine_kidney_lines, ]
gemcitabine_test_kidney <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
gemcitabine_kidney_least_auc <- auc(gemcitabine_test_kidney$least_sensitive, new_ic50)

gemcitabine_test_scaled_large_intestine <- gemcitabine_rna_seq_test_scaled[gemcitabine_large_intestine_lines, ]
gemcitabine_test_large_intestine <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
gemcitabine_large_intestine_least_auc <- auc(gemcitabine_test_large_intestine$least_sensitive, new_ic50)

gemcitabine_test_scaled_lung <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_lines, ]
gemcitabine_test_lung <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_least_auc <- auc(gemcitabine_test_lung$least_sensitive, new_ic50)

gemcitabine_test_scaled_lung_NSCLC <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_NSCLC_lines, ]
gemcitabine_test_lung_NSCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung_NSCLC), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_NSCLC_least_auc <- auc(gemcitabine_test_lung_NSCLC$least_sensitive, new_ic50)

gemcitabine_test_scaled_lung_SCLC <- gemcitabine_rna_seq_test_scaled[gemcitabine_lung_SCLC_lines, ]
gemcitabine_test_lung_SCLC <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_lung_SCLC), s = 'lambda.1se', interval = 'conf')
gemcitabine_lung_SCLC_least_auc <- auc(gemcitabine_test_lung_SCLC$least_sensitive, new_ic50)

gemcitabine_test_scaled_nervous_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_nervous_system_lines, ]
gemcitabine_test_nervous_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_nervous_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_nervous_system_least_auc <- auc(gemcitabine_test_nervous_system$least_sensitive, new_ic50)

gemcitabine_test_scaled_neuroblastoma <- gemcitabine_rna_seq_test_scaled[gemcitabine_neuroblastoma_lines, ]
gemcitabine_test_neuroblastoma <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_neuroblastoma), s = 'lambda.1se', interval = 'conf')
gemcitabine_neuroblastoma_least_auc <- auc(gemcitabine_test_neuroblastoma$least_sensitive, new_ic50)

gemcitabine_test_scaled_pancreas <- gemcitabine_rna_seq_test_scaled[gemcitabine_pancreas_lines, ]
gemcitabine_test_pancreas <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
gemcitabine_pancreas_least_auc <- auc(gemcitabine_test_pancreas$least_sensitive, new_ic50)

gemcitabine_test_scaled_skin <- gemcitabine_rna_seq_test_scaled[gemcitabine_skin_lines, ]
gemcitabine_test_skin <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
gemcitabine_skin_least_auc <- auc(gemcitabine_test_skin$least_sensitive, new_ic50)

gemcitabine_test_scaled_soft_tissue <- gemcitabine_rna_seq_test_scaled[gemcitabine_soft_tissue_lines, ]
gemcitabine_test_soft_tissue <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
gemcitabine_soft_tissue_least_auc <- auc(gemcitabine_test_soft_tissue$least_sensitive, new_ic50)

gemcitabine_test_scaled_thyroid <- gemcitabine_rna_seq_test_scaled[gemcitabine_thyroid_lines, ]
gemcitabine_test_thyroid <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
gemcitabine_thyroid_least_auc <- auc(gemcitabine_test_thyroid$least_sensitive, new_ic50)

gemcitabine_test_scaled_urogenital_system <- gemcitabine_rna_seq_test_scaled[gemcitabine_urogenital_system_lines, ]
gemcitabine_test_urogenital_system <- gemcitabine_test[which(gemcitabine_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_test_scaled_urogenital_system), s = 'lambda.1se', interval = 'conf')
gemcitabine_urogenital_system_least_auc <- auc(gemcitabine_test_urogenital_system$least_sensitive, new_ic50)

gemcitabine_least_auc <- c(gemcitabine_aero_dig_tract_least_auc, gemcitabine_bone_least_auc, 
                         gemcitabine_breast_least_auc, gemcitabine_digestive_system_least_auc, 
                         gemcitabine_kidney_least_auc, gemcitabine_large_intestine_least_auc, 
                         gemcitabine_lung_least_auc, gemcitabine_lung_NSCLC_least_auc, gemcitabine_lung_SCLC_least_auc, 
                         gemcitabine_nervous_system_least_auc, gemcitabine_neuroblastoma_least_auc, 
                         gemcitabine_pancreas_least_auc, gemcitabine_skin_least_auc, 
                         gemcitabine_soft_tissue_least_auc, gemcitabine_thyroid_least_auc, 
                         gemcitabine_urogenital_system_least_auc)

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

methotrexate_test_scaled_aero_dig_tract <- methotrexate_rna_seq_test_scaled[methotrexate_aero_dig_tract_lines, ]
methotrexate_test_aero_dig_tract <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_aero_dig_tract), s = 'lambda.min', interval = 'conf')
methotrexate_aero_dig_tract_most_auc <- auc(methotrexate_test_aero_dig_tract$most_sensitive, new_ic50)

methotrexate_test_scaled_bone <- methotrexate_rna_seq_test_scaled[methotrexate_bone_lines, ]
methotrexate_test_bone <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_bone), s = 'lambda.min', interval = 'conf')
methotrexate_bone_most_auc <- auc(methotrexate_test_bone$most_sensitive, new_ic50)

methotrexate_test_scaled_breast <- methotrexate_rna_seq_test_scaled[methotrexate_breast_lines, ]
methotrexate_test_breast <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_breast), s = 'lambda.min', interval = 'conf')
methotrexate_breast_most_auc <- auc(methotrexate_test_breast$most_sensitive, new_ic50)

methotrexate_test_scaled_digestive_system <- methotrexate_rna_seq_test_scaled[methotrexate_digestive_system_lines, ]
methotrexate_test_digestive_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_digestive_system), s = 'lambda.min', interval = 'conf')
methotrexate_digestive_system_most_auc <- auc(methotrexate_test_digestive_system$most_sensitive, new_ic50)

methotrexate_test_scaled_kidney <- methotrexate_rna_seq_test_scaled[methotrexate_kidney_lines, ]
methotrexate_test_kidney <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_kidney), s = 'lambda.min', interval = 'conf')
methotrexate_kidney_most_auc <- auc(methotrexate_test_kidney$most_sensitive, new_ic50)

methotrexate_test_scaled_large_intestine <- methotrexate_rna_seq_test_scaled[methotrexate_large_intestine_lines, ]
methotrexate_test_large_intestine <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_large_intestine), s = 'lambda.min', interval = 'conf')
methotrexate_large_intestine_most_auc <- auc(methotrexate_test_large_intestine$most_sensitive, new_ic50)

methotrexate_test_scaled_lung <- methotrexate_rna_seq_test_scaled[methotrexate_lung_lines, ]
methotrexate_test_lung <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung), s = 'lambda.min', interval = 'conf')
methotrexate_lung_most_auc <- auc(methotrexate_test_lung$most_sensitive, new_ic50)

methotrexate_test_scaled_lung_NSCLC <- methotrexate_rna_seq_test_scaled[methotrexate_lung_NSCLC_lines, ]
methotrexate_test_lung_NSCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung_NSCLC), s = 'lambda.min', interval = 'conf')
methotrexate_lung_NSCLC_most_auc <- auc(methotrexate_test_lung_NSCLC$most_sensitive, new_ic50)

methotrexate_test_scaled_lung_SCLC <- methotrexate_rna_seq_test_scaled[methotrexate_lung_SCLC_lines, ]
methotrexate_test_lung_SCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung_SCLC), s = 'lambda.min', interval = 'conf')
methotrexate_lung_SCLC_most_auc <- auc(methotrexate_test_lung_SCLC$most_sensitive, new_ic50)

methotrexate_test_scaled_nervous_system <- methotrexate_rna_seq_test_scaled[methotrexate_nervous_system_lines, ]
methotrexate_test_nervous_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_nervous_system), s = 'lambda.min', interval = 'conf')
methotrexate_nervous_system_most_auc <- auc(methotrexate_test_nervous_system$most_sensitive, new_ic50)

methotrexate_test_scaled_neuroblastoma <- methotrexate_rna_seq_test_scaled[methotrexate_neuroblastoma_lines, ]
methotrexate_test_neuroblastoma <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_neuroblastoma), s = 'lambda.min', interval = 'conf')
methotrexate_neuroblastoma_most_auc <- auc(methotrexate_test_neuroblastoma$most_sensitive, new_ic50)

methotrexate_test_scaled_pancreas <- methotrexate_rna_seq_test_scaled[methotrexate_pancreas_lines, ]
methotrexate_test_pancreas <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pancreas), s = 'lambda.min', interval = 'conf')
methotrexate_pancreas_most_auc <- auc(methotrexate_test_pancreas$most_sensitive, new_ic50)

methotrexate_test_scaled_skin <- methotrexate_rna_seq_test_scaled[methotrexate_skin_lines, ]
methotrexate_test_skin <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_skin), s = 'lambda.min', interval = 'conf')
methotrexate_skin_most_auc <- auc(methotrexate_test_skin$most_sensitive, new_ic50)

methotrexate_test_scaled_soft_tissue <- methotrexate_rna_seq_test_scaled[methotrexate_soft_tissue_lines, ]
methotrexate_test_soft_tissue <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_soft_tissue), s = 'lambda.min', interval = 'conf')
methotrexate_soft_tissue_most_auc <- auc(methotrexate_test_soft_tissue$most_sensitive, new_ic50)

methotrexate_test_scaled_thyroid <- methotrexate_rna_seq_test_scaled[methotrexate_thyroid_lines, ]
methotrexate_test_thyroid <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_thyroid), s = 'lambda.min', interval = 'conf')
methotrexate_thyroid_most_auc <- auc(methotrexate_test_thyroid$most_sensitive, new_ic50)

methotrexate_test_scaled_urogenital_system <- methotrexate_rna_seq_test_scaled[methotrexate_urogenital_system_lines, ]
methotrexate_test_urogenital_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_urogenital_system), s = 'lambda.min', interval = 'conf')
methotrexate_urogenital_system_most_auc <- auc(methotrexate_test_urogenital_system$most_sensitive, new_ic50)

methotrexate_most_auc <- c(methotrexate_aero_dig_tract_most_auc, methotrexate_bone_most_auc, 
                        methotrexate_breast_most_auc, methotrexate_digestive_system_most_auc, 
                        methotrexate_kidney_most_auc, methotrexate_large_intestine_most_auc, 
                        methotrexate_lung_most_auc, methotrexate_lung_NSCLC_most_auc, methotrexate_lung_SCLC_most_auc, 
                        methotrexate_nervous_system_most_auc, methotrexate_neuroblastoma_most_auc, 
                        methotrexate_pancreas_most_auc, methotrexate_skin_most_auc, 
                        methotrexate_soft_tissue_most_auc, methotrexate_thyroid_most_auc, 
                        methotrexate_urogenital_system_most_auc)

methotrexate_test_scaled_aero_dig_tract <- methotrexate_rna_seq_test_scaled[methotrexate_aero_dig_tract_lines, ]
methotrexate_test_aero_dig_tract <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'aero_dig_tract'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_aero_dig_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_aero_dig_tract_least_auc <- auc(methotrexate_test_aero_dig_tract$least_sensitive, new_ic50)

methotrexate_test_scaled_bone <- methotrexate_rna_seq_test_scaled[methotrexate_bone_lines, ]
methotrexate_test_bone <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'bone'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
methotrexate_bone_least_auc <- auc(methotrexate_test_bone$least_sensitive, new_ic50)

methotrexate_test_scaled_breast <- methotrexate_rna_seq_test_scaled[methotrexate_breast_lines, ]
methotrexate_test_breast <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'breast'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
methotrexate_breast_least_auc <- auc(methotrexate_test_breast$least_sensitive, new_ic50)

methotrexate_test_scaled_digestive_system <- methotrexate_rna_seq_test_scaled[methotrexate_digestive_system_lines, ]
methotrexate_test_digestive_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'digestive_system'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_digestive_system), s = 'lambda.1se', interval = 'conf')
methotrexate_digestive_system_least_auc <- auc(methotrexate_test_digestive_system$least_sensitive, new_ic50)

methotrexate_test_scaled_kidney <- methotrexate_rna_seq_test_scaled[methotrexate_kidney_lines, ]
methotrexate_test_kidney <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'kidney'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
methotrexate_kidney_least_auc <- auc(methotrexate_test_kidney$least_sensitive, new_ic50)

methotrexate_test_scaled_large_intestine <- methotrexate_rna_seq_test_scaled[methotrexate_large_intestine_lines, ]
methotrexate_test_large_intestine <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'large_intestine'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
methotrexate_large_intestine_least_auc <- auc(methotrexate_test_large_intestine$least_sensitive, new_ic50)

methotrexate_test_scaled_lung <- methotrexate_rna_seq_test_scaled[methotrexate_lung_lines, ]
methotrexate_test_lung <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
methotrexate_lung_least_auc <- auc(methotrexate_test_lung$least_sensitive, new_ic50)

methotrexate_test_scaled_lung_NSCLC <- methotrexate_rna_seq_test_scaled[methotrexate_lung_NSCLC_lines, ]
methotrexate_test_lung_NSCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_NSCLC'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung_NSCLC), s = 'lambda.1se', interval = 'conf')
methotrexate_lung_NSCLC_least_auc <- auc(methotrexate_test_lung_NSCLC$least_sensitive, new_ic50)

methotrexate_test_scaled_lung_SCLC <- methotrexate_rna_seq_test_scaled[methotrexate_lung_SCLC_lines, ]
methotrexate_test_lung_SCLC <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'lung_SCLC'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung_SCLC), s = 'lambda.1se', interval = 'conf')
methotrexate_lung_SCLC_least_auc <- auc(methotrexate_test_lung_SCLC$least_sensitive, new_ic50)

methotrexate_test_scaled_nervous_system <- methotrexate_rna_seq_test_scaled[methotrexate_nervous_system_lines, ]
methotrexate_test_nervous_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'nervous_system'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_nervous_system), s = 'lambda.1se', interval = 'conf')
methotrexate_nervous_system_least_auc <- auc(methotrexate_test_nervous_system$least_sensitive, new_ic50)

methotrexate_test_scaled_neuroblastoma <- methotrexate_rna_seq_test_scaled[methotrexate_neuroblastoma_lines, ]
methotrexate_test_neuroblastoma <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'neuroblastoma'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_neuroblastoma), s = 'lambda.1se', interval = 'conf')
methotrexate_neuroblastoma_least_auc <- auc(methotrexate_test_neuroblastoma$least_sensitive, new_ic50)

methotrexate_test_scaled_pancreas <- methotrexate_rna_seq_test_scaled[methotrexate_pancreas_lines, ]
methotrexate_test_pancreas <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'pancreas'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
methotrexate_pancreas_least_auc <- auc(methotrexate_test_pancreas$least_sensitive, new_ic50)

methotrexate_test_scaled_skin <- methotrexate_rna_seq_test_scaled[methotrexate_skin_lines, ]
methotrexate_test_skin <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'skin'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
methotrexate_skin_least_auc <- auc(methotrexate_test_skin$least_sensitive, new_ic50)

methotrexate_test_scaled_soft_tissue <- methotrexate_rna_seq_test_scaled[methotrexate_soft_tissue_lines, ]
methotrexate_test_soft_tissue <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'soft_tissue'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
methotrexate_soft_tissue_least_auc <- auc(methotrexate_test_soft_tissue$least_sensitive, new_ic50)

methotrexate_test_scaled_thyroid <- methotrexate_rna_seq_test_scaled[methotrexate_thyroid_lines, ]
methotrexate_test_thyroid <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'thyroid'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
methotrexate_thyroid_least_auc <- auc(methotrexate_test_thyroid$least_sensitive, new_ic50)

methotrexate_test_scaled_urogenital_system <- methotrexate_rna_seq_test_scaled[methotrexate_urogenital_system_lines, ]
methotrexate_test_urogenital_system <- methotrexate_test[which(methotrexate_test$Cell_line_tissue_type == 'urogenital_system'), ]
new_ic50 <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_urogenital_system), s = 'lambda.1se', interval = 'conf')
methotrexate_urogenital_system_least_auc <- auc(methotrexate_test_urogenital_system$least_sensitive, new_ic50)

methotrexate_least_auc <- c(methotrexate_aero_dig_tract_least_auc, methotrexate_bone_least_auc, 
                         methotrexate_breast_least_auc, methotrexate_digestive_system_least_auc, 
                         methotrexate_kidney_least_auc, methotrexate_large_intestine_least_auc, 
                         methotrexate_lung_least_auc, methotrexate_lung_NSCLC_least_auc, methotrexate_lung_SCLC_least_auc, 
                         methotrexate_nervous_system_least_auc, methotrexate_neuroblastoma_least_auc, 
                         methotrexate_pancreas_least_auc, methotrexate_skin_least_auc, 
                         methotrexate_soft_tissue_least_auc, methotrexate_thyroid_least_auc, 
                         methotrexate_urogenital_system_least_auc)

all_gdsc_auc_by_type <- data.frame(cisplatin_most_auc, cisplatin_least_auc, 
                                   etoposide_most_auc, etoposide_least_auc, 
                                   gemcitabine_most_auc, gemcitabine_least_auc, 
                                   methotrexate_most_auc, methotrexate_least_auc)

#all_gdsc_auc_by_type <- data.frame(t(all_gdsc_auc_by_type))
all_gdsc_auc_by_type <- rbind(all_gdsc_auc_by_type, overall_auc)
#all_gdsc_auc_by_type <- all_gdsc_auc_by_type[c(9, 1:8), ]
colnames(all_gdsc_auc_by_type) <- c('cisplatin_most', 'cisplatin_least', 
                                    'etoposide_most', 'etoposide_least', 'gemcitabine_most', 
                                    'gemcitabine_least', 'methotrexate_most', 'methotrexate_least')
rownames(all_gdsc_auc_by_type) <- c('aero_dig_tract', 'bone', 'breast', 'digestive_system', 'kidney', 
                                    'large_intestine', 'lung', 'lung_NSCLC', 'lung_SCLC', 
                                    'nervous_system', 'neuroblastoma', 'pancreas', 'skin', 
                                    'soft_tissue', 'thyroid', 'urogenital_system', 'overall')
for (i in 1:nrow(all_gdsc_auc_by_type)) {
  for (j in 1:ncol(all_gdsc_auc_by_type)) {
    if (all_gdsc_auc_by_type[i,j] == 'NaN') {
      all_gdsc_auc_by_type[i,j] <- 0.50
    }
  }
}
for (i in 1:nrow(all_gdsc_auc_by_type)) {
  for (j in 1:ncol(all_gdsc_auc_by_type)) {
    if (all_gdsc_auc_by_type[i,j] < 0.50) {
      all_gdsc_auc_by_type[i,j] <- 0.50
    }
  }
}
all_gdsc_auc_by_type <- round(all_gdsc_auc_by_type, digits = 2)
all_gdsc_auc_by_type <- data.frame(t(all_gdsc_auc_by_type))
all_gdsc_auc_by_type <- all_gdsc_auc_by_type[, c(17, 1:16)]
colors <- colorRampPalette(c("dodgerblue", "white", "red"))(100)

png(filename = 'Images/GDSC_AUC_heatmap.png', width = 960)
heatmap.2(as.matrix(all_gdsc_auc_by_type), trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, density.info = 'none', key.xlab = 'AUC', key.title = '', cexRow = 0.9, cexCol = 0.9, cellnote = all_gdsc_auc_by_type, notecol = 'black', colsep = 1, sepwidth = c(0.1,0.1), srtCol = 45, margins = c(5,8))
dev.off()

### load gene lists for each model ---
cisplatin_most_model_min <- read.csv('GLM_Models/cisplatin_most_model_min.csv')
cisplatin_least_model_min <- read.csv('GLM_Models/cisplatin_least_model_min.csv')
cisplatin_most_model_1se <- read.csv('GLM_Models/cisplatin_most_model_1se.csv')
cisplatin_least_model_1se <- read.csv('GLM_Models/cisplatin_least_model_1se.csv')

etoposide_most_model_min <- read.csv('GLM_Models/etoposide_most_model_min.csv')
etoposide_least_model_min <- read.csv('GLM_Models/etoposide_least_model_min.csv')
etoposide_most_model_1se <- read.csv('GLM_Models/etoposide_most_model_1se.csv')
etoposide_least_model_1se <- read.csv('GLM_Models/etoposide_least_model_1se.csv')

gemcitabine_most_model_min <- read.csv('GLM_Models/gemcitabine_most_model_min.csv')
gemcitabine_least_model_min <- read.csv('GLM_Models/gemcitabine_least_model_min.csv')
gemcitabine_most_model_1se <- read.csv('GLM_Models/gemcitabine_most_model_1se.csv')
gemcitabine_least_model_1se <- read.csv('GLM_Models/gemcitabine_least_model_1se.csv')

methotrexate_most_model_min <- read.csv('GLM_Models/methotrexate_most_model_min.csv')
methotrexate_least_model_min <- read.csv('GLM_Models/methotrexate_least_model_min.csv')
methotrexate_most_model_1se <- read.csv('GLM_Models/methotrexate_most_model_1se.csv')
methotrexate_least_model_1se <- read.csv('GLM_Models/methotrexate_least_model_1se.csv')

## create table containing number of genes in each model and overall AUC ----
# the number of genes used in min models
min_model_genes <- c(nrow(cisplatin_most_model_min), nrow(cisplatin_least_model_min), nrow(etoposide_most_model_min), nrow(etoposide_least_model_min), 
                     nrow(gemcitabine_most_model_min), nrow(gemcitabine_least_model_min), nrow(methotrexate_most_model_min), nrow(methotrexate_least_model_min))

# the number of genes used in 1se models
one_se_model_genes <- c(nrow(cisplatin_most_model_1se), nrow(cisplatin_least_model_1se), nrow(etoposide_most_model_1se), nrow(etoposide_least_model_1se), 
                        nrow(gemcitabine_most_model_1se), nrow(gemcitabine_least_model_1se), nrow(methotrexate_most_model_1se), nrow(methotrexate_least_model_1se))

# overall AUC values for min models
auc_min_models <- c(cisplatin_most_test_gdsc_auc_min, cisplatin_least_test_gdsc_auc_min, etoposide_most_test_gdsc_auc_min, etoposide_least_test_gdsc_auc_min, 
                    gemcitabine_most_test_gdsc_auc_min, gemcitabine_least_test_gdsc_auc_min, methotrexate_most_test_gdsc_auc_min, methotrexate_least_test_gdsc_auc_min)

# overall AUC values for 1se models
auc_1se_models <- c(cisplatin_most_test_gdsc_auc_1se, cisplatin_least_test_gdsc_auc_1se, etoposide_most_test_gdsc_auc_1se, etoposide_least_test_gdsc_auc_1se, 
                    gemcitabine_most_test_gdsc_auc_1se, gemcitabine_least_test_gdsc_auc_1se, methotrexate_most_test_gdsc_auc_1se, methotrexate_least_test_gdsc_auc_1se)
# put together
genes_and_auc_df <- data.frame(min_model_genes, auc_min_models, one_se_model_genes, auc_1se_models)
rownames(genes_and_auc_df) <- c('Cisplatin_most_sensitive', 'Cisplatin_least_sensitive', 'Etoposide_most_sensitive', 'Etoposide_least_sensitive', 
                                'Gemcitabine_most_sensitive', 'Gemcitabine_least_sensitive', 'Methotrexate_most_sensitive', 'Methotrexate_least_sensitive')
colnames(genes_and_auc_df) <- c('Number_genes_min_model', 'AUC_min_model', 'Number_genes_1se_model', 'AUC_1se_model')
genes_and_auc_df$Number_genes_lost <- genes_and_auc_df$Number_genes_min_model - genes_and_auc_df$Number_genes_1se_model
# creating table
genes_auc_formattable <- formattable(genes_and_auc_df, list(
  AUC_min_model = color_tile('white', 'red'),
  AUC_1se_model = color_tile('white', 'red'),
  area(col= c(Number_genes_min_model, Number_genes_1se_model, Number_genes_lost)) ~ normalize_bar('springgreen')
))
# and exporting it
webshot::install_phantomjs()
# from: https://github.com/renkun-ken/formattable/issues/26
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

export_formattable(genes_auc_formattable, file = 'Images/Genes_AUC_table.png')

### get auc for individual cancer types -------

# add TCGA classes to all testing data
# CISPLATIN
cisplatin_idx       <- gdsc_cell_line_info$COSMIC_ID %in% cisplatin_test$COSMIC_ID
gdsc_info_cisplatin <- gdsc_cell_line_info[cisplatin_idx, ]
cisplatin_test      <- merge(cisplatin_test, gdsc_info_cisplatin, by = 'COSMIC_ID')

# ETOPOSIDE
etoposide_idx       <- gdsc_cell_line_info$COSMIC_ID %in% etoposide_test$COSMIC_ID
gdsc_info_etoposide <- gdsc_cell_line_info[etoposide_idx, ]
etoposide_test      <- merge(etoposide_test, gdsc_info_etoposide, by = 'COSMIC_ID')

# GEMCITABINE
gemcitabine_idx       <- gdsc_cell_line_info$COSMIC_ID %in% gemcitabine_test$COSMIC_ID
gdsc_info_gemcitabine <- gdsc_cell_line_info[gemcitabine_idx, ]
gemcitabine_test      <- merge(gemcitabine_test, gdsc_info_gemcitabine, by = 'COSMIC_ID')

# METHOTREXATE
methotrexate_idx        <- gdsc_cell_line_info$COSMIC_ID %in% methotrexate_test$COSMIC_ID
gdsc_info_methotrexate  <- gdsc_cell_line_info[methotrexate_idx, ]
methotrexate_test       <- merge(methotrexate_test, gdsc_info_methotrexate, by = 'COSMIC_ID')

# AUC by TCGA class
## CISPLATIN
table(cisplatin_test$TCGA_class)

# BLCA 5
cisplatin_blca_idx <- cisplatin_test$TCGA_class == 'BLCA' #5

cisplatin_blca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_blca_idx]
actual_cisplatin_blca_most         <- cisplatin_test$most_sensitive[cisplatin_blca_idx]
cisplatin_blca_most_min_auc        <- auc(actual_cisplatin_blca_most, cisplatin_blca_most_sensitive_min)

cisplatin_blca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_blca_idx]
actual_cisplatin_blca_least <- cisplatin_test$least_sensitive[cisplatin_blca_idx]
cisplatin_blca_least_min_auc <- auc(actual_cisplatin_blca_least, cisplatin_blca_least_sensitive_min)

cisplatin_blca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_blca_idx]
actual_cisplatin_blca_most <- cisplatin_test$most_sensitive[cisplatin_blca_idx]
cisplatin_blca_most_1se_auc <- auc(actual_cisplatin_blca_most, cisplatin_blca_most_sensitive_1se)

cisplatin_blca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_blca_idx]
actual_cisplatin_blca_least <- cisplatin_test$least_sensitive[cisplatin_blca_idx]
cisplatin_blca_least_1se_auc <- auc(actual_cisplatin_blca_least, cisplatin_blca_least_sensitive_1se)

# BRCA
cisplatin_brca_idx <- cisplatin_test$TCGA_class == 'BRCA' #28

cisplatin_brca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_brca_idx]
actual_cisplatin_brca_most         <- cisplatin_test$most_sensitive[cisplatin_brca_idx]
cisplatin_brca_most_min_auc        <- auc(actual_cisplatin_brca_most, cisplatin_brca_most_sensitive_min)

cisplatin_brca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_brca_idx]
actual_cisplatin_brca_least <- cisplatin_test$least_sensitive[cisplatin_brca_idx]
cisplatin_brca_least_min_auc <- auc(actual_cisplatin_brca_least, cisplatin_brca_least_sensitive_min)

cisplatin_brca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_brca_idx]
actual_cisplatin_brca_most <- cisplatin_test$most_sensitive[cisplatin_brca_idx]
cisplatin_brca_most_1se_auc <- auc(actual_cisplatin_brca_most, cisplatin_brca_most_sensitive_1se)

cisplatin_brca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_brca_idx]
actual_cisplatin_brca_least <- cisplatin_test$least_sensitive[cisplatin_brca_idx]
cisplatin_brca_least_1se_auc <- auc(actual_cisplatin_brca_least, cisplatin_brca_least_sensitive_1se)

# CESC
cisplatin_cesc_idx <- cisplatin_test$TCGA_class == 'CESC' #4

cisplatin_cesc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_cesc_idx]
actual_cisplatin_cesc_most         <- cisplatin_test$most_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_most_min_auc        <- auc(actual_cisplatin_cesc_most, cisplatin_cesc_most_sensitive_min)

cisplatin_cesc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_cesc_idx]
actual_cisplatin_cesc_least <- cisplatin_test$least_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_least_min_auc <- auc(actual_cisplatin_cesc_least, cisplatin_cesc_least_sensitive_min)

cisplatin_cesc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_cesc_idx]
actual_cisplatin_cesc_most <- cisplatin_test$most_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_most_1se_auc <- auc(actual_cisplatin_cesc_most, cisplatin_cesc_most_sensitive_1se)

cisplatin_cesc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_cesc_idx]
actual_cisplatin_cesc_least <- cisplatin_test$least_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_least_1se_auc <- auc(actual_cisplatin_cesc_least, cisplatin_cesc_least_sensitive_1se)

# COREAD
cisplatin_coread_idx <- cisplatin_test$TCGA_class == 'COREAD' #18

cisplatin_coread_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_coread_idx]
actual_cisplatin_coread_most         <- cisplatin_test$most_sensitive[cisplatin_coread_idx]
cisplatin_coread_most_min_auc        <- auc(actual_cisplatin_coread_most, cisplatin_coread_most_sensitive_min)

cisplatin_coread_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_coread_idx]
actual_cisplatin_coread_least <- cisplatin_test$least_sensitive[cisplatin_coread_idx]
cisplatin_coread_least_min_auc <- auc(actual_cisplatin_coread_least, cisplatin_coread_least_sensitive_min)

cisplatin_coread_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_coread_idx]
actual_cisplatin_coread_most <- cisplatin_test$most_sensitive[cisplatin_coread_idx]
cisplatin_coread_most_1se_auc <- auc(actual_cisplatin_coread_most, cisplatin_coread_most_sensitive_1se)

cisplatin_coread_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_coread_idx]
actual_cisplatin_coread_least <- cisplatin_test$least_sensitive[cisplatin_coread_idx]
cisplatin_coread_least_1se_auc <- auc(actual_cisplatin_coread_least, cisplatin_coread_least_sensitive_1se)

# ESCA
cisplatin_esca_idx <- cisplatin_test$TCGA_class == 'ESCA' #13

cisplatin_esca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_esca_idx]
actual_cisplatin_esca_most         <- cisplatin_test$most_sensitive[cisplatin_esca_idx]
cisplatin_esca_most_min_auc        <- auc(actual_cisplatin_esca_most, cisplatin_esca_most_sensitive_min)

cisplatin_esca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_esca_idx]
actual_cisplatin_esca_least <- cisplatin_test$least_sensitive[cisplatin_esca_idx]
cisplatin_esca_least_min_auc <- auc(actual_cisplatin_esca_least, cisplatin_esca_least_sensitive_min)

cisplatin_esca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_esca_idx]
actual_cisplatin_esca_most <- cisplatin_test$most_sensitive[cisplatin_esca_idx]
cisplatin_esca_most_1se_auc <- auc(actual_cisplatin_esca_most, cisplatin_esca_most_sensitive_1se)

cisplatin_esca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_esca_idx]
actual_cisplatin_esca_least <- cisplatin_test$least_sensitive[cisplatin_esca_idx]
cisplatin_esca_least_1se_auc <- auc(actual_cisplatin_esca_least, cisplatin_esca_least_sensitive_1se)

# HNSC
cisplatin_hnsc_idx <- cisplatin_test$TCGA_class == 'HNSC' #19

cisplatin_hnsc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_most         <- cisplatin_test$most_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_most_min_auc        <- auc(actual_cisplatin_hnsc_most, cisplatin_hnsc_most_sensitive_min)

cisplatin_hnsc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_least <- cisplatin_test$least_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_least_min_auc <- auc(actual_cisplatin_hnsc_least, cisplatin_hnsc_least_sensitive_min)

cisplatin_hnsc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_most <- cisplatin_test$most_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_most_1se_auc <- auc(actual_cisplatin_hnsc_most, cisplatin_hnsc_most_sensitive_1se)

cisplatin_hnsc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_least <- cisplatin_test$least_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_least_1se_auc <- auc(actual_cisplatin_hnsc_least, cisplatin_hnsc_least_sensitive_1se)

# KIRC
cisplatin_kirc_idx <- cisplatin_test$TCGA_class == 'KIRC' #15

cisplatin_kirc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_kirc_idx]
actual_cisplatin_kirc_most         <- cisplatin_test$most_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_most_min_auc        <- auc(actual_cisplatin_kirc_most, cisplatin_kirc_most_sensitive_min)

cisplatin_kirc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_kirc_idx]
actual_cisplatin_kirc_least <- cisplatin_test$least_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_least_min_auc <- auc(actual_cisplatin_kirc_least, cisplatin_kirc_least_sensitive_min)

cisplatin_kirc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_kirc_idx]
actual_cisplatin_kirc_most <- cisplatin_test$most_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_most_1se_auc <- auc(actual_cisplatin_kirc_most, cisplatin_kirc_most_sensitive_1se)

cisplatin_kirc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_kirc_idx]
actual_cisplatin_kirc_least <- cisplatin_test$least_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_least_1se_auc <- auc(actual_cisplatin_kirc_least, cisplatin_kirc_least_sensitive_1se)

# LIHC
cisplatin_lihc_idx <- cisplatin_test$TCGA_class == 'LIHC' #7

cisplatin_lihc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_lihc_idx]
actual_cisplatin_lihc_most         <- cisplatin_test$most_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_most_min_auc        <- auc(actual_cisplatin_lihc_most, cisplatin_lihc_most_sensitive_min)

cisplatin_lihc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_lihc_idx]
actual_cisplatin_lihc_least <- cisplatin_test$least_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_least_min_auc <- auc(actual_cisplatin_lihc_least, cisplatin_lihc_least_sensitive_min)

cisplatin_lihc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_lihc_idx]
actual_cisplatin_lihc_most <- cisplatin_test$most_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_most_1se_auc <- auc(actual_cisplatin_lihc_most, cisplatin_lihc_most_sensitive_1se)

cisplatin_lihc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_lihc_idx]
actual_cisplatin_lihc_least <- cisplatin_test$least_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_least_1se_auc <- auc(actual_cisplatin_lihc_least, cisplatin_lihc_least_sensitive_1se)

# LUAD
cisplatin_luad_idx <- cisplatin_test$TCGA_class == 'LUAD' #26

cisplatin_luad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_luad_idx]
actual_cisplatin_luad_most         <- cisplatin_test$most_sensitive[cisplatin_luad_idx]
cisplatin_luad_most_min_auc        <- auc(actual_cisplatin_luad_most, cisplatin_luad_most_sensitive_min)

cisplatin_luad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_luad_idx]
actual_cisplatin_luad_least <- cisplatin_test$least_sensitive[cisplatin_luad_idx]
cisplatin_luad_least_min_auc <- auc(actual_cisplatin_luad_least, cisplatin_luad_least_sensitive_min)

cisplatin_luad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_luad_idx]
actual_cisplatin_luad_most <- cisplatin_test$most_sensitive[cisplatin_luad_idx]
cisplatin_luad_most_1se_auc <- auc(actual_cisplatin_luad_most, cisplatin_luad_most_sensitive_1se)

cisplatin_luad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_luad_idx]
actual_cisplatin_luad_least <- cisplatin_test$least_sensitive[cisplatin_luad_idx]
cisplatin_luad_least_1se_auc <- auc(actual_cisplatin_luad_least, cisplatin_luad_least_sensitive_1se)

# LUSC
cisplatin_lusc_idx <- cisplatin_test$TCGA_class == 'LUSC' #5

cisplatin_lusc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_lusc_idx]
actual_cisplatin_lusc_most         <- cisplatin_test$most_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_most_min_auc        <- auc(actual_cisplatin_lusc_most, cisplatin_lusc_most_sensitive_min)

cisplatin_lusc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_lusc_idx]
actual_cisplatin_lusc_least <- cisplatin_test$least_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_least_min_auc <- auc(actual_cisplatin_lusc_least, cisplatin_lusc_least_sensitive_min)

cisplatin_lusc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_lusc_idx]
actual_cisplatin_lusc_most <- cisplatin_test$most_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_most_1se_auc <- auc(actual_cisplatin_lusc_most, cisplatin_lusc_most_sensitive_1se)

cisplatin_lusc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_lusc_idx]
actual_cisplatin_lusc_least <- cisplatin_test$least_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_least_1se_auc <- auc(actual_cisplatin_lusc_least, cisplatin_lusc_least_sensitive_1se)

# MB
cisplatin_mb_idx <- cisplatin_test$TCGA_class == 'MB' #3

cisplatin_mb_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_mb_idx]
actual_cisplatin_mb_most         <- cisplatin_test$most_sensitive[cisplatin_mb_idx]
cisplatin_mb_most_min_auc        <- auc(actual_cisplatin_mb_most, cisplatin_mb_most_sensitive_min)

cisplatin_mb_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_mb_idx]
actual_cisplatin_mb_least <- cisplatin_test$least_sensitive[cisplatin_mb_idx]
cisplatin_mb_least_min_auc <- auc(actual_cisplatin_mb_least, cisplatin_mb_least_sensitive_min)

cisplatin_mb_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_mb_idx]
actual_cisplatin_mb_most <- cisplatin_test$most_sensitive[cisplatin_mb_idx]
cisplatin_mb_most_1se_auc <- auc(actual_cisplatin_mb_most, cisplatin_mb_most_sensitive_1se)

cisplatin_mb_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_mb_idx]
actual_cisplatin_mb_least <- cisplatin_test$least_sensitive[cisplatin_mb_idx]
cisplatin_mb_least_1se_auc <- auc(actual_cisplatin_mb_least, cisplatin_mb_least_sensitive_1se)

# MESO
cisplatin_meso_idx <- cisplatin_test$TCGA_class == 'MESO' #5

cisplatin_meso_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_meso_idx]
actual_cisplatin_meso_most         <- cisplatin_test$most_sensitive[cisplatin_meso_idx]
cisplatin_meso_most_min_auc        <- auc(actual_cisplatin_meso_most, cisplatin_meso_most_sensitive_min)

cisplatin_meso_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_meso_idx]
actual_cisplatin_meso_least <- cisplatin_test$least_sensitive[cisplatin_meso_idx]
cisplatin_meso_least_min_auc <- auc(actual_cisplatin_meso_least, cisplatin_meso_least_sensitive_min)

cisplatin_meso_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_meso_idx]
actual_cisplatin_meso_most <- cisplatin_test$most_sensitive[cisplatin_meso_idx]
cisplatin_meso_most_1se_auc <- auc(actual_cisplatin_meso_most, cisplatin_meso_most_sensitive_1se)

cisplatin_meso_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_meso_idx]
actual_cisplatin_meso_least <- cisplatin_test$least_sensitive[cisplatin_meso_idx]
cisplatin_meso_least_1se_auc <- auc(actual_cisplatin_meso_least, cisplatin_meso_least_sensitive_1se)

# NB
cisplatin_nb_idx <- cisplatin_test$TCGA_class == 'NB' #13

cisplatin_nb_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_nb_idx]
actual_cisplatin_nb_most         <- cisplatin_test$most_sensitive[cisplatin_nb_idx]
cisplatin_nb_most_min_auc        <- auc(actual_cisplatin_nb_most, cisplatin_nb_most_sensitive_min)

cisplatin_nb_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_nb_idx]
actual_cisplatin_nb_least <- cisplatin_test$least_sensitive[cisplatin_nb_idx]
cisplatin_nb_least_min_auc <- auc(actual_cisplatin_nb_least, cisplatin_nb_least_sensitive_min)

cisplatin_nb_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_nb_idx]
actual_cisplatin_nb_most <- cisplatin_test$most_sensitive[cisplatin_nb_idx]
cisplatin_nb_most_1se_auc <- auc(actual_cisplatin_nb_most, cisplatin_nb_most_sensitive_1se)

cisplatin_nb_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_nb_idx]
actual_cisplatin_nb_least <- cisplatin_test$least_sensitive[cisplatin_nb_idx]
cisplatin_nb_least_1se_auc <- auc(actual_cisplatin_nb_least, cisplatin_nb_least_sensitive_1se)

# OV
cisplatin_ov_idx <- cisplatin_test$TCGA_class == 'OV' #17

cisplatin_ov_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_ov_idx]
actual_cisplatin_ov_most         <- cisplatin_test$most_sensitive[cisplatin_ov_idx]
cisplatin_ov_most_min_auc        <- auc(actual_cisplatin_ov_most, cisplatin_ov_most_sensitive_min)

cisplatin_ov_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_ov_idx]
actual_cisplatin_ov_least <- cisplatin_test$least_sensitive[cisplatin_ov_idx]
cisplatin_ov_least_min_auc <- auc(actual_cisplatin_ov_least, cisplatin_ov_least_sensitive_min)

cisplatin_ov_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_ov_idx]
actual_cisplatin_ov_most <- cisplatin_test$most_sensitive[cisplatin_ov_idx]
cisplatin_ov_most_1se_auc <- auc(actual_cisplatin_ov_most, cisplatin_ov_most_sensitive_1se)

cisplatin_ov_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_ov_idx]
actual_cisplatin_ov_least <- cisplatin_test$least_sensitive[cisplatin_ov_idx]
cisplatin_ov_least_1se_auc <- auc(actual_cisplatin_ov_least, cisplatin_ov_least_sensitive_1se)

# PAAD
cisplatin_paad_idx <- cisplatin_test$TCGA_class == 'PAAD' #13

cisplatin_paad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_paad_idx]
actual_cisplatin_paad_most         <- cisplatin_test$most_sensitive[cisplatin_paad_idx]
cisplatin_paad_most_min_auc        <- auc(actual_cisplatin_paad_most, cisplatin_paad_most_sensitive_min)

cisplatin_paad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_paad_idx]
actual_cisplatin_paad_least <- cisplatin_test$least_sensitive[cisplatin_paad_idx]
cisplatin_paad_least_min_auc <- auc(actual_cisplatin_paad_least, cisplatin_paad_least_sensitive_min)

cisplatin_paad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_paad_idx]
actual_cisplatin_paad_most <- cisplatin_test$most_sensitive[cisplatin_paad_idx]
cisplatin_paad_most_1se_auc <- auc(actual_cisplatin_paad_most, cisplatin_paad_most_sensitive_1se)

cisplatin_paad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_paad_idx]
actual_cisplatin_paad_least <- cisplatin_test$least_sensitive[cisplatin_paad_idx]
cisplatin_paad_least_1se_auc <- auc(actual_cisplatin_paad_least, cisplatin_paad_least_sensitive_1se)

# SCLC
cisplatin_sclc_idx <- cisplatin_test$TCGA_class == 'SCLC' #31

cisplatin_sclc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_sclc_idx]
actual_cisplatin_sclc_most         <- cisplatin_test$most_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_most_min_auc        <- auc(actual_cisplatin_sclc_most, cisplatin_sclc_most_sensitive_min)

cisplatin_sclc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_sclc_idx]
actual_cisplatin_sclc_least <- cisplatin_test$least_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_least_min_auc <- auc(actual_cisplatin_sclc_least, cisplatin_sclc_least_sensitive_min)

cisplatin_sclc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_sclc_idx]
actual_cisplatin_sclc_most <- cisplatin_test$most_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_most_1se_auc <- auc(actual_cisplatin_sclc_most, cisplatin_sclc_most_sensitive_1se)

cisplatin_sclc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_sclc_idx]
actual_cisplatin_sclc_least <- cisplatin_test$least_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_least_1se_auc <- auc(actual_cisplatin_sclc_least, cisplatin_sclc_least_sensitive_1se)

# SKCM
cisplatin_skcm_idx <- cisplatin_test$TCGA_class == 'SKCM' #26

cisplatin_skcm_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_skcm_idx]
actual_cisplatin_skcm_most         <- cisplatin_test$most_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_most_min_auc        <- auc(actual_cisplatin_skcm_most, cisplatin_skcm_most_sensitive_min)

cisplatin_skcm_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_skcm_idx]
actual_cisplatin_skcm_least <- cisplatin_test$least_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_least_min_auc <- auc(actual_cisplatin_skcm_least, cisplatin_skcm_least_sensitive_min)

cisplatin_skcm_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_skcm_idx]
actual_cisplatin_skcm_most <- cisplatin_test$most_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_most_1se_auc <- auc(actual_cisplatin_skcm_most, cisplatin_skcm_most_sensitive_1se)

cisplatin_skcm_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_skcm_idx]
actual_cisplatin_skcm_least <- cisplatin_test$least_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_least_1se_auc <- auc(actual_cisplatin_skcm_least, cisplatin_skcm_least_sensitive_1se)

# STAD
cisplatin_stad_idx <- cisplatin_test$TCGA_class == 'STAD' #8

cisplatin_stad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_stad_idx]
actual_cisplatin_stad_most         <- cisplatin_test$most_sensitive[cisplatin_stad_idx]
cisplatin_stad_most_min_auc        <- auc(actual_cisplatin_stad_most, cisplatin_stad_most_sensitive_min)

cisplatin_stad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_stad_idx]
actual_cisplatin_stad_least <- cisplatin_test$least_sensitive[cisplatin_stad_idx]
cisplatin_stad_least_min_auc <- auc(actual_cisplatin_stad_least, cisplatin_stad_least_sensitive_min)

cisplatin_stad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_stad_idx]
actual_cisplatin_stad_most <- cisplatin_test$most_sensitive[cisplatin_stad_idx]
cisplatin_stad_most_1se_auc <- auc(actual_cisplatin_stad_most, cisplatin_stad_most_sensitive_1se)

cisplatin_stad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_stad_idx]
actual_cisplatin_stad_least <- cisplatin_test$least_sensitive[cisplatin_stad_idx]
cisplatin_stad_least_1se_auc <- auc(actual_cisplatin_stad_least, cisplatin_stad_least_sensitive_1se)

# UCEC
cisplatin_ucec_idx <- cisplatin_test$TCGA_class == 'UCEC' #4

cisplatin_ucec_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_ucec_idx]
actual_cisplatin_ucec_most         <- cisplatin_test$most_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_most_min_auc        <- auc(actual_cisplatin_ucec_most, cisplatin_ucec_most_sensitive_min)

cisplatin_ucec_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_ucec_idx]
actual_cisplatin_ucec_least <- cisplatin_test$least_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_least_min_auc <- auc(actual_cisplatin_ucec_least, cisplatin_ucec_least_sensitive_min)

cisplatin_ucec_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_ucec_idx]
actual_cisplatin_ucec_most <- cisplatin_test$most_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_most_1se_auc <- auc(actual_cisplatin_ucec_most, cisplatin_ucec_most_sensitive_1se)

cisplatin_ucec_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_ucec_idx]
actual_cisplatin_ucec_least <- cisplatin_test$least_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_least_1se_auc <- auc(actual_cisplatin_ucec_least, cisplatin_ucec_least_sensitive_1se)

# UNCLASSIFIED
cisplatin_unclassified_idx <- cisplatin_test$TCGA_class == 'UNCLASSIFIED' #52

cisplatin_unclassified_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_most         <- cisplatin_test$most_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_most_min_auc        <- auc(actual_cisplatin_unclassified_most, cisplatin_unclassified_most_sensitive_min)

cisplatin_unclassified_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_least <- cisplatin_test$least_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_least_min_auc <- auc(actual_cisplatin_unclassified_least, cisplatin_unclassified_least_sensitive_min)

cisplatin_unclassified_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_most <- cisplatin_test$most_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_most_1se_auc <- auc(actual_cisplatin_unclassified_most, cisplatin_unclassified_most_sensitive_1se)

cisplatin_unclassified_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_least <- cisplatin_test$least_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_least_1se_auc <- auc(actual_cisplatin_unclassified_least, cisplatin_unclassified_least_sensitive_1se)

# ETOPOSIDE
etoposide_blca_idx <- etoposide_test$TCGA_class == 'BLCA' #5

etoposide_blca_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_blca_idx]
actual_etoposide_blca_most         <- etoposide_test$most_sensitive[etoposide_blca_idx]
etoposide_blca_most_min_auc        <- auc(actual_etoposide_blca_most, etoposide_blca_most_sensitive_min)

etoposide_blca_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_blca_idx]
actual_etoposide_blca_least <- etoposide_test$least_sensitive[etoposide_blca_idx]
etoposide_blca_least_min_auc <- auc(actual_etoposide_blca_least, etoposide_blca_least_sensitive_min)

etoposide_blca_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_blca_idx]
actual_etoposide_blca_most <- etoposide_test$most_sensitive[etoposide_blca_idx]
etoposide_blca_most_1se_auc <- auc(actual_etoposide_blca_most, etoposide_blca_most_sensitive_1se)

etoposide_blca_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_blca_idx]
actual_etoposide_blca_least <- etoposide_test$least_sensitive[etoposide_blca_idx]
etoposide_blca_least_1se_auc <- auc(actual_etoposide_blca_least, etoposide_blca_least_sensitive_1se)

# BRCA
etoposide_brca_idx <- etoposide_test$TCGA_class == 'BRCA' #28

etoposide_brca_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_brca_idx]
actual_etoposide_brca_most         <- etoposide_test$most_sensitive[etoposide_brca_idx]
etoposide_brca_most_min_auc        <- auc(actual_etoposide_brca_most, etoposide_brca_most_sensitive_min)

etoposide_brca_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_brca_idx]
actual_etoposide_brca_least <- etoposide_test$least_sensitive[etoposide_brca_idx]
etoposide_brca_least_min_auc <- auc(actual_etoposide_brca_least, etoposide_brca_least_sensitive_min)

etoposide_brca_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_brca_idx]
actual_etoposide_brca_most <- etoposide_test$most_sensitive[etoposide_brca_idx]
etoposide_brca_most_1se_auc <- auc(actual_etoposide_brca_most, etoposide_brca_most_sensitive_1se)

etoposide_brca_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_brca_idx]
actual_etoposide_brca_least <- etoposide_test$least_sensitive[etoposide_brca_idx]
etoposide_brca_least_1se_auc <- auc(actual_etoposide_brca_least, etoposide_brca_least_sensitive_1se)

# CESC
etoposide_cesc_idx <- etoposide_test$TCGA_class == 'CESC' #4

etoposide_cesc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_cesc_idx]
actual_etoposide_cesc_most         <- etoposide_test$most_sensitive[etoposide_cesc_idx]
etoposide_cesc_most_min_auc        <- auc(actual_etoposide_cesc_most, etoposide_cesc_most_sensitive_min)

etoposide_cesc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_cesc_idx]
actual_etoposide_cesc_least <- etoposide_test$least_sensitive[etoposide_cesc_idx]
etoposide_cesc_least_min_auc <- auc(actual_etoposide_cesc_least, etoposide_cesc_least_sensitive_min)

etoposide_cesc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_cesc_idx]
actual_etoposide_cesc_most <- etoposide_test$most_sensitive[etoposide_cesc_idx]
etoposide_cesc_most_1se_auc <- auc(actual_etoposide_cesc_most, etoposide_cesc_most_sensitive_1se)

etoposide_cesc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_cesc_idx]
actual_etoposide_cesc_least <- etoposide_test$least_sensitive[etoposide_cesc_idx]
etoposide_cesc_least_1se_auc <- auc(actual_etoposide_cesc_least, etoposide_cesc_least_sensitive_1se)

# COREAD
etoposide_coread_idx <- etoposide_test$TCGA_class == 'COREAD' #18

etoposide_coread_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_coread_idx]
actual_etoposide_coread_most         <- etoposide_test$most_sensitive[etoposide_coread_idx]
etoposide_coread_most_min_auc        <- auc(actual_etoposide_coread_most, etoposide_coread_most_sensitive_min)

etoposide_coread_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_coread_idx]
actual_etoposide_coread_least <- etoposide_test$least_sensitive[etoposide_coread_idx]
etoposide_coread_least_min_auc <- auc(actual_etoposide_coread_least, etoposide_coread_least_sensitive_min)

etoposide_coread_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_coread_idx]
actual_etoposide_coread_most <- etoposide_test$most_sensitive[etoposide_coread_idx]
etoposide_coread_most_1se_auc <- auc(actual_etoposide_coread_most, etoposide_coread_most_sensitive_1se)

etoposide_coread_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_coread_idx]
actual_etoposide_coread_least <- etoposide_test$least_sensitive[etoposide_coread_idx]
etoposide_coread_least_1se_auc <- auc(actual_etoposide_coread_least, etoposide_coread_least_sensitive_1se)

# ESCA
etoposide_esca_idx <- etoposide_test$TCGA_class == 'ESCA' #13

etoposide_esca_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_esca_idx]
actual_etoposide_esca_most         <- etoposide_test$most_sensitive[etoposide_esca_idx]
etoposide_esca_most_min_auc        <- auc(actual_etoposide_esca_most, etoposide_esca_most_sensitive_min)

etoposide_esca_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_esca_idx]
actual_etoposide_esca_least <- etoposide_test$least_sensitive[etoposide_esca_idx]
etoposide_esca_least_min_auc <- auc(actual_etoposide_esca_least, etoposide_esca_least_sensitive_min)

etoposide_esca_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_esca_idx]
actual_etoposide_esca_most <- etoposide_test$most_sensitive[etoposide_esca_idx]
etoposide_esca_most_1se_auc <- auc(actual_etoposide_esca_most, etoposide_esca_most_sensitive_1se)

etoposide_esca_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_esca_idx]
actual_etoposide_esca_least <- etoposide_test$least_sensitive[etoposide_esca_idx]
etoposide_esca_least_1se_auc <- auc(actual_etoposide_esca_least, etoposide_esca_least_sensitive_1se)

# HNSC
etoposide_hnsc_idx <- etoposide_test$TCGA_class == 'HNSC' #19

etoposide_hnsc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_hnsc_idx]
actual_etoposide_hnsc_most         <- etoposide_test$most_sensitive[etoposide_hnsc_idx]
etoposide_hnsc_most_min_auc        <- auc(actual_etoposide_hnsc_most, etoposide_hnsc_most_sensitive_min)

etoposide_hnsc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_hnsc_idx]
actual_etoposide_hnsc_least <- etoposide_test$least_sensitive[etoposide_hnsc_idx]
etoposide_hnsc_least_min_auc <- auc(actual_etoposide_hnsc_least, etoposide_hnsc_least_sensitive_min)

etoposide_hnsc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_hnsc_idx]
actual_etoposide_hnsc_most <- etoposide_test$most_sensitive[etoposide_hnsc_idx]
etoposide_hnsc_most_1se_auc <- auc(actual_etoposide_hnsc_most, etoposide_hnsc_most_sensitive_1se)

etoposide_hnsc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_hnsc_idx]
actual_etoposide_hnsc_least <- etoposide_test$least_sensitive[etoposide_hnsc_idx]
etoposide_hnsc_least_1se_auc <- auc(actual_etoposide_hnsc_least, etoposide_hnsc_least_sensitive_1se)

# KIRC
etoposide_kirc_idx <- etoposide_test$TCGA_class == 'KIRC' #15

etoposide_kirc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_kirc_idx]
actual_etoposide_kirc_most         <- etoposide_test$most_sensitive[etoposide_kirc_idx]
etoposide_kirc_most_min_auc        <- auc(actual_etoposide_kirc_most, etoposide_kirc_most_sensitive_min)

etoposide_kirc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_kirc_idx]
actual_etoposide_kirc_least <- etoposide_test$least_sensitive[etoposide_kirc_idx]
etoposide_kirc_least_min_auc <- auc(actual_etoposide_kirc_least, etoposide_kirc_least_sensitive_min)

etoposide_kirc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_kirc_idx]
actual_etoposide_kirc_most <- etoposide_test$most_sensitive[etoposide_kirc_idx]
etoposide_kirc_most_1se_auc <- auc(actual_etoposide_kirc_most, etoposide_kirc_most_sensitive_1se)

etoposide_kirc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_kirc_idx]
actual_etoposide_kirc_least <- etoposide_test$least_sensitive[etoposide_kirc_idx]
etoposide_kirc_least_1se_auc <- auc(actual_etoposide_kirc_least, etoposide_kirc_least_sensitive_1se)

# LIHC
etoposide_lihc_idx <- etoposide_test$TCGA_class == 'LIHC' #7

etoposide_lihc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_lihc_idx]
actual_etoposide_lihc_most         <- etoposide_test$most_sensitive[etoposide_lihc_idx]
etoposide_lihc_most_min_auc        <- auc(actual_etoposide_lihc_most, etoposide_lihc_most_sensitive_min)

etoposide_lihc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_lihc_idx]
actual_etoposide_lihc_least <- etoposide_test$least_sensitive[etoposide_lihc_idx]
etoposide_lihc_least_min_auc <- auc(actual_etoposide_lihc_least, etoposide_lihc_least_sensitive_min)

etoposide_lihc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_lihc_idx]
actual_etoposide_lihc_most <- etoposide_test$most_sensitive[etoposide_lihc_idx]
etoposide_lihc_most_1se_auc <- auc(actual_etoposide_lihc_most, etoposide_lihc_most_sensitive_1se)

etoposide_lihc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_lihc_idx]
actual_etoposide_lihc_least <- etoposide_test$least_sensitive[etoposide_lihc_idx]
etoposide_lihc_least_1se_auc <- auc(actual_etoposide_lihc_least, etoposide_lihc_least_sensitive_1se)

# LUAD
etoposide_luad_idx <- etoposide_test$TCGA_class == 'LUAD' #26

etoposide_luad_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_luad_idx]
actual_etoposide_luad_most         <- etoposide_test$most_sensitive[etoposide_luad_idx]
etoposide_luad_most_min_auc        <- auc(actual_etoposide_luad_most, etoposide_luad_most_sensitive_min)

etoposide_luad_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_luad_idx]
actual_etoposide_luad_least <- etoposide_test$least_sensitive[etoposide_luad_idx]
etoposide_luad_least_min_auc <- auc(actual_etoposide_luad_least, etoposide_luad_least_sensitive_min)

etoposide_luad_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_luad_idx]
actual_etoposide_luad_most <- etoposide_test$most_sensitive[etoposide_luad_idx]
etoposide_luad_most_1se_auc <- auc(actual_etoposide_luad_most, etoposide_luad_most_sensitive_1se)

etoposide_luad_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_luad_idx]
actual_etoposide_luad_least <- etoposide_test$least_sensitive[etoposide_luad_idx]
etoposide_luad_least_1se_auc <- auc(actual_etoposide_luad_least, etoposide_luad_least_sensitive_1se)

# LUSC
etoposide_lusc_idx <- etoposide_test$TCGA_class == 'LUSC' #5

etoposide_lusc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_lusc_idx]
actual_etoposide_lusc_most         <- etoposide_test$most_sensitive[etoposide_lusc_idx]
etoposide_lusc_most_min_auc        <- auc(actual_etoposide_lusc_most, etoposide_lusc_most_sensitive_min)

etoposide_lusc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_lusc_idx]
actual_etoposide_lusc_least <- etoposide_test$least_sensitive[etoposide_lusc_idx]
etoposide_lusc_least_min_auc <- auc(actual_etoposide_lusc_least, etoposide_lusc_least_sensitive_min)

etoposide_lusc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_lusc_idx]
actual_etoposide_lusc_most <- etoposide_test$most_sensitive[etoposide_lusc_idx]
etoposide_lusc_most_1se_auc <- auc(actual_etoposide_lusc_most, etoposide_lusc_most_sensitive_1se)

etoposide_lusc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_lusc_idx]
actual_etoposide_lusc_least <- etoposide_test$least_sensitive[etoposide_lusc_idx]
etoposide_lusc_least_1se_auc <- auc(actual_etoposide_lusc_least, etoposide_lusc_least_sensitive_1se)

# MB
etoposide_mb_idx <- etoposide_test$TCGA_class == 'MB' #3

etoposide_mb_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_mb_idx]
actual_etoposide_mb_most         <- etoposide_test$most_sensitive[etoposide_mb_idx]
etoposide_mb_most_min_auc        <- auc(actual_etoposide_mb_most, etoposide_mb_most_sensitive_min)

etoposide_mb_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_mb_idx]
actual_etoposide_mb_least <- etoposide_test$least_sensitive[etoposide_mb_idx]
etoposide_mb_least_min_auc <- auc(actual_etoposide_mb_least, etoposide_mb_least_sensitive_min)

etoposide_mb_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_mb_idx]
actual_etoposide_mb_most <- etoposide_test$most_sensitive[etoposide_mb_idx]
etoposide_mb_most_1se_auc <- auc(actual_etoposide_mb_most, etoposide_mb_most_sensitive_1se)

etoposide_mb_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_mb_idx]
actual_etoposide_mb_least <- etoposide_test$least_sensitive[etoposide_mb_idx]
etoposide_mb_least_1se_auc <- auc(actual_etoposide_mb_least, etoposide_mb_least_sensitive_1se)

# MESO
etoposide_meso_idx <- etoposide_test$TCGA_class == 'MESO' #5

etoposide_meso_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_meso_idx]
actual_etoposide_meso_most         <- etoposide_test$most_sensitive[etoposide_meso_idx]
etoposide_meso_most_min_auc        <- auc(actual_etoposide_meso_most, etoposide_meso_most_sensitive_min)

etoposide_meso_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_meso_idx]
actual_etoposide_meso_least <- etoposide_test$least_sensitive[etoposide_meso_idx]
etoposide_meso_least_min_auc <- auc(actual_etoposide_meso_least, etoposide_meso_least_sensitive_min)

etoposide_meso_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_meso_idx]
actual_etoposide_meso_most <- etoposide_test$most_sensitive[etoposide_meso_idx]
etoposide_meso_most_1se_auc <- auc(actual_etoposide_meso_most, etoposide_meso_most_sensitive_1se)

etoposide_meso_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_meso_idx]
actual_etoposide_meso_least <- etoposide_test$least_sensitive[etoposide_meso_idx]
etoposide_meso_least_1se_auc <- auc(actual_etoposide_meso_least, etoposide_meso_least_sensitive_1se)

# NB
etoposide_nb_idx <- etoposide_test$TCGA_class == 'NB' #13

etoposide_nb_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_nb_idx]
actual_etoposide_nb_most         <- etoposide_test$most_sensitive[etoposide_nb_idx]
etoposide_nb_most_min_auc        <- auc(actual_etoposide_nb_most, etoposide_nb_most_sensitive_min)

etoposide_nb_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_nb_idx]
actual_etoposide_nb_least <- etoposide_test$least_sensitive[etoposide_nb_idx]
etoposide_nb_least_min_auc <- auc(actual_etoposide_nb_least, etoposide_nb_least_sensitive_min)

etoposide_nb_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_nb_idx]
actual_etoposide_nb_most <- etoposide_test$most_sensitive[etoposide_nb_idx]
etoposide_nb_most_1se_auc <- auc(actual_etoposide_nb_most, etoposide_nb_most_sensitive_1se)

etoposide_nb_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_nb_idx]
actual_etoposide_nb_least <- etoposide_test$least_sensitive[etoposide_nb_idx]
etoposide_nb_least_1se_auc <- auc(actual_etoposide_nb_least, etoposide_nb_least_sensitive_1se)

# OV
etoposide_ov_idx <- etoposide_test$TCGA_class == 'OV' #17

etoposide_ov_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_ov_idx]
actual_etoposide_ov_most         <- etoposide_test$most_sensitive[etoposide_ov_idx]
etoposide_ov_most_min_auc        <- auc(actual_etoposide_ov_most, etoposide_ov_most_sensitive_min)

etoposide_ov_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_ov_idx]
actual_etoposide_ov_least <- etoposide_test$least_sensitive[etoposide_ov_idx]
etoposide_ov_least_min_auc <- auc(actual_etoposide_ov_least, etoposide_ov_least_sensitive_min)

etoposide_ov_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_ov_idx]
actual_etoposide_ov_most <- etoposide_test$most_sensitive[etoposide_ov_idx]
etoposide_ov_most_1se_auc <- auc(actual_etoposide_ov_most, etoposide_ov_most_sensitive_1se)

etoposide_ov_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_ov_idx]
actual_etoposide_ov_least <- etoposide_test$least_sensitive[etoposide_ov_idx]
etoposide_ov_least_1se_auc <- auc(actual_etoposide_ov_least, etoposide_ov_least_sensitive_1se)

# PAAD
etoposide_paad_idx <- etoposide_test$TCGA_class == 'PAAD' #13

etoposide_paad_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_paad_idx]
actual_etoposide_paad_most         <- etoposide_test$most_sensitive[etoposide_paad_idx]
etoposide_paad_most_min_auc        <- auc(actual_etoposide_paad_most, etoposide_paad_most_sensitive_min)

etoposide_paad_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_paad_idx]
actual_etoposide_paad_least <- etoposide_test$least_sensitive[etoposide_paad_idx]
etoposide_paad_least_min_auc <- auc(actual_etoposide_paad_least, etoposide_paad_least_sensitive_min)

etoposide_paad_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_paad_idx]
actual_etoposide_paad_most <- etoposide_test$most_sensitive[etoposide_paad_idx]
etoposide_paad_most_1se_auc <- auc(actual_etoposide_paad_most, etoposide_paad_most_sensitive_1se)

etoposide_paad_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_paad_idx]
actual_etoposide_paad_least <- etoposide_test$least_sensitive[etoposide_paad_idx]
etoposide_paad_least_1se_auc <- auc(actual_etoposide_paad_least, etoposide_paad_least_sensitive_1se)

# SCLC
etoposide_sclc_idx <- etoposide_test$TCGA_class == 'SCLC' #31

etoposide_sclc_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_sclc_idx]
actual_etoposide_sclc_most         <- etoposide_test$most_sensitive[etoposide_sclc_idx]
etoposide_sclc_most_min_auc        <- auc(actual_etoposide_sclc_most, etoposide_sclc_most_sensitive_min)

etoposide_sclc_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_sclc_idx]
actual_etoposide_sclc_least <- etoposide_test$least_sensitive[etoposide_sclc_idx]
etoposide_sclc_least_min_auc <- auc(actual_etoposide_sclc_least, etoposide_sclc_least_sensitive_min)

etoposide_sclc_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_sclc_idx]
actual_etoposide_sclc_most <- etoposide_test$most_sensitive[etoposide_sclc_idx]
etoposide_sclc_most_1se_auc <- auc(actual_etoposide_sclc_most, etoposide_sclc_most_sensitive_1se)

etoposide_sclc_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_sclc_idx]
actual_etoposide_sclc_least <- etoposide_test$least_sensitive[etoposide_sclc_idx]
etoposide_sclc_least_1se_auc <- auc(actual_etoposide_sclc_least, etoposide_sclc_least_sensitive_1se)

# SKCM
etoposide_skcm_idx <- etoposide_test$TCGA_class == 'SKCM' #26

etoposide_skcm_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_skcm_idx]
actual_etoposide_skcm_most         <- etoposide_test$most_sensitive[etoposide_skcm_idx]
etoposide_skcm_most_min_auc        <- auc(actual_etoposide_skcm_most, etoposide_skcm_most_sensitive_min)

etoposide_skcm_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_skcm_idx]
actual_etoposide_skcm_least <- etoposide_test$least_sensitive[etoposide_skcm_idx]
etoposide_skcm_least_min_auc <- auc(actual_etoposide_skcm_least, etoposide_skcm_least_sensitive_min)

etoposide_skcm_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_skcm_idx]
actual_etoposide_skcm_most <- etoposide_test$most_sensitive[etoposide_skcm_idx]
etoposide_skcm_most_1se_auc <- auc(actual_etoposide_skcm_most, etoposide_skcm_most_sensitive_1se)

etoposide_skcm_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_skcm_idx]
actual_etoposide_skcm_least <- etoposide_test$least_sensitive[etoposide_skcm_idx]
etoposide_skcm_least_1se_auc <- auc(actual_etoposide_skcm_least, etoposide_skcm_least_sensitive_1se)

# STAD
etoposide_stad_idx <- etoposide_test$TCGA_class == 'STAD' #8

etoposide_stad_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_stad_idx]
actual_etoposide_stad_most         <- etoposide_test$most_sensitive[etoposide_stad_idx]
etoposide_stad_most_min_auc        <- auc(actual_etoposide_stad_most, etoposide_stad_most_sensitive_min)

etoposide_stad_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_stad_idx]
actual_etoposide_stad_least <- etoposide_test$least_sensitive[etoposide_stad_idx]
etoposide_stad_least_min_auc <- auc(actual_etoposide_stad_least, etoposide_stad_least_sensitive_min)

etoposide_stad_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_stad_idx]
actual_etoposide_stad_most <- etoposide_test$most_sensitive[etoposide_stad_idx]
etoposide_stad_most_1se_auc <- auc(actual_etoposide_stad_most, etoposide_stad_most_sensitive_1se)

etoposide_stad_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_stad_idx]
actual_etoposide_stad_least <- etoposide_test$least_sensitive[etoposide_stad_idx]
etoposide_stad_least_1se_auc <- auc(actual_etoposide_stad_least, etoposide_stad_least_sensitive_1se)

# UCEC
etoposide_ucec_idx <- etoposide_test$TCGA_class == 'UCEC' #4

etoposide_ucec_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_ucec_idx]
actual_etoposide_ucec_most         <- etoposide_test$most_sensitive[etoposide_ucec_idx]
etoposide_ucec_most_min_auc        <- auc(actual_etoposide_ucec_most, etoposide_ucec_most_sensitive_min)

etoposide_ucec_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_ucec_idx]
actual_etoposide_ucec_least <- etoposide_test$least_sensitive[etoposide_ucec_idx]
etoposide_ucec_least_min_auc <- auc(actual_etoposide_ucec_least, etoposide_ucec_least_sensitive_min)

etoposide_ucec_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_ucec_idx]
actual_etoposide_ucec_most <- etoposide_test$most_sensitive[etoposide_ucec_idx]
etoposide_ucec_most_1se_auc <- auc(actual_etoposide_ucec_most, etoposide_ucec_most_sensitive_1se)

etoposide_ucec_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_ucec_idx]
actual_etoposide_ucec_least <- etoposide_test$least_sensitive[etoposide_ucec_idx]
etoposide_ucec_least_1se_auc <- auc(actual_etoposide_ucec_least, etoposide_ucec_least_sensitive_1se)

# UNCLASSIFIED
etoposide_unclassified_idx <- etoposide_test$TCGA_class == 'UNCLASSIFIED' #52

etoposide_unclassified_most_sensitive_min  <- new_etoposide_most_sensitive_min[etoposide_unclassified_idx]
actual_etoposide_unclassified_most         <- etoposide_test$most_sensitive[etoposide_unclassified_idx]
etoposide_unclassified_most_min_auc        <- auc(actual_etoposide_unclassified_most, etoposide_unclassified_most_sensitive_min)

etoposide_unclassified_least_sensitive_min <- new_etoposide_least_sensitive_min[etoposide_unclassified_idx]
actual_etoposide_unclassified_least <- etoposide_test$least_sensitive[etoposide_unclassified_idx]
etoposide_unclassified_least_min_auc <- auc(actual_etoposide_unclassified_least, etoposide_unclassified_least_sensitive_min)

etoposide_unclassified_most_sensitive_1se <- new_etoposide_most_sensitive_1se[etoposide_unclassified_idx]
actual_etoposide_unclassified_most <- etoposide_test$most_sensitive[etoposide_unclassified_idx]
etoposide_unclassified_most_1se_auc <- auc(actual_etoposide_unclassified_most, etoposide_unclassified_most_sensitive_1se)

etoposide_unclassified_least_sensitive_1se <- new_etoposide_least_sensitive_1se[etoposide_unclassified_idx]
actual_etoposide_unclassified_least <- etoposide_test$least_sensitive[etoposide_unclassified_idx]
etoposide_unclassified_least_1se_auc <- auc(actual_etoposide_unclassified_least, etoposide_unclassified_least_sensitive_1se)

# GEMICTABINE
gemcitabine_blca_idx <- gemcitabine_test$TCGA_class == 'BLCA' #5

gemcitabine_blca_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_blca_idx]
actual_gemcitabine_blca_most         <- gemcitabine_test$most_sensitive[gemcitabine_blca_idx]
gemcitabine_blca_most_min_auc        <- auc(actual_gemcitabine_blca_most, gemcitabine_blca_most_sensitive_min)

gemcitabine_blca_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_blca_idx]
actual_gemcitabine_blca_least <- gemcitabine_test$least_sensitive[gemcitabine_blca_idx]
gemcitabine_blca_least_min_auc <- auc(actual_gemcitabine_blca_least, gemcitabine_blca_least_sensitive_min)

gemcitabine_blca_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_blca_idx]
actual_gemcitabine_blca_most <- gemcitabine_test$most_sensitive[gemcitabine_blca_idx]
gemcitabine_blca_most_1se_auc <- auc(actual_gemcitabine_blca_most, gemcitabine_blca_most_sensitive_1se)

gemcitabine_blca_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_blca_idx]
actual_gemcitabine_blca_least <- gemcitabine_test$least_sensitive[gemcitabine_blca_idx]
gemcitabine_blca_least_1se_auc <- auc(actual_gemcitabine_blca_least, gemcitabine_blca_least_sensitive_1se)

# BRCA
gemcitabine_brca_idx <- gemcitabine_test$TCGA_class == 'BRCA' #28

gemcitabine_brca_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_brca_idx]
actual_gemcitabine_brca_most         <- gemcitabine_test$most_sensitive[gemcitabine_brca_idx]
gemcitabine_brca_most_min_auc        <- auc(actual_gemcitabine_brca_most, gemcitabine_brca_most_sensitive_min)

gemcitabine_brca_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_brca_idx]
actual_gemcitabine_brca_least <- gemcitabine_test$least_sensitive[gemcitabine_brca_idx]
gemcitabine_brca_least_min_auc <- auc(actual_gemcitabine_brca_least, gemcitabine_brca_least_sensitive_min)

gemcitabine_brca_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_brca_idx]
actual_gemcitabine_brca_most <- gemcitabine_test$most_sensitive[gemcitabine_brca_idx]
gemcitabine_brca_most_1se_auc <- auc(actual_gemcitabine_brca_most, gemcitabine_brca_most_sensitive_1se)

gemcitabine_brca_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_brca_idx]
actual_gemcitabine_brca_least <- gemcitabine_test$least_sensitive[gemcitabine_brca_idx]
gemcitabine_brca_least_1se_auc <- auc(actual_gemcitabine_brca_least, gemcitabine_brca_least_sensitive_1se)

# CESC
gemcitabine_cesc_idx <- gemcitabine_test$TCGA_class == 'CESC' #4

gemcitabine_cesc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_cesc_idx]
actual_gemcitabine_cesc_most         <- gemcitabine_test$most_sensitive[gemcitabine_cesc_idx]
gemcitabine_cesc_most_min_auc        <- auc(actual_gemcitabine_cesc_most, gemcitabine_cesc_most_sensitive_min)

gemcitabine_cesc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_cesc_idx]
actual_gemcitabine_cesc_least <- gemcitabine_test$least_sensitive[gemcitabine_cesc_idx]
gemcitabine_cesc_least_min_auc <- auc(actual_gemcitabine_cesc_least, gemcitabine_cesc_least_sensitive_min)

gemcitabine_cesc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_cesc_idx]
actual_gemcitabine_cesc_most <- gemcitabine_test$most_sensitive[gemcitabine_cesc_idx]
gemcitabine_cesc_most_1se_auc <- auc(actual_gemcitabine_cesc_most, gemcitabine_cesc_most_sensitive_1se)

gemcitabine_cesc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_cesc_idx]
actual_gemcitabine_cesc_least <- gemcitabine_test$least_sensitive[gemcitabine_cesc_idx]
gemcitabine_cesc_least_1se_auc <- auc(actual_gemcitabine_cesc_least, gemcitabine_cesc_least_sensitive_1se)

# COREAD
gemcitabine_coread_idx <- gemcitabine_test$TCGA_class == 'COREAD' #18

gemcitabine_coread_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_coread_idx]
actual_gemcitabine_coread_most         <- gemcitabine_test$most_sensitive[gemcitabine_coread_idx]
gemcitabine_coread_most_min_auc        <- auc(actual_gemcitabine_coread_most, gemcitabine_coread_most_sensitive_min)

gemcitabine_coread_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_coread_idx]
actual_gemcitabine_coread_least <- gemcitabine_test$least_sensitive[gemcitabine_coread_idx]
gemcitabine_coread_least_min_auc <- auc(actual_gemcitabine_coread_least, gemcitabine_coread_least_sensitive_min)

gemcitabine_coread_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_coread_idx]
actual_gemcitabine_coread_most <- gemcitabine_test$most_sensitive[gemcitabine_coread_idx]
gemcitabine_coread_most_1se_auc <- auc(actual_gemcitabine_coread_most, gemcitabine_coread_most_sensitive_1se)

gemcitabine_coread_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_coread_idx]
actual_gemcitabine_coread_least <- gemcitabine_test$least_sensitive[gemcitabine_coread_idx]
gemcitabine_coread_least_1se_auc <- auc(actual_gemcitabine_coread_least, gemcitabine_coread_least_sensitive_1se)

# ESCA
gemcitabine_esca_idx <- gemcitabine_test$TCGA_class == 'ESCA' #13

gemcitabine_esca_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_esca_idx]
actual_gemcitabine_esca_most         <- gemcitabine_test$most_sensitive[gemcitabine_esca_idx]
gemcitabine_esca_most_min_auc        <- auc(actual_gemcitabine_esca_most, gemcitabine_esca_most_sensitive_min)

gemcitabine_esca_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_esca_idx]
actual_gemcitabine_esca_least <- gemcitabine_test$least_sensitive[gemcitabine_esca_idx]
gemcitabine_esca_least_min_auc <- auc(actual_gemcitabine_esca_least, gemcitabine_esca_least_sensitive_min)

gemcitabine_esca_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_esca_idx]
actual_gemcitabine_esca_most <- gemcitabine_test$most_sensitive[gemcitabine_esca_idx]
gemcitabine_esca_most_1se_auc <- auc(actual_gemcitabine_esca_most, gemcitabine_esca_most_sensitive_1se)

gemcitabine_esca_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_esca_idx]
actual_gemcitabine_esca_least <- gemcitabine_test$least_sensitive[gemcitabine_esca_idx]
gemcitabine_esca_least_1se_auc <- auc(actual_gemcitabine_esca_least, gemcitabine_esca_least_sensitive_1se)

# HNSC
gemcitabine_hnsc_idx <- gemcitabine_test$TCGA_class == 'HNSC' #19

gemcitabine_hnsc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_hnsc_idx]
actual_gemcitabine_hnsc_most         <- gemcitabine_test$most_sensitive[gemcitabine_hnsc_idx]
gemcitabine_hnsc_most_min_auc        <- auc(actual_gemcitabine_hnsc_most, gemcitabine_hnsc_most_sensitive_min)

gemcitabine_hnsc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_hnsc_idx]
actual_gemcitabine_hnsc_least <- gemcitabine_test$least_sensitive[gemcitabine_hnsc_idx]
gemcitabine_hnsc_least_min_auc <- auc(actual_gemcitabine_hnsc_least, gemcitabine_hnsc_least_sensitive_min)

gemcitabine_hnsc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_hnsc_idx]
actual_gemcitabine_hnsc_most <- gemcitabine_test$most_sensitive[gemcitabine_hnsc_idx]
gemcitabine_hnsc_most_1se_auc <- auc(actual_gemcitabine_hnsc_most, gemcitabine_hnsc_most_sensitive_1se)

gemcitabine_hnsc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_hnsc_idx]
actual_gemcitabine_hnsc_least <- gemcitabine_test$least_sensitive[gemcitabine_hnsc_idx]
gemcitabine_hnsc_least_1se_auc <- auc(actual_gemcitabine_hnsc_least, gemcitabine_hnsc_least_sensitive_1se)

# KIRC
gemcitabine_kirc_idx <- gemcitabine_test$TCGA_class == 'KIRC' #15

gemcitabine_kirc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_kirc_idx]
actual_gemcitabine_kirc_most         <- gemcitabine_test$most_sensitive[gemcitabine_kirc_idx]
gemcitabine_kirc_most_min_auc        <- auc(actual_gemcitabine_kirc_most, gemcitabine_kirc_most_sensitive_min)

gemcitabine_kirc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_kirc_idx]
actual_gemcitabine_kirc_least <- gemcitabine_test$least_sensitive[gemcitabine_kirc_idx]
gemcitabine_kirc_least_min_auc <- auc(actual_gemcitabine_kirc_least, gemcitabine_kirc_least_sensitive_min)

gemcitabine_kirc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_kirc_idx]
actual_gemcitabine_kirc_most <- gemcitabine_test$most_sensitive[gemcitabine_kirc_idx]
gemcitabine_kirc_most_1se_auc <- auc(actual_gemcitabine_kirc_most, gemcitabine_kirc_most_sensitive_1se)

gemcitabine_kirc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_kirc_idx]
actual_gemcitabine_kirc_least <- gemcitabine_test$least_sensitive[gemcitabine_kirc_idx]
gemcitabine_kirc_least_1se_auc <- auc(actual_gemcitabine_kirc_least, gemcitabine_kirc_least_sensitive_1se)

# LIHC
gemcitabine_lihc_idx <- gemcitabine_test$TCGA_class == 'LIHC' #7

gemcitabine_lihc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_lihc_idx]
actual_gemcitabine_lihc_most         <- gemcitabine_test$most_sensitive[gemcitabine_lihc_idx]
gemcitabine_lihc_most_min_auc        <- auc(actual_gemcitabine_lihc_most, gemcitabine_lihc_most_sensitive_min)

gemcitabine_lihc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_lihc_idx]
actual_gemcitabine_lihc_least <- gemcitabine_test$least_sensitive[gemcitabine_lihc_idx]
gemcitabine_lihc_least_min_auc <- auc(actual_gemcitabine_lihc_least, gemcitabine_lihc_least_sensitive_min)

gemcitabine_lihc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_lihc_idx]
actual_gemcitabine_lihc_most <- gemcitabine_test$most_sensitive[gemcitabine_lihc_idx]
gemcitabine_lihc_most_1se_auc <- auc(actual_gemcitabine_lihc_most, gemcitabine_lihc_most_sensitive_1se)

gemcitabine_lihc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_lihc_idx]
actual_gemcitabine_lihc_least <- gemcitabine_test$least_sensitive[gemcitabine_lihc_idx]
gemcitabine_lihc_least_1se_auc <- auc(actual_gemcitabine_lihc_least, gemcitabine_lihc_least_sensitive_1se)

# LUAD
gemcitabine_luad_idx <- gemcitabine_test$TCGA_class == 'LUAD' #26

gemcitabine_luad_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_luad_idx]
actual_gemcitabine_luad_most         <- gemcitabine_test$most_sensitive[gemcitabine_luad_idx]
gemcitabine_luad_most_min_auc        <- auc(actual_gemcitabine_luad_most, gemcitabine_luad_most_sensitive_min)

gemcitabine_luad_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_luad_idx]
actual_gemcitabine_luad_least <- gemcitabine_test$least_sensitive[gemcitabine_luad_idx]
gemcitabine_luad_least_min_auc <- auc(actual_gemcitabine_luad_least, gemcitabine_luad_least_sensitive_min)

gemcitabine_luad_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_luad_idx]
actual_gemcitabine_luad_most <- gemcitabine_test$most_sensitive[gemcitabine_luad_idx]
gemcitabine_luad_most_1se_auc <- auc(actual_gemcitabine_luad_most, gemcitabine_luad_most_sensitive_1se)

gemcitabine_luad_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_luad_idx]
actual_gemcitabine_luad_least <- gemcitabine_test$least_sensitive[gemcitabine_luad_idx]
gemcitabine_luad_least_1se_auc <- auc(actual_gemcitabine_luad_least, gemcitabine_luad_least_sensitive_1se)

# LUSC
gemcitabine_lusc_idx <- gemcitabine_test$TCGA_class == 'LUSC' #5

gemcitabine_lusc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_lusc_idx]
actual_gemcitabine_lusc_most         <- gemcitabine_test$most_sensitive[gemcitabine_lusc_idx]
gemcitabine_lusc_most_min_auc        <- auc(actual_gemcitabine_lusc_most, gemcitabine_lusc_most_sensitive_min)

gemcitabine_lusc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_lusc_idx]
actual_gemcitabine_lusc_least <- gemcitabine_test$least_sensitive[gemcitabine_lusc_idx]
gemcitabine_lusc_least_min_auc <- auc(actual_gemcitabine_lusc_least, gemcitabine_lusc_least_sensitive_min)

gemcitabine_lusc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_lusc_idx]
actual_gemcitabine_lusc_most <- gemcitabine_test$most_sensitive[gemcitabine_lusc_idx]
gemcitabine_lusc_most_1se_auc <- auc(actual_gemcitabine_lusc_most, gemcitabine_lusc_most_sensitive_1se)

gemcitabine_lusc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_lusc_idx]
actual_gemcitabine_lusc_least <- gemcitabine_test$least_sensitive[gemcitabine_lusc_idx]
gemcitabine_lusc_least_1se_auc <- auc(actual_gemcitabine_lusc_least, gemcitabine_lusc_least_sensitive_1se)

# MB
gemcitabine_mb_idx <- gemcitabine_test$TCGA_class == 'MB' #3

gemcitabine_mb_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_mb_idx]
actual_gemcitabine_mb_most         <- gemcitabine_test$most_sensitive[gemcitabine_mb_idx]
gemcitabine_mb_most_min_auc        <- auc(actual_gemcitabine_mb_most, gemcitabine_mb_most_sensitive_min)

gemcitabine_mb_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_mb_idx]
actual_gemcitabine_mb_least <- gemcitabine_test$least_sensitive[gemcitabine_mb_idx]
gemcitabine_mb_least_min_auc <- auc(actual_gemcitabine_mb_least, gemcitabine_mb_least_sensitive_min)

gemcitabine_mb_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_mb_idx]
actual_gemcitabine_mb_most <- gemcitabine_test$most_sensitive[gemcitabine_mb_idx]
gemcitabine_mb_most_1se_auc <- auc(actual_gemcitabine_mb_most, gemcitabine_mb_most_sensitive_1se)

gemcitabine_mb_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_mb_idx]
actual_gemcitabine_mb_least <- gemcitabine_test$least_sensitive[gemcitabine_mb_idx]
gemcitabine_mb_least_1se_auc <- auc(actual_gemcitabine_mb_least, gemcitabine_mb_least_sensitive_1se)

# MESO
gemcitabine_meso_idx <- gemcitabine_test$TCGA_class == 'MESO' #5

gemcitabine_meso_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_meso_idx]
actual_gemcitabine_meso_most         <- gemcitabine_test$most_sensitive[gemcitabine_meso_idx]
gemcitabine_meso_most_min_auc        <- auc(actual_gemcitabine_meso_most, gemcitabine_meso_most_sensitive_min)

gemcitabine_meso_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_meso_idx]
actual_gemcitabine_meso_least <- gemcitabine_test$least_sensitive[gemcitabine_meso_idx]
gemcitabine_meso_least_min_auc <- auc(actual_gemcitabine_meso_least, gemcitabine_meso_least_sensitive_min)

gemcitabine_meso_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_meso_idx]
actual_gemcitabine_meso_most <- gemcitabine_test$most_sensitive[gemcitabine_meso_idx]
gemcitabine_meso_most_1se_auc <- auc(actual_gemcitabine_meso_most, gemcitabine_meso_most_sensitive_1se)

gemcitabine_meso_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_meso_idx]
actual_gemcitabine_meso_least <- gemcitabine_test$least_sensitive[gemcitabine_meso_idx]
gemcitabine_meso_least_1se_auc <- auc(actual_gemcitabine_meso_least, gemcitabine_meso_least_sensitive_1se)

# NB
gemcitabine_nb_idx <- gemcitabine_test$TCGA_class == 'NB' #13

gemcitabine_nb_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_nb_idx]
actual_gemcitabine_nb_most         <- gemcitabine_test$most_sensitive[gemcitabine_nb_idx]
gemcitabine_nb_most_min_auc        <- auc(actual_gemcitabine_nb_most, gemcitabine_nb_most_sensitive_min)

gemcitabine_nb_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_nb_idx]
actual_gemcitabine_nb_least <- gemcitabine_test$least_sensitive[gemcitabine_nb_idx]
gemcitabine_nb_least_min_auc <- auc(actual_gemcitabine_nb_least, gemcitabine_nb_least_sensitive_min)

gemcitabine_nb_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_nb_idx]
actual_gemcitabine_nb_most <- gemcitabine_test$most_sensitive[gemcitabine_nb_idx]
gemcitabine_nb_most_1se_auc <- auc(actual_gemcitabine_nb_most, gemcitabine_nb_most_sensitive_1se)

gemcitabine_nb_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_nb_idx]
actual_gemcitabine_nb_least <- gemcitabine_test$least_sensitive[gemcitabine_nb_idx]
gemcitabine_nb_least_1se_auc <- auc(actual_gemcitabine_nb_least, gemcitabine_nb_least_sensitive_1se)

# OV
gemcitabine_ov_idx <- gemcitabine_test$TCGA_class == 'OV' #17

gemcitabine_ov_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_ov_idx]
actual_gemcitabine_ov_most         <- gemcitabine_test$most_sensitive[gemcitabine_ov_idx]
gemcitabine_ov_most_min_auc        <- auc(actual_gemcitabine_ov_most, gemcitabine_ov_most_sensitive_min)

gemcitabine_ov_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_ov_idx]
actual_gemcitabine_ov_least <- gemcitabine_test$least_sensitive[gemcitabine_ov_idx]
gemcitabine_ov_least_min_auc <- auc(actual_gemcitabine_ov_least, gemcitabine_ov_least_sensitive_min)

gemcitabine_ov_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_ov_idx]
actual_gemcitabine_ov_most <- gemcitabine_test$most_sensitive[gemcitabine_ov_idx]
gemcitabine_ov_most_1se_auc <- auc(actual_gemcitabine_ov_most, gemcitabine_ov_most_sensitive_1se)

gemcitabine_ov_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_ov_idx]
actual_gemcitabine_ov_least <- gemcitabine_test$least_sensitive[gemcitabine_ov_idx]
gemcitabine_ov_least_1se_auc <- auc(actual_gemcitabine_ov_least, gemcitabine_ov_least_sensitive_1se)

# PAAD
gemcitabine_paad_idx <- gemcitabine_test$TCGA_class == 'PAAD' #13

gemcitabine_paad_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_paad_idx]
actual_gemcitabine_paad_most         <- gemcitabine_test$most_sensitive[gemcitabine_paad_idx]
gemcitabine_paad_most_min_auc        <- auc(actual_gemcitabine_paad_most, gemcitabine_paad_most_sensitive_min)

gemcitabine_paad_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_paad_idx]
actual_gemcitabine_paad_least <- gemcitabine_test$least_sensitive[gemcitabine_paad_idx]
gemcitabine_paad_least_min_auc <- auc(actual_gemcitabine_paad_least, gemcitabine_paad_least_sensitive_min)

gemcitabine_paad_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_paad_idx]
actual_gemcitabine_paad_most <- gemcitabine_test$most_sensitive[gemcitabine_paad_idx]
gemcitabine_paad_most_1se_auc <- auc(actual_gemcitabine_paad_most, gemcitabine_paad_most_sensitive_1se)

gemcitabine_paad_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_paad_idx]
actual_gemcitabine_paad_least <- gemcitabine_test$least_sensitive[gemcitabine_paad_idx]
gemcitabine_paad_least_1se_auc <- auc(actual_gemcitabine_paad_least, gemcitabine_paad_least_sensitive_1se)

# SCLC
gemcitabine_sclc_idx <- gemcitabine_test$TCGA_class == 'SCLC' #31

gemcitabine_sclc_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_sclc_idx]
actual_gemcitabine_sclc_most         <- gemcitabine_test$most_sensitive[gemcitabine_sclc_idx]
gemcitabine_sclc_most_min_auc        <- auc(actual_gemcitabine_sclc_most, gemcitabine_sclc_most_sensitive_min)

gemcitabine_sclc_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_sclc_idx]
actual_gemcitabine_sclc_least <- gemcitabine_test$least_sensitive[gemcitabine_sclc_idx]
gemcitabine_sclc_least_min_auc <- auc(actual_gemcitabine_sclc_least, gemcitabine_sclc_least_sensitive_min)

gemcitabine_sclc_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_sclc_idx]
actual_gemcitabine_sclc_most <- gemcitabine_test$most_sensitive[gemcitabine_sclc_idx]
gemcitabine_sclc_most_1se_auc <- auc(actual_gemcitabine_sclc_most, gemcitabine_sclc_most_sensitive_1se)

gemcitabine_sclc_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_sclc_idx]
actual_gemcitabine_sclc_least <- gemcitabine_test$least_sensitive[gemcitabine_sclc_idx]
gemcitabine_sclc_least_1se_auc <- auc(actual_gemcitabine_sclc_least, gemcitabine_sclc_least_sensitive_1se)

# SKCM
gemcitabine_skcm_idx <- gemcitabine_test$TCGA_class == 'SKCM' #26

gemcitabine_skcm_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_skcm_idx]
actual_gemcitabine_skcm_most         <- gemcitabine_test$most_sensitive[gemcitabine_skcm_idx]
gemcitabine_skcm_most_min_auc        <- auc(actual_gemcitabine_skcm_most, gemcitabine_skcm_most_sensitive_min)

gemcitabine_skcm_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_skcm_idx]
actual_gemcitabine_skcm_least <- gemcitabine_test$least_sensitive[gemcitabine_skcm_idx]
gemcitabine_skcm_least_min_auc <- auc(actual_gemcitabine_skcm_least, gemcitabine_skcm_least_sensitive_min)

gemcitabine_skcm_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_skcm_idx]
actual_gemcitabine_skcm_most <- gemcitabine_test$most_sensitive[gemcitabine_skcm_idx]
gemcitabine_skcm_most_1se_auc <- auc(actual_gemcitabine_skcm_most, gemcitabine_skcm_most_sensitive_1se)

gemcitabine_skcm_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_skcm_idx]
actual_gemcitabine_skcm_least <- gemcitabine_test$least_sensitive[gemcitabine_skcm_idx]
gemcitabine_skcm_least_1se_auc <- auc(actual_gemcitabine_skcm_least, gemcitabine_skcm_least_sensitive_1se)

# STAD
gemcitabine_stad_idx <- gemcitabine_test$TCGA_class == 'STAD' #8

gemcitabine_stad_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_stad_idx]
actual_gemcitabine_stad_most         <- gemcitabine_test$most_sensitive[gemcitabine_stad_idx]
gemcitabine_stad_most_min_auc        <- auc(actual_gemcitabine_stad_most, gemcitabine_stad_most_sensitive_min)

gemcitabine_stad_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_stad_idx]
actual_gemcitabine_stad_least <- gemcitabine_test$least_sensitive[gemcitabine_stad_idx]
gemcitabine_stad_least_min_auc <- auc(actual_gemcitabine_stad_least, gemcitabine_stad_least_sensitive_min)

gemcitabine_stad_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_stad_idx]
actual_gemcitabine_stad_most <- gemcitabine_test$most_sensitive[gemcitabine_stad_idx]
gemcitabine_stad_most_1se_auc <- auc(actual_gemcitabine_stad_most, gemcitabine_stad_most_sensitive_1se)

gemcitabine_stad_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_stad_idx]
actual_gemcitabine_stad_least <- gemcitabine_test$least_sensitive[gemcitabine_stad_idx]
gemcitabine_stad_least_1se_auc <- auc(actual_gemcitabine_stad_least, gemcitabine_stad_least_sensitive_1se)

# UCEC
gemcitabine_ucec_idx <- gemcitabine_test$TCGA_class == 'UCEC' #4

gemcitabine_ucec_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_ucec_idx]
actual_gemcitabine_ucec_most         <- gemcitabine_test$most_sensitive[gemcitabine_ucec_idx]
gemcitabine_ucec_most_min_auc        <- auc(actual_gemcitabine_ucec_most, gemcitabine_ucec_most_sensitive_min)

gemcitabine_ucec_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_ucec_idx]
actual_gemcitabine_ucec_least <- gemcitabine_test$least_sensitive[gemcitabine_ucec_idx]
gemcitabine_ucec_least_min_auc <- auc(actual_gemcitabine_ucec_least, gemcitabine_ucec_least_sensitive_min)

gemcitabine_ucec_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_ucec_idx]
actual_gemcitabine_ucec_most <- gemcitabine_test$most_sensitive[gemcitabine_ucec_idx]
gemcitabine_ucec_most_1se_auc <- auc(actual_gemcitabine_ucec_most, gemcitabine_ucec_most_sensitive_1se)

gemcitabine_ucec_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_ucec_idx]
actual_gemcitabine_ucec_least <- gemcitabine_test$least_sensitive[gemcitabine_ucec_idx]
gemcitabine_ucec_least_1se_auc <- auc(actual_gemcitabine_ucec_least, gemcitabine_ucec_least_sensitive_1se)

# UNCLASSIFIED
gemcitabine_unclassified_idx <- gemcitabine_test$TCGA_class == 'UNCLASSIFIED' #52

gemcitabine_unclassified_most_sensitive_min  <- new_gemcitabine_most_sensitive_min[gemcitabine_unclassified_idx]
actual_gemcitabine_unclassified_most         <- gemcitabine_test$most_sensitive[gemcitabine_unclassified_idx]
gemcitabine_unclassified_most_min_auc        <- auc(actual_gemcitabine_unclassified_most, gemcitabine_unclassified_most_sensitive_min)

gemcitabine_unclassified_least_sensitive_min <- new_gemcitabine_least_sensitive_min[gemcitabine_unclassified_idx]
actual_gemcitabine_unclassified_least <- gemcitabine_test$least_sensitive[gemcitabine_unclassified_idx]
gemcitabine_unclassified_least_min_auc <- auc(actual_gemcitabine_unclassified_least, gemcitabine_unclassified_least_sensitive_min)

gemcitabine_unclassified_most_sensitive_1se <- new_gemcitabine_most_sensitive_1se[gemcitabine_unclassified_idx]
actual_gemcitabine_unclassified_most <- gemcitabine_test$most_sensitive[gemcitabine_unclassified_idx]
gemcitabine_unclassified_most_1se_auc <- auc(actual_gemcitabine_unclassified_most, gemcitabine_unclassified_most_sensitive_1se)

gemcitabine_unclassified_least_sensitive_1se <- new_gemcitabine_least_sensitive_1se[gemcitabine_unclassified_idx]
actual_gemcitabine_unclassified_least <- gemcitabine_test$least_sensitive[gemcitabine_unclassified_idx]
gemcitabine_unclassified_least_1se_auc <- auc(actual_gemcitabine_unclassified_least, gemcitabine_unclassified_least_sensitive_1se)

# METHOTREXATE
methotrexate_blca_idx <- methotrexate_test$TCGA_class == 'BLCA' #5

methotrexate_blca_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_blca_idx]
actual_methotrexate_blca_most         <- methotrexate_test$most_sensitive[methotrexate_blca_idx]
methotrexate_blca_most_min_auc        <- auc(actual_methotrexate_blca_most, methotrexate_blca_most_sensitive_min)

methotrexate_blca_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_blca_idx]
actual_methotrexate_blca_least <- methotrexate_test$least_sensitive[methotrexate_blca_idx]
methotrexate_blca_least_min_auc <- auc(actual_methotrexate_blca_least, methotrexate_blca_least_sensitive_min)

methotrexate_blca_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_blca_idx]
actual_methotrexate_blca_most <- methotrexate_test$most_sensitive[methotrexate_blca_idx]
methotrexate_blca_most_1se_auc <- auc(actual_methotrexate_blca_most, methotrexate_blca_most_sensitive_1se)

methotrexate_blca_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_blca_idx]
actual_methotrexate_blca_least <- methotrexate_test$least_sensitive[methotrexate_blca_idx]
methotrexate_blca_least_1se_auc <- auc(actual_methotrexate_blca_least, methotrexate_blca_least_sensitive_1se)

# BRCA
methotrexate_brca_idx <- methotrexate_test$TCGA_class == 'BRCA' #28

methotrexate_brca_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_brca_idx]
actual_methotrexate_brca_most         <- methotrexate_test$most_sensitive[methotrexate_brca_idx]
methotrexate_brca_most_min_auc        <- auc(actual_methotrexate_brca_most, methotrexate_brca_most_sensitive_min)

methotrexate_brca_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_brca_idx]
actual_methotrexate_brca_least <- methotrexate_test$least_sensitive[methotrexate_brca_idx]
methotrexate_brca_least_min_auc <- auc(actual_methotrexate_brca_least, methotrexate_brca_least_sensitive_min)

methotrexate_brca_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_brca_idx]
actual_methotrexate_brca_most <- methotrexate_test$most_sensitive[methotrexate_brca_idx]
methotrexate_brca_most_1se_auc <- auc(actual_methotrexate_brca_most, methotrexate_brca_most_sensitive_1se)

methotrexate_brca_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_brca_idx]
actual_methotrexate_brca_least <- methotrexate_test$least_sensitive[methotrexate_brca_idx]
methotrexate_brca_least_1se_auc <- auc(actual_methotrexate_brca_least, methotrexate_brca_least_sensitive_1se)

# CESC
methotrexate_cesc_idx <- methotrexate_test$TCGA_class == 'CESC' #4

methotrexate_cesc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_cesc_idx]
actual_methotrexate_cesc_most         <- methotrexate_test$most_sensitive[methotrexate_cesc_idx]
methotrexate_cesc_most_min_auc        <- auc(actual_methotrexate_cesc_most, methotrexate_cesc_most_sensitive_min)

methotrexate_cesc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_cesc_idx]
actual_methotrexate_cesc_least <- methotrexate_test$least_sensitive[methotrexate_cesc_idx]
methotrexate_cesc_least_min_auc <- auc(actual_methotrexate_cesc_least, methotrexate_cesc_least_sensitive_min)

methotrexate_cesc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_cesc_idx]
actual_methotrexate_cesc_most <- methotrexate_test$most_sensitive[methotrexate_cesc_idx]
methotrexate_cesc_most_1se_auc <- auc(actual_methotrexate_cesc_most, methotrexate_cesc_most_sensitive_1se)

methotrexate_cesc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_cesc_idx]
actual_methotrexate_cesc_least <- methotrexate_test$least_sensitive[methotrexate_cesc_idx]
methotrexate_cesc_least_1se_auc <- auc(actual_methotrexate_cesc_least, methotrexate_cesc_least_sensitive_1se)

# COREAD
methotrexate_coread_idx <- methotrexate_test$TCGA_class == 'COREAD' #18

methotrexate_coread_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_coread_idx]
actual_methotrexate_coread_most         <- methotrexate_test$most_sensitive[methotrexate_coread_idx]
methotrexate_coread_most_min_auc        <- auc(actual_methotrexate_coread_most, methotrexate_coread_most_sensitive_min)

methotrexate_coread_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_coread_idx]
actual_methotrexate_coread_least <- methotrexate_test$least_sensitive[methotrexate_coread_idx]
methotrexate_coread_least_min_auc <- auc(actual_methotrexate_coread_least, methotrexate_coread_least_sensitive_min)

methotrexate_coread_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_coread_idx]
actual_methotrexate_coread_most <- methotrexate_test$most_sensitive[methotrexate_coread_idx]
methotrexate_coread_most_1se_auc <- auc(actual_methotrexate_coread_most, methotrexate_coread_most_sensitive_1se)

methotrexate_coread_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_coread_idx]
actual_methotrexate_coread_least <- methotrexate_test$least_sensitive[methotrexate_coread_idx]
methotrexate_coread_least_1se_auc <- auc(actual_methotrexate_coread_least, methotrexate_coread_least_sensitive_1se)

# ESCA
methotrexate_esca_idx <- methotrexate_test$TCGA_class == 'ESCA' #13

methotrexate_esca_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_esca_idx]
actual_methotrexate_esca_most         <- methotrexate_test$most_sensitive[methotrexate_esca_idx]
methotrexate_esca_most_min_auc        <- auc(actual_methotrexate_esca_most, methotrexate_esca_most_sensitive_min)

methotrexate_esca_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_esca_idx]
actual_methotrexate_esca_least <- methotrexate_test$least_sensitive[methotrexate_esca_idx]
methotrexate_esca_least_min_auc <- auc(actual_methotrexate_esca_least, methotrexate_esca_least_sensitive_min)

methotrexate_esca_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_esca_idx]
actual_methotrexate_esca_most <- methotrexate_test$most_sensitive[methotrexate_esca_idx]
methotrexate_esca_most_1se_auc <- auc(actual_methotrexate_esca_most, methotrexate_esca_most_sensitive_1se)

methotrexate_esca_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_esca_idx]
actual_methotrexate_esca_least <- methotrexate_test$least_sensitive[methotrexate_esca_idx]
methotrexate_esca_least_1se_auc <- auc(actual_methotrexate_esca_least, methotrexate_esca_least_sensitive_1se)

# HNSC
methotrexate_hnsc_idx <- methotrexate_test$TCGA_class == 'HNSC' #19

methotrexate_hnsc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_hnsc_idx]
actual_methotrexate_hnsc_most         <- methotrexate_test$most_sensitive[methotrexate_hnsc_idx]
methotrexate_hnsc_most_min_auc        <- auc(actual_methotrexate_hnsc_most, methotrexate_hnsc_most_sensitive_min)

methotrexate_hnsc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_hnsc_idx]
actual_methotrexate_hnsc_least <- methotrexate_test$least_sensitive[methotrexate_hnsc_idx]
methotrexate_hnsc_least_min_auc <- auc(actual_methotrexate_hnsc_least, methotrexate_hnsc_least_sensitive_min)

methotrexate_hnsc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_hnsc_idx]
actual_methotrexate_hnsc_most <- methotrexate_test$most_sensitive[methotrexate_hnsc_idx]
methotrexate_hnsc_most_1se_auc <- auc(actual_methotrexate_hnsc_most, methotrexate_hnsc_most_sensitive_1se)

methotrexate_hnsc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_hnsc_idx]
actual_methotrexate_hnsc_least <- methotrexate_test$least_sensitive[methotrexate_hnsc_idx]
methotrexate_hnsc_least_1se_auc <- auc(actual_methotrexate_hnsc_least, methotrexate_hnsc_least_sensitive_1se)

# KIRC
methotrexate_kirc_idx <- methotrexate_test$TCGA_class == 'KIRC' #15

methotrexate_kirc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_kirc_idx]
actual_methotrexate_kirc_most         <- methotrexate_test$most_sensitive[methotrexate_kirc_idx]
methotrexate_kirc_most_min_auc        <- auc(actual_methotrexate_kirc_most, methotrexate_kirc_most_sensitive_min)

methotrexate_kirc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_kirc_idx]
actual_methotrexate_kirc_least <- methotrexate_test$least_sensitive[methotrexate_kirc_idx]
methotrexate_kirc_least_min_auc <- auc(actual_methotrexate_kirc_least, methotrexate_kirc_least_sensitive_min)

methotrexate_kirc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_kirc_idx]
actual_methotrexate_kirc_most <- methotrexate_test$most_sensitive[methotrexate_kirc_idx]
methotrexate_kirc_most_1se_auc <- auc(actual_methotrexate_kirc_most, methotrexate_kirc_most_sensitive_1se)

methotrexate_kirc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_kirc_idx]
actual_methotrexate_kirc_least <- methotrexate_test$least_sensitive[methotrexate_kirc_idx]
methotrexate_kirc_least_1se_auc <- auc(actual_methotrexate_kirc_least, methotrexate_kirc_least_sensitive_1se)

# LIHC
methotrexate_lihc_idx <- methotrexate_test$TCGA_class == 'LIHC' #7

methotrexate_lihc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_lihc_idx]
actual_methotrexate_lihc_most         <- methotrexate_test$most_sensitive[methotrexate_lihc_idx]
methotrexate_lihc_most_min_auc        <- auc(actual_methotrexate_lihc_most, methotrexate_lihc_most_sensitive_min)

methotrexate_lihc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_lihc_idx]
actual_methotrexate_lihc_least <- methotrexate_test$least_sensitive[methotrexate_lihc_idx]
methotrexate_lihc_least_min_auc <- auc(actual_methotrexate_lihc_least, methotrexate_lihc_least_sensitive_min)

methotrexate_lihc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_lihc_idx]
actual_methotrexate_lihc_most <- methotrexate_test$most_sensitive[methotrexate_lihc_idx]
methotrexate_lihc_most_1se_auc <- auc(actual_methotrexate_lihc_most, methotrexate_lihc_most_sensitive_1se)

methotrexate_lihc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_lihc_idx]
actual_methotrexate_lihc_least <- methotrexate_test$least_sensitive[methotrexate_lihc_idx]
methotrexate_lihc_least_1se_auc <- auc(actual_methotrexate_lihc_least, methotrexate_lihc_least_sensitive_1se)

# LUAD
methotrexate_luad_idx <- methotrexate_test$TCGA_class == 'LUAD' #26

methotrexate_luad_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_luad_idx]
actual_methotrexate_luad_most         <- methotrexate_test$most_sensitive[methotrexate_luad_idx]
methotrexate_luad_most_min_auc        <- auc(actual_methotrexate_luad_most, methotrexate_luad_most_sensitive_min)

methotrexate_luad_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_luad_idx]
actual_methotrexate_luad_least <- methotrexate_test$least_sensitive[methotrexate_luad_idx]
methotrexate_luad_least_min_auc <- auc(actual_methotrexate_luad_least, methotrexate_luad_least_sensitive_min)

methotrexate_luad_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_luad_idx]
actual_methotrexate_luad_most <- methotrexate_test$most_sensitive[methotrexate_luad_idx]
methotrexate_luad_most_1se_auc <- auc(actual_methotrexate_luad_most, methotrexate_luad_most_sensitive_1se)

methotrexate_luad_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_luad_idx]
actual_methotrexate_luad_least <- methotrexate_test$least_sensitive[methotrexate_luad_idx]
methotrexate_luad_least_1se_auc <- auc(actual_methotrexate_luad_least, methotrexate_luad_least_sensitive_1se)

# LUSC
methotrexate_lusc_idx <- methotrexate_test$TCGA_class == 'LUSC' #5

methotrexate_lusc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_lusc_idx]
actual_methotrexate_lusc_most         <- methotrexate_test$most_sensitive[methotrexate_lusc_idx]
methotrexate_lusc_most_min_auc        <- auc(actual_methotrexate_lusc_most, methotrexate_lusc_most_sensitive_min)

methotrexate_lusc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_lusc_idx]
actual_methotrexate_lusc_least <- methotrexate_test$least_sensitive[methotrexate_lusc_idx]
methotrexate_lusc_least_min_auc <- auc(actual_methotrexate_lusc_least, methotrexate_lusc_least_sensitive_min)

methotrexate_lusc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_lusc_idx]
actual_methotrexate_lusc_most <- methotrexate_test$most_sensitive[methotrexate_lusc_idx]
methotrexate_lusc_most_1se_auc <- auc(actual_methotrexate_lusc_most, methotrexate_lusc_most_sensitive_1se)

methotrexate_lusc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_lusc_idx]
actual_methotrexate_lusc_least <- methotrexate_test$least_sensitive[methotrexate_lusc_idx]
methotrexate_lusc_least_1se_auc <- auc(actual_methotrexate_lusc_least, methotrexate_lusc_least_sensitive_1se)

# MB
methotrexate_mb_idx <- methotrexate_test$TCGA_class == 'MB' #3

methotrexate_mb_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_mb_idx]
actual_methotrexate_mb_most         <- methotrexate_test$most_sensitive[methotrexate_mb_idx]
methotrexate_mb_most_min_auc        <- auc(actual_methotrexate_mb_most, methotrexate_mb_most_sensitive_min)

methotrexate_mb_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_mb_idx]
actual_methotrexate_mb_least <- methotrexate_test$least_sensitive[methotrexate_mb_idx]
methotrexate_mb_least_min_auc <- auc(actual_methotrexate_mb_least, methotrexate_mb_least_sensitive_min)

methotrexate_mb_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_mb_idx]
actual_methotrexate_mb_most <- methotrexate_test$most_sensitive[methotrexate_mb_idx]
methotrexate_mb_most_1se_auc <- auc(actual_methotrexate_mb_most, methotrexate_mb_most_sensitive_1se)

methotrexate_mb_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_mb_idx]
actual_methotrexate_mb_least <- methotrexate_test$least_sensitive[methotrexate_mb_idx]
methotrexate_mb_least_1se_auc <- auc(actual_methotrexate_mb_least, methotrexate_mb_least_sensitive_1se)

# MESO
methotrexate_meso_idx <- methotrexate_test$TCGA_class == 'MESO' #5

methotrexate_meso_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_meso_idx]
actual_methotrexate_meso_most         <- methotrexate_test$most_sensitive[methotrexate_meso_idx]
methotrexate_meso_most_min_auc        <- auc(actual_methotrexate_meso_most, methotrexate_meso_most_sensitive_min)

methotrexate_meso_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_meso_idx]
actual_methotrexate_meso_least <- methotrexate_test$least_sensitive[methotrexate_meso_idx]
methotrexate_meso_least_min_auc <- auc(actual_methotrexate_meso_least, methotrexate_meso_least_sensitive_min)

methotrexate_meso_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_meso_idx]
actual_methotrexate_meso_most <- methotrexate_test$most_sensitive[methotrexate_meso_idx]
methotrexate_meso_most_1se_auc <- auc(actual_methotrexate_meso_most, methotrexate_meso_most_sensitive_1se)

methotrexate_meso_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_meso_idx]
actual_methotrexate_meso_least <- methotrexate_test$least_sensitive[methotrexate_meso_idx]
methotrexate_meso_least_1se_auc <- auc(actual_methotrexate_meso_least, methotrexate_meso_least_sensitive_1se)

# NB
methotrexate_nb_idx <- methotrexate_test$TCGA_class == 'NB' #13

methotrexate_nb_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_nb_idx]
actual_methotrexate_nb_most         <- methotrexate_test$most_sensitive[methotrexate_nb_idx]
methotrexate_nb_most_min_auc        <- auc(actual_methotrexate_nb_most, methotrexate_nb_most_sensitive_min)

methotrexate_nb_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_nb_idx]
actual_methotrexate_nb_least <- methotrexate_test$least_sensitive[methotrexate_nb_idx]
methotrexate_nb_least_min_auc <- auc(actual_methotrexate_nb_least, methotrexate_nb_least_sensitive_min)

methotrexate_nb_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_nb_idx]
actual_methotrexate_nb_most <- methotrexate_test$most_sensitive[methotrexate_nb_idx]
methotrexate_nb_most_1se_auc <- auc(actual_methotrexate_nb_most, methotrexate_nb_most_sensitive_1se)

methotrexate_nb_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_nb_idx]
actual_methotrexate_nb_least <- methotrexate_test$least_sensitive[methotrexate_nb_idx]
methotrexate_nb_least_1se_auc <- auc(actual_methotrexate_nb_least, methotrexate_nb_least_sensitive_1se)

# OV
methotrexate_ov_idx <- methotrexate_test$TCGA_class == 'OV' #17

methotrexate_ov_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_ov_idx]
actual_methotrexate_ov_most         <- methotrexate_test$most_sensitive[methotrexate_ov_idx]
methotrexate_ov_most_min_auc        <- auc(actual_methotrexate_ov_most, methotrexate_ov_most_sensitive_min)

methotrexate_ov_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_ov_idx]
actual_methotrexate_ov_least <- methotrexate_test$least_sensitive[methotrexate_ov_idx]
methotrexate_ov_least_min_auc <- auc(actual_methotrexate_ov_least, methotrexate_ov_least_sensitive_min)

methotrexate_ov_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_ov_idx]
actual_methotrexate_ov_most <- methotrexate_test$most_sensitive[methotrexate_ov_idx]
methotrexate_ov_most_1se_auc <- auc(actual_methotrexate_ov_most, methotrexate_ov_most_sensitive_1se)

methotrexate_ov_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_ov_idx]
actual_methotrexate_ov_least <- methotrexate_test$least_sensitive[methotrexate_ov_idx]
methotrexate_ov_least_1se_auc <- auc(actual_methotrexate_ov_least, methotrexate_ov_least_sensitive_1se)

# PAAD
methotrexate_paad_idx <- methotrexate_test$TCGA_class == 'PAAD' #13

methotrexate_paad_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_paad_idx]
actual_methotrexate_paad_most         <- methotrexate_test$most_sensitive[methotrexate_paad_idx]
methotrexate_paad_most_min_auc        <- auc(actual_methotrexate_paad_most, methotrexate_paad_most_sensitive_min)

methotrexate_paad_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_paad_idx]
actual_methotrexate_paad_least <- methotrexate_test$least_sensitive[methotrexate_paad_idx]
methotrexate_paad_least_min_auc <- auc(actual_methotrexate_paad_least, methotrexate_paad_least_sensitive_min)

methotrexate_paad_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_paad_idx]
actual_methotrexate_paad_most <- methotrexate_test$most_sensitive[methotrexate_paad_idx]
methotrexate_paad_most_1se_auc <- auc(actual_methotrexate_paad_most, methotrexate_paad_most_sensitive_1se)

methotrexate_paad_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_paad_idx]
actual_methotrexate_paad_least <- methotrexate_test$least_sensitive[methotrexate_paad_idx]
methotrexate_paad_least_1se_auc <- auc(actual_methotrexate_paad_least, methotrexate_paad_least_sensitive_1se)

# SCLC
methotrexate_sclc_idx <- methotrexate_test$TCGA_class == 'SCLC' #31

methotrexate_sclc_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_sclc_idx]
actual_methotrexate_sclc_most         <- methotrexate_test$most_sensitive[methotrexate_sclc_idx]
methotrexate_sclc_most_min_auc        <- auc(actual_methotrexate_sclc_most, methotrexate_sclc_most_sensitive_min)

methotrexate_sclc_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_sclc_idx]
actual_methotrexate_sclc_least <- methotrexate_test$least_sensitive[methotrexate_sclc_idx]
methotrexate_sclc_least_min_auc <- auc(actual_methotrexate_sclc_least, methotrexate_sclc_least_sensitive_min)

methotrexate_sclc_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_sclc_idx]
actual_methotrexate_sclc_most <- methotrexate_test$most_sensitive[methotrexate_sclc_idx]
methotrexate_sclc_most_1se_auc <- auc(actual_methotrexate_sclc_most, methotrexate_sclc_most_sensitive_1se)

methotrexate_sclc_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_sclc_idx]
actual_methotrexate_sclc_least <- methotrexate_test$least_sensitive[methotrexate_sclc_idx]
methotrexate_sclc_least_1se_auc <- auc(actual_methotrexate_sclc_least, methotrexate_sclc_least_sensitive_1se)

# SKCM
methotrexate_skcm_idx <- methotrexate_test$TCGA_class == 'SKCM' #26

methotrexate_skcm_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_skcm_idx]
actual_methotrexate_skcm_most         <- methotrexate_test$most_sensitive[methotrexate_skcm_idx]
methotrexate_skcm_most_min_auc        <- auc(actual_methotrexate_skcm_most, methotrexate_skcm_most_sensitive_min)

methotrexate_skcm_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_skcm_idx]
actual_methotrexate_skcm_least <- methotrexate_test$least_sensitive[methotrexate_skcm_idx]
methotrexate_skcm_least_min_auc <- auc(actual_methotrexate_skcm_least, methotrexate_skcm_least_sensitive_min)

methotrexate_skcm_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_skcm_idx]
actual_methotrexate_skcm_most <- methotrexate_test$most_sensitive[methotrexate_skcm_idx]
methotrexate_skcm_most_1se_auc <- auc(actual_methotrexate_skcm_most, methotrexate_skcm_most_sensitive_1se)

methotrexate_skcm_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_skcm_idx]
actual_methotrexate_skcm_least <- methotrexate_test$least_sensitive[methotrexate_skcm_idx]
methotrexate_skcm_least_1se_auc <- auc(actual_methotrexate_skcm_least, methotrexate_skcm_least_sensitive_1se)

# STAD
methotrexate_stad_idx <- methotrexate_test$TCGA_class == 'STAD' #8

methotrexate_stad_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_stad_idx]
actual_methotrexate_stad_most         <- methotrexate_test$most_sensitive[methotrexate_stad_idx]
methotrexate_stad_most_min_auc        <- auc(actual_methotrexate_stad_most, methotrexate_stad_most_sensitive_min)

methotrexate_stad_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_stad_idx]
actual_methotrexate_stad_least <- methotrexate_test$least_sensitive[methotrexate_stad_idx]
methotrexate_stad_least_min_auc <- auc(actual_methotrexate_stad_least, methotrexate_stad_least_sensitive_min)

methotrexate_stad_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_stad_idx]
actual_methotrexate_stad_most <- methotrexate_test$most_sensitive[methotrexate_stad_idx]
methotrexate_stad_most_1se_auc <- auc(actual_methotrexate_stad_most, methotrexate_stad_most_sensitive_1se)

methotrexate_stad_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_stad_idx]
actual_methotrexate_stad_least <- methotrexate_test$least_sensitive[methotrexate_stad_idx]
methotrexate_stad_least_1se_auc <- auc(actual_methotrexate_stad_least, methotrexate_stad_least_sensitive_1se)

# UCEC
methotrexate_ucec_idx <- methotrexate_test$TCGA_class == 'UCEC' #4

methotrexate_ucec_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_ucec_idx]
actual_methotrexate_ucec_most         <- methotrexate_test$most_sensitive[methotrexate_ucec_idx]
methotrexate_ucec_most_min_auc        <- auc(actual_methotrexate_ucec_most, methotrexate_ucec_most_sensitive_min)

methotrexate_ucec_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_ucec_idx]
actual_methotrexate_ucec_least <- methotrexate_test$least_sensitive[methotrexate_ucec_idx]
methotrexate_ucec_least_min_auc <- auc(actual_methotrexate_ucec_least, methotrexate_ucec_least_sensitive_min)

methotrexate_ucec_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_ucec_idx]
actual_methotrexate_ucec_most <- methotrexate_test$most_sensitive[methotrexate_ucec_idx]
methotrexate_ucec_most_1se_auc <- auc(actual_methotrexate_ucec_most, methotrexate_ucec_most_sensitive_1se)

methotrexate_ucec_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_ucec_idx]
actual_methotrexate_ucec_least <- methotrexate_test$least_sensitive[methotrexate_ucec_idx]
methotrexate_ucec_least_1se_auc <- auc(actual_methotrexate_ucec_least, methotrexate_ucec_least_sensitive_1se)

# UNCLASSIFIED
methotrexate_unclassified_idx <- methotrexate_test$TCGA_class == 'UNCLASSIFIED' #52

methotrexate_unclassified_most_sensitive_min  <- new_methotrexate_most_sensitive_min[methotrexate_unclassified_idx]
actual_methotrexate_unclassified_most         <- methotrexate_test$most_sensitive[methotrexate_unclassified_idx]
methotrexate_unclassified_most_min_auc        <- auc(actual_methotrexate_unclassified_most, methotrexate_unclassified_most_sensitive_min)

methotrexate_unclassified_least_sensitive_min <- new_methotrexate_least_sensitive_min[methotrexate_unclassified_idx]
actual_methotrexate_unclassified_least <- methotrexate_test$least_sensitive[methotrexate_unclassified_idx]
methotrexate_unclassified_least_min_auc <- auc(actual_methotrexate_unclassified_least, methotrexate_unclassified_least_sensitive_min)

methotrexate_unclassified_most_sensitive_1se <- new_methotrexate_most_sensitive_1se[methotrexate_unclassified_idx]
actual_methotrexate_unclassified_most <- methotrexate_test$most_sensitive[methotrexate_unclassified_idx]
methotrexate_unclassified_most_1se_auc <- auc(actual_methotrexate_unclassified_most, methotrexate_unclassified_most_sensitive_1se)

methotrexate_unclassified_least_sensitive_1se <- new_methotrexate_least_sensitive_1se[methotrexate_unclassified_idx]
actual_methotrexate_unclassified_least <- methotrexate_test$least_sensitive[methotrexate_unclassified_idx]
methotrexate_unclassified_least_1se_auc <- auc(actual_methotrexate_unclassified_least, methotrexate_unclassified_least_sensitive_1se)

### put it together and graph -----
cisplatin_most_min_auc <- c(cisplatin_blca_most_min_auc, cisplatin_brca_most_min_auc, cisplatin_cesc_most_min_auc, 
                            cisplatin_coread_most_min_auc, cisplatin_esca_most_min_auc, cisplatin_hnsc_most_min_auc, 
                            cisplatin_kirc_most_min_auc, cisplatin_lihc_most_min_auc, cisplatin_luad_most_min_auc, 
                            cisplatin_lusc_most_min_auc, cisplatin_mb_most_min_auc, cisplatin_meso_most_min_auc, 
                            cisplatin_nb_most_min_auc, cisplatin_ov_most_min_auc, cisplatin_paad_most_min_auc, 
                            cisplatin_sclc_most_min_auc, cisplatin_skcm_most_min_auc, cisplatin_stad_most_min_auc, 
                            cisplatin_ucec_most_min_auc, cisplatin_unclassified_most_min_auc)

cisplatin_least_min_auc <- c(cisplatin_blca_least_min_auc, cisplatin_brca_least_min_auc, cisplatin_cesc_least_min_auc, 
                             cisplatin_coread_least_min_auc, cisplatin_esca_least_min_auc, cisplatin_hnsc_least_min_auc, 
                             cisplatin_kirc_least_min_auc, cisplatin_lihc_least_min_auc, cisplatin_luad_least_min_auc, 
                             cisplatin_lusc_least_min_auc, cisplatin_mb_least_min_auc, cisplatin_meso_least_min_auc, 
                             cisplatin_nb_least_min_auc, cisplatin_ov_least_min_auc, cisplatin_paad_least_min_auc, 
                             cisplatin_sclc_least_min_auc, cisplatin_skcm_least_min_auc, cisplatin_stad_least_min_auc, 
                             cisplatin_ucec_least_min_auc, cisplatin_unclassified_least_min_auc)

cisplatin_most_1se_auc <- c(cisplatin_blca_most_1se_auc, cisplatin_brca_most_1se_auc, cisplatin_cesc_most_1se_auc, 
                            cisplatin_coread_most_1se_auc, cisplatin_esca_most_1se_auc, cisplatin_hnsc_most_1se_auc, 
                            cisplatin_kirc_most_1se_auc, cisplatin_lihc_most_1se_auc, cisplatin_luad_most_1se_auc, 
                            cisplatin_lusc_most_1se_auc, cisplatin_mb_most_1se_auc, cisplatin_meso_most_1se_auc, 
                            cisplatin_nb_most_1se_auc, cisplatin_ov_most_1se_auc, cisplatin_paad_most_1se_auc, 
                            cisplatin_sclc_most_1se_auc, cisplatin_skcm_most_1se_auc, cisplatin_stad_most_1se_auc, 
                            cisplatin_ucec_most_1se_auc, cisplatin_unclassified_most_1se_auc)

cisplatin_least_1se_auc <- c(cisplatin_blca_least_1se_auc, cisplatin_brca_least_1se_auc, cisplatin_cesc_least_1se_auc, 
                             cisplatin_coread_least_1se_auc, cisplatin_esca_least_1se_auc, cisplatin_hnsc_least_1se_auc, 
                             cisplatin_kirc_least_1se_auc, cisplatin_lihc_least_1se_auc, cisplatin_luad_least_1se_auc, 
                             cisplatin_lusc_least_1se_auc, cisplatin_mb_least_1se_auc, cisplatin_meso_least_1se_auc, 
                             cisplatin_nb_least_1se_auc, cisplatin_ov_least_1se_auc, cisplatin_paad_least_1se_auc, 
                             cisplatin_sclc_least_1se_auc, cisplatin_skcm_least_1se_auc, cisplatin_stad_least_1se_auc, 
                             cisplatin_ucec_least_1se_auc, cisplatin_unclassified_least_1se_auc)

etoposide_most_min_auc <- c(etoposide_blca_most_min_auc, etoposide_brca_most_min_auc, etoposide_cesc_most_min_auc, 
                            etoposide_coread_most_min_auc, etoposide_esca_most_min_auc, etoposide_hnsc_most_min_auc, 
                            etoposide_kirc_most_min_auc, etoposide_lihc_most_min_auc, etoposide_luad_most_min_auc, 
                            etoposide_lusc_most_min_auc, etoposide_mb_most_min_auc, etoposide_meso_most_min_auc, 
                            etoposide_nb_most_min_auc, etoposide_ov_most_min_auc, etoposide_paad_most_min_auc, 
                            etoposide_sclc_most_min_auc, etoposide_skcm_most_min_auc, etoposide_stad_most_min_auc, 
                            etoposide_ucec_most_min_auc, etoposide_unclassified_most_min_auc)

etoposide_least_min_auc <- c(etoposide_blca_least_min_auc, etoposide_brca_least_min_auc, etoposide_cesc_least_min_auc, 
                             etoposide_coread_least_min_auc, etoposide_esca_least_min_auc, etoposide_hnsc_least_min_auc, 
                             etoposide_kirc_least_min_auc, etoposide_lihc_least_min_auc, etoposide_luad_least_min_auc, 
                             etoposide_lusc_least_min_auc, etoposide_mb_least_min_auc, etoposide_meso_least_min_auc, 
                             etoposide_nb_least_min_auc, etoposide_ov_least_min_auc, etoposide_paad_least_min_auc, 
                             etoposide_sclc_least_min_auc, etoposide_skcm_least_min_auc, etoposide_stad_least_min_auc, 
                             etoposide_ucec_least_min_auc, etoposide_unclassified_least_min_auc)

etoposide_most_1se_auc <- c(etoposide_blca_most_1se_auc, etoposide_brca_most_1se_auc, etoposide_cesc_most_1se_auc, 
                            etoposide_coread_most_1se_auc, etoposide_esca_most_1se_auc, etoposide_hnsc_most_1se_auc, 
                            etoposide_kirc_most_1se_auc, etoposide_lihc_most_1se_auc, etoposide_luad_most_1se_auc, 
                            etoposide_lusc_most_1se_auc, etoposide_mb_most_1se_auc, etoposide_meso_most_1se_auc, 
                            etoposide_nb_most_1se_auc, etoposide_ov_most_1se_auc, etoposide_paad_most_1se_auc, 
                            etoposide_sclc_most_1se_auc, etoposide_skcm_most_1se_auc, etoposide_stad_most_1se_auc, 
                            etoposide_ucec_most_1se_auc, etoposide_unclassified_most_1se_auc)

etoposide_least_1se_auc <- c(etoposide_blca_least_1se_auc, etoposide_brca_least_1se_auc, etoposide_cesc_least_1se_auc, 
                             etoposide_coread_least_1se_auc, etoposide_esca_least_1se_auc, etoposide_hnsc_least_1se_auc, 
                             etoposide_kirc_least_1se_auc, etoposide_lihc_least_1se_auc, etoposide_luad_least_1se_auc, 
                             etoposide_lusc_least_1se_auc, etoposide_mb_least_1se_auc, etoposide_meso_least_1se_auc, 
                             etoposide_nb_least_1se_auc, etoposide_ov_least_1se_auc, etoposide_paad_least_1se_auc, 
                             etoposide_sclc_least_1se_auc, etoposide_skcm_least_1se_auc, etoposide_stad_least_1se_auc, 
                             etoposide_ucec_least_1se_auc, etoposide_unclassified_least_1se_auc)

gemcitabine_most_min_auc <- c(gemcitabine_blca_most_min_auc, gemcitabine_brca_most_min_auc, gemcitabine_cesc_most_min_auc, 
                            gemcitabine_coread_most_min_auc, gemcitabine_esca_most_min_auc, gemcitabine_hnsc_most_min_auc, 
                            gemcitabine_kirc_most_min_auc, gemcitabine_lihc_most_min_auc, gemcitabine_luad_most_min_auc, 
                            gemcitabine_lusc_most_min_auc, gemcitabine_mb_most_min_auc, gemcitabine_meso_most_min_auc, 
                            gemcitabine_nb_most_min_auc, gemcitabine_ov_most_min_auc, gemcitabine_paad_most_min_auc, 
                            gemcitabine_sclc_most_min_auc, gemcitabine_skcm_most_min_auc, gemcitabine_stad_most_min_auc, 
                            gemcitabine_ucec_most_min_auc, gemcitabine_unclassified_most_min_auc)

gemcitabine_least_min_auc <- c(gemcitabine_blca_least_min_auc, gemcitabine_brca_least_min_auc, gemcitabine_cesc_least_min_auc, 
                             gemcitabine_coread_least_min_auc, gemcitabine_esca_least_min_auc, gemcitabine_hnsc_least_min_auc, 
                             gemcitabine_kirc_least_min_auc, gemcitabine_lihc_least_min_auc, gemcitabine_luad_least_min_auc, 
                             gemcitabine_lusc_least_min_auc, gemcitabine_mb_least_min_auc, gemcitabine_meso_least_min_auc, 
                             gemcitabine_nb_least_min_auc, gemcitabine_ov_least_min_auc, gemcitabine_paad_least_min_auc, 
                             gemcitabine_sclc_least_min_auc, gemcitabine_skcm_least_min_auc, gemcitabine_stad_least_min_auc, 
                             gemcitabine_ucec_least_min_auc, gemcitabine_unclassified_least_min_auc)

gemcitabine_most_1se_auc <- c(gemcitabine_blca_most_1se_auc, gemcitabine_brca_most_1se_auc, gemcitabine_cesc_most_1se_auc, 
                            gemcitabine_coread_most_1se_auc, gemcitabine_esca_most_1se_auc, gemcitabine_hnsc_most_1se_auc, 
                            gemcitabine_kirc_most_1se_auc, gemcitabine_lihc_most_1se_auc, gemcitabine_luad_most_1se_auc, 
                            gemcitabine_lusc_most_1se_auc, gemcitabine_mb_most_1se_auc, gemcitabine_meso_most_1se_auc, 
                            gemcitabine_nb_most_1se_auc, gemcitabine_ov_most_1se_auc, gemcitabine_paad_most_1se_auc, 
                            gemcitabine_sclc_most_1se_auc, gemcitabine_skcm_most_1se_auc, gemcitabine_stad_most_1se_auc, 
                            gemcitabine_ucec_most_1se_auc, gemcitabine_unclassified_most_1se_auc)

gemcitabine_least_1se_auc <- c(gemcitabine_blca_least_1se_auc, gemcitabine_brca_least_1se_auc, gemcitabine_cesc_least_1se_auc, 
                             gemcitabine_coread_least_1se_auc, gemcitabine_esca_least_1se_auc, gemcitabine_hnsc_least_1se_auc, 
                             gemcitabine_kirc_least_1se_auc, gemcitabine_lihc_least_1se_auc, gemcitabine_luad_least_1se_auc, 
                             gemcitabine_lusc_least_1se_auc, gemcitabine_mb_least_1se_auc, gemcitabine_meso_least_1se_auc, 
                             gemcitabine_nb_least_1se_auc, gemcitabine_ov_least_1se_auc, gemcitabine_paad_least_1se_auc, 
                             gemcitabine_sclc_least_1se_auc, gemcitabine_skcm_least_1se_auc, gemcitabine_stad_least_1se_auc, 
                             gemcitabine_ucec_least_1se_auc, gemcitabine_unclassified_least_1se_auc)

methotrexate_most_min_auc <- c(methotrexate_blca_most_min_auc, methotrexate_brca_most_min_auc, methotrexate_cesc_most_min_auc, 
                            methotrexate_coread_most_min_auc, methotrexate_esca_most_min_auc, methotrexate_hnsc_most_min_auc, 
                            methotrexate_kirc_most_min_auc, methotrexate_lihc_most_min_auc, methotrexate_luad_most_min_auc, 
                            methotrexate_lusc_most_min_auc, methotrexate_mb_most_min_auc, methotrexate_meso_most_min_auc, 
                            methotrexate_nb_most_min_auc, methotrexate_ov_most_min_auc, methotrexate_paad_most_min_auc, 
                            methotrexate_sclc_most_min_auc, methotrexate_skcm_most_min_auc, methotrexate_stad_most_min_auc, 
                            methotrexate_ucec_most_min_auc, methotrexate_unclassified_most_min_auc)

methotrexate_least_min_auc <- c(methotrexate_blca_least_min_auc, methotrexate_brca_least_min_auc, methotrexate_cesc_least_min_auc, 
                             methotrexate_coread_least_min_auc, methotrexate_esca_least_min_auc, methotrexate_hnsc_least_min_auc, 
                             methotrexate_kirc_least_min_auc, methotrexate_lihc_least_min_auc, methotrexate_luad_least_min_auc, 
                             methotrexate_lusc_least_min_auc, methotrexate_mb_least_min_auc, methotrexate_meso_least_min_auc, 
                             methotrexate_nb_least_min_auc, methotrexate_ov_least_min_auc, methotrexate_paad_least_min_auc, 
                             methotrexate_sclc_least_min_auc, methotrexate_skcm_least_min_auc, methotrexate_stad_least_min_auc, 
                             methotrexate_ucec_least_min_auc, methotrexate_unclassified_least_min_auc)

methotrexate_most_1se_auc <- c(methotrexate_blca_most_1se_auc, methotrexate_brca_most_1se_auc, methotrexate_cesc_most_1se_auc, 
                            methotrexate_coread_most_1se_auc, methotrexate_esca_most_1se_auc, methotrexate_hnsc_most_1se_auc, 
                            methotrexate_kirc_most_1se_auc, methotrexate_lihc_most_1se_auc, methotrexate_luad_most_1se_auc, 
                            methotrexate_lusc_most_1se_auc, methotrexate_mb_most_1se_auc, methotrexate_meso_most_1se_auc, 
                            methotrexate_nb_most_1se_auc, methotrexate_ov_most_1se_auc, methotrexate_paad_most_1se_auc, 
                            methotrexate_sclc_most_1se_auc, methotrexate_skcm_most_1se_auc, methotrexate_stad_most_1se_auc, 
                            methotrexate_ucec_most_1se_auc, methotrexate_unclassified_most_1se_auc)

methotrexate_least_1se_auc <- c(methotrexate_blca_least_1se_auc, methotrexate_brca_least_1se_auc, methotrexate_cesc_least_1se_auc, 
                             methotrexate_coread_least_1se_auc, methotrexate_esca_least_1se_auc, methotrexate_hnsc_least_1se_auc, 
                             methotrexate_kirc_least_1se_auc, methotrexate_lihc_least_1se_auc, methotrexate_luad_least_1se_auc, 
                             methotrexate_lusc_least_1se_auc, methotrexate_mb_least_1se_auc, methotrexate_meso_least_1se_auc, 
                             methotrexate_nb_least_1se_auc, methotrexate_ov_least_1se_auc, methotrexate_paad_least_1se_auc, 
                             methotrexate_sclc_least_1se_auc, methotrexate_skcm_least_1se_auc, methotrexate_stad_least_1se_auc, 
                             methotrexate_ucec_least_1se_auc, methotrexate_unclassified_least_1se_auc)


all_models_auc_cancer_type <- data.frame(cisplatin_most_min_auc, cisplatin_most_1se_auc, 
                                         cisplatin_least_min_auc, cisplatin_least_1se_auc,
                                         etoposide_most_min_auc, etoposide_most_1se_auc, 
                                         etoposide_least_min_auc, etoposide_least_1se_auc, 
                                         gemcitabine_most_min_auc, gemcitabine_most_1se_auc, 
                                         gemcitabine_least_min_auc, gemcitabine_least_1se_auc, 
                                         methotrexate_most_min_auc, methotrexate_most_1se_auc, 
                                         methotrexate_least_min_auc, methotrexate_least_1se_auc)
for (i in 1:nrow(all_models_auc_cancer_type)) {
  for (j in 1:ncol(all_models_auc_cancer_type)) {
    if (all_models_auc_cancer_type[i,j] == 'NaN') {
      all_models_auc_cancer_type[i,j] <- 0.5
    }
  }
}

for (i in 1:nrow(all_models_auc_cancer_type)) {
  for (j in 1:ncol(all_models_auc_cancer_type)) {
    if (all_models_auc_cancer_type[i,j] < 0.5) {
      all_models_auc_cancer_type[i,j] <- 0.5
    }
  }
}

rownames(all_models_auc_cancer_type) <- c('BLCA', 'BRCA', 'CESC', 'COREAD', 'ESCA', 'HNSC', 
                                          'KIRC', 'LIHC', 'LUAD', 'LUSC', 'MB', 'MESO', 
                                          'NB', 'OV', 'PAAD', 'SCLC', 'SKCM', 'STAD', 
                                          'UCEC', 'UNCLASSIFIED')
all_models_auc_cancer_type <- round(all_models_auc_cancer_type, digits = 2)
all_mat <- as.matrix(all_models_auc_cancer_type)
colors <- colorRampPalette(c('blue', 'white', 'red'))(100)
heatmap.2(all_mat, trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, key = FALSE, 
          cexRow = 0.9, cexCol = 0.7, margins = c(8,8))

### graphing auc curves for gdsc testing results ---------
#pred <- prediction(new_bleo_most_sensitive_min[bleo_unclassified_idx], bleomycin_test$most_sensitive[bleo_unclassified_idx])
#pred2 <- prediction(new_bleo_most_sensitive_1se[bleo_unclassified_idx], bleomycin_test$most_sensitive[bleo_unclassified_idx])
#pred3 <- prediction(new_bleo_least_sensitive_min[bleo_unclassified_idx], bleomycin_test$least_sensitive[bleo_unclassified_idx])
#pred4 <- prediction(new_bleo_least_sensitive_1se[bleo_unclassified_idx], bleomycin_test$least_sensitive[bleo_unclassified_idx])

#perf <- performance(pred, 'tpr', 'fpr')
#perf2 <- performance(pred2, 'tpr', 'fpr')
#perf3 <- performance(pred3, 'tpr', 'fpr')
#perf4 <- performance(pred4, 'tpr', 'fpr')

#png(file = '')
#plot(perf, col = 'red')
#plot(perf2, col = 'blue', add = TRUE)
#plot(perf3, col = 'red', lty = 2, add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
#abline(0,1, lty = 3)
#auc <- performance(pred4, 'auc')
#auc <- unlist(slot(auc, 'y.values'))
#minauc <- min(round(auc, digits = 2))
#maxauc <- max(round(auc, digits = 2))
#minauct <- paste(c('minAUC = '), minauc, sep = '')
#maxauct <- paste(c("max(AUC) = "),maxauc,sep="")
#legend(0.7,0.6,c('most_min (0.76)', 'most_1se (0.75)', 'least_min (0.92)', 'least_1se (0.90)'), lty = c(1,1,2,2), col = c('red', 'blue', 'red', 'blue'), border="white",cex=0.9,box.col = "white")
#dev.off()


### testing on TCGA classes treated with drug ----
# BLCA W CISPLATIN (49)
blca_clinical <- read.csv('Processed_Clinical_Data/blca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(blca_clinical$most_sensitive)
blca_clinical <- blca_clinical[!na_idx, ]
table(blca_clinical$drug_name)
blca_clinical_cisplatin <- blca_clinical[which(blca_clinical$drug_name == 'cisplatin' | blca_clinical$drug_name == 'Cisplatin' | 
                                                 blca_clinical$drug_name == 'Cisplatnin'), ]

blca_clinical_cisplatin$most_sensitive  <- ifelse(blca_clinical_cisplatin$PFS < quantile(blca_clinical_cisplatin$PFS, probs = .20), 1, 0)
blca_clinical_cisplatin$least_sensitive <- ifelse(blca_clinical_cisplatin$PFS > quantile(blca_clinical_cisplatin$PFS, probs = .80), 1, 0)

blca_gene <- read.csv('Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(blca_gene) <- gsub('\\.', '-', colnames(blca_gene))
blca_matching_idx <- blca_clinical_cisplatin$submitter_id.samples %in% colnames(blca_gene)
blca_clinical_cisplatin_short <- blca_clinical_cisplatin[blca_matching_idx, ]
blca_matching_idx <- colnames(blca_gene) %in% blca_clinical_cisplatin_short$submitter_id.samples
blca_gene_short <- blca_gene[, blca_matching_idx]
blca_gene_short <- t(blca_gene_short)
blca_gene_short_scaled <- apply(blca_gene_short, 2, scale)

rm(blca_tcga_most_min_auc)
rm(blca_tcga_most_1se_auc)
rm(blca_tcga_least_min_auc)
rm(blca_tcga_least_1se_auc)

new_blca_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_most_min_auc <- auc(blca_clinical_cisplatin_short$most_sensitive, new_blca_tcga_cisplatin_most_sensitive_min)
blca_tcga_most_min_auc <- round(blca_tcga_most_min_auc, digits = 2)

new_blca_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_most_1se_auc <- auc(blca_clinical_cisplatin_short$most_sensitive, new_blca_tcga_cisplatin_most_sensitive_1se)
blca_tcga_most_1se_auc <- round(blca_tcga_most_1se_auc, digits = 2)

new_blca_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_least_min_auc <- auc(blca_clinical_cisplatin_short$least_sensitive, new_blca_tcga_cisplatin_least_sensitive_min)
blca_tcga_least_min_auc <- round(blca_tcga_least_min_auc, digits = 2)

new_blca_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_least_1se_auc <- auc(blca_clinical_cisplatin_short$least_sensitive, new_blca_tcga_cisplatin_least_sensitive_1se)
blca_tcga_least_1se_auc <- round(blca_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_blca_tcga_cisplatin_most_sensitive_min, blca_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_blca_tcga_cisplatin_most_sensitive_1se, blca_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_blca_tcga_cisplatin_least_sensitive_min, blca_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_blca_tcga_cisplatin_least_sensitive_1se, blca_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/blca_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 1, main = 'BLCA treated with cisplatin')
plot(perf3, col = 'blue', add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.55)', 'least_sensitive (AUC = 0.69)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
dev.off()

# BLCA W GEMCITABINE (30)
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

new_blca_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_most_min_auc <- auc(blca_clinical_gemcitabine_short$most_sensitive, new_blca_tcga_gemcitabine_most_sensitive_min)
blca_tcga_most_min_auc <- round(blca_tcga_most_min_auc, digits = 2)

new_blca_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
blca_tcga_most_1se_auc <- auc(blca_clinical_gemcitabine_short$most_sensitive, new_blca_tcga_gemcitabine_most_sensitive_1se)
blca_tcga_most_1se_auc <- round(blca_tcga_most_1se_auc, digits = 2)

new_blca_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
blca_tcga_least_min_auc <- auc(blca_clinical_gemcitabine_short$least_sensitive, new_blca_tcga_gemcitabine_least_sensitive_min)
blca_tcga_least_min_auc <- round(blca_tcga_least_min_auc, digits = 2)

new_blca_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = blca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
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
png(filename = 'Images/blca_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 1, main = 'BLCA treated with gemcitabine (GDSC)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.72)', 'least_sensitive (AUC = 0.56)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
dev.off()

# BRCA W DOXORUBICIN (47)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_doxorubicin <- brca_clinical[which(brca_clinical$drug_name == 'doxorubicin' | brca_clinical$drug_name == 'Doxorubicin' | 
                                                 brca_clinical$drug_name == 'DOXORUBICIN'), ]

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

new_brca_tcga_doxorubicin_most_sensitive_min <- predict(dox_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_doxorubicin_short$most_sensitive, new_brca_tcga_doxorubicin_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_doxorubicin_most_sensitive_1se <- predict(dox_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_doxorubicin_short$most_sensitive, new_brca_tcga_doxorubicin_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_doxorubicin_least_sensitive_min <- predict(dox_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_doxorubicin_short$least_sensitive, new_brca_tcga_doxorubicin_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_doxorubicin_least_sensitive_1se <- predict(dox_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
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
png(filename = 'Images/brca_doxorubicin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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

new_brca_tcga_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_gemcitabine_short$most_sensitive, new_brca_tcga_gemcitabine_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_gemcitabine_short$most_sensitive, new_brca_tcga_gemcitabine_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_gemcitabine_short$least_sensitive, new_brca_tcga_gemcitabine_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
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
png(filename = 'Images/brca_gemcitabine_tcga_gdsc_auc.png', height = 480, width = 800)
#plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 1, main = 'BRCA treated with gemcitabine (GDSC)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 1.0)', 'least_sensitive (AUC = 1.0)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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

new_brca_tcga_methotrexate_most_sensitive_min <- predict(methotrexate_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_most_min_auc <- auc(brca_clinical_methotrexate_short$most_sensitive, new_brca_tcga_methotrexate_most_sensitive_min)
brca_tcga_most_min_auc <- round(brca_tcga_most_min_auc, digits = 2)

new_brca_tcga_methotrexate_most_sensitive_1se <- predict(methotrexate_most_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
brca_tcga_most_1se_auc <- auc(brca_clinical_methotrexate_short$most_sensitive, new_brca_tcga_methotrexate_most_sensitive_1se)
brca_tcga_most_1se_auc <- round(brca_tcga_most_1se_auc, digits = 2)

new_brca_tcga_methotrexate_least_sensitive_min <- predict(methotrexate_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
brca_tcga_least_min_auc <- auc(brca_clinical_methotrexate_short$least_sensitive, new_brca_tcga_methotrexate_least_sensitive_min)
brca_tcga_least_min_auc <- round(brca_tcga_least_min_auc, digits = 2)

new_brca_tcga_methotrexate_least_sensitive_1se <- predict(methotrexate_least_fit_elnet, newx = brca_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
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
png(filename = 'Images/brca_methotrexate_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red', main = 'BRCA treated with methotrexate (GDSC)')
#plot(perf2, col = 'red', lty = 2)
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 1.0)', 'least_sensitive (AUC = 0.64)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
dev.off()

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

rm(cesc_tcga_most_min_auc)
rm(cesc_tcga_most_1se_auc)
rm(cesc_tcga_least_min_auc)
rm(cesc_tcga_least_1se_auc)

new_cesc_tcga_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_most_min_auc <- auc(cesc_clinical_cisplatin_short$most_sensitive, new_cesc_tcga_cisplatin_most_sensitive_min)
cesc_tcga_most_min_auc <- round(cesc_tcga_most_min_auc, digits = 2)

new_cesc_tcga_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_most_1se_auc <- auc(cesc_clinical_cisplatin_short$most_sensitive, new_cesc_tcga_cisplatin_most_sensitive_1se)
cesc_tcga_most_1se_auc <- round(cesc_tcga_most_1se_auc, digits = 2)

new_cesc_tcga_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.min', type = 'response', na.action = 'na.pass')
cesc_tcga_least_min_auc <- auc(cesc_clinical_cisplatin_short$least_sensitive, new_cesc_tcga_cisplatin_least_sensitive_min)
cesc_tcga_least_min_auc <- round(cesc_tcga_least_min_auc, digits = 2)

new_cesc_tcga_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = cesc_gene_short_scaled, s = 'lambda.1se', type = 'response', na.action = 'na.pass')
cesc_tcga_least_1se_auc <- auc(cesc_clinical_cisplatin_short$least_sensitive, new_cesc_tcga_cisplatin_least_sensitive_1se)
cesc_tcga_least_1se_auc <- round(cesc_tcga_least_1se_auc, digits = 2)

rm(pred)
rm(pred2)
rm(pred3)
rm(pred4)

pred <- prediction(new_cesc_tcga_cisplatin_most_sensitive_min, cesc_clinical_cisplatin_short$most_sensitive)
pred2 <- prediction(new_cesc_tcga_cisplatin_most_sensitive_1se, cesc_clinical_cisplatin_short$most_sensitive)
pred3 <- prediction(new_cesc_tcga_cisplatin_least_sensitive_min, cesc_clinical_cisplatin_short$least_sensitive)
pred4 <- prediction(new_cesc_tcga_cisplatin_least_sensitive_1se, cesc_clinical_cisplatin_short$least_sensitive)

rm(perf)
rm(perf2)
rm(perf3)
rm(perf4)

perf <- performance(pred, 'tpr', 'fpr')
perf2 <- performance(pred2, 'tpr', 'fpr')
perf3 <- performance(pred3, 'tpr', 'fpr')
perf4 <- performance(pred4, 'tpr', 'fpr')

# change these to have the order make sense
png(filename = 'Images/cesc_cisplatin_tcga_gdsc_auc.png', height = 480, width = 800)
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'ESCA treated with cisplatin')
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf, col = 'red', main = 'LUAD treated with etoposide (GDSC)')
plot(perf2, col = 'red', lty = 1, main = 'LUAD treated with etoposide (GDSC)')
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (NA)', 'most_1se (1.0)', 'least_min (0.00)', 'least_1se (0.00)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'OV treated with cisplatin')
plot(perf3, col = 'blue', add = TRUE)
#plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.83)', 'least_sensitive (AUC = 0.80)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
etoposide_ccle <- read.csv('Processed_Clinical_Data/etoposide_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
fluorouracil_ccle <- read.csv('Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
gemcitabine_ccle <- read.csv('Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
methotrexate_ccle <- read.csv('Processed_Clinical_Data/methotrexate_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)

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
etoposide_lines <- etoposide_ccle$Cell.line.name #649
etoposide_lines <- as.character(etoposide_lines)
fluorouracil_lines <- fluorouracil_ccle$Cell.line.name #652
fluorouracil_lines <- as.character(fluorouracil_lines)
gemcitabine_lines <- gemcitabine_ccle$Cell.line.name #589
gemcitabine_lines <- as.character(gemcitabine_lines)
methotrexate_lines <- methotrexate_ccle$Cell.line.name #622
methotrexate_lines <- as.character(methotrexate_lines)


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

etoposide_rna_seq_train     <- ccle_train[intersect(etoposide_lines, rownames(ccle_train)), ]
# 310 x 14209
etoposide_rna_seq_test      <- ccle_test[intersect(etoposide_lines, rownames(ccle_test)), ]
# 307 x 14209

fluorouracil_rna_seq_train     <- ccle_train[intersect(fluorouracil_lines, rownames(ccle_train)), ]
# 313 x 14209
fluorouracil_rna_seq_test      <- ccle_test[intersect(fluorouracil_lines, rownames(ccle_test)), ]
# 307 x 14209

gemcitabine_rna_seq_train     <- ccle_train[intersect(gemcitabine_lines, rownames(ccle_train)), ]
# 276 x 14209
gemcitabine_rna_seq_test      <- ccle_test[intersect(gemcitabine_lines, rownames(ccle_test)), ]
# 286 x 14209

methotrexate_rna_seq_train     <- ccle_train[intersect(methotrexate_lines, rownames(ccle_train)), ]
# 299 x 14209
methotrexate_rna_seq_test      <- ccle_test[intersect(methotrexate_lines, rownames(ccle_test)), ]
# 295 x 14209


carboplatin_train        <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_train)), ]
carboplatin_test         <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_test)), ]

cyclophosphamide_train        <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_train)), ]
cyclophosphamide_test         <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_test)), ]

docetaxel_train        <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_train)), ]
docetaxel_test         <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_test)), ]

etoposide_train        <- etoposide_ccle[which(etoposide_ccle$Cell.line.name %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide_ccle[which(etoposide_ccle$Cell.line.name %in% rownames(etoposide_rna_seq_test)), ]

fluorouracil_train        <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_train)), ]
fluorouracil_test         <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_test)), ]

gemcitabine_train        <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test         <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_test)), ]

methotrexate_train        <- methotrexate_ccle[which(methotrexate_ccle$Cell.line.name %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test         <- methotrexate_ccle[which(methotrexate_ccle$Cell.line.name %in% rownames(methotrexate_rna_seq_test)), ]

carboplatin_rna_seq_train_scaled          <- apply(carboplatin_rna_seq_train, 2, scale)
carboplatin_rna_seq_test_scaled           <- as.data.frame(apply(carboplatin_rna_seq_test, 2, scale))

cyclophosphamide_rna_seq_train_scaled          <- apply(cyclophosphamide_rna_seq_train, 2, scale)
cyclophosphamide_rna_seq_test_scaled           <- as.data.frame(apply(cyclophosphamide_rna_seq_test, 2, scale))

docetaxel_rna_seq_train_scaled          <- apply(docetaxel_rna_seq_train, 2, scale)
docetaxel_rna_seq_test_scaled           <- as.data.frame(apply(docetaxel_rna_seq_test, 2, scale))

etoposide_rna_seq_train_scaled          <- apply(etoposide_rna_seq_train, 2, scale)
etoposide_rna_seq_test_scaled           <- as.data.frame(apply(etoposide_rna_seq_test, 2, scale))

fluorouracil_rna_seq_train_scaled          <- apply(fluorouracil_rna_seq_train, 2, scale)
fluorouracil_rna_seq_test_scaled           <- as.data.frame(apply(fluorouracil_rna_seq_test, 2, scale))

gemcitabine_rna_seq_train_scaled          <- apply(gemcitabine_rna_seq_train, 2, scale)
gemcitabine_rna_seq_test_scaled           <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))

methotrexate_rna_seq_train_scaled          <- apply(methotrexate_rna_seq_train, 2, scale)
methotrexate_rna_seq_test_scaled           <- as.data.frame(apply(methotrexate_rna_seq_test, 2, scale))


### load models ----
carboplatin_ccle_most_fit_elnet      <- readRDS('GLM_Models/carboplatin_ccle_most_model.rds')
carboplatin_ccle_least_fit_elnet     <- readRDS('GLM_Models/carboplatin_ccle_least_model.rds')

cyclophosphamide_ccle_most_fit_elnet      <- readRDS('GLM_Models/cyclophosphamide_ccle_most_model.rds')
cyclophosphamide_ccle_least_fit_elnet     <- readRDS('GLM_Models/cyclophosphamide_ccle_least_model.rds')

docetaxel_ccle_most_fit_elnet   <- readRDS('GLM_Models/docetaxel_ccle_most_model.rds')
docetaxel_ccle_least_fit_elnet  <- readRDS('GLM_Models/docetaxel_ccle_least_model.rds')

etoposide_ccle_most_fit_elnet <- readRDS('GLM_Models/etoposide_ccle_most_model.rds')
etoposide_ccle_least_fit_elnet <- readRDS('GLM_Models/etoposide_ccle_least_model.rds')

fluorouracil_ccle_most_fit_elnet <- readRDS('GLM_Models/fluorouracil_ccle_most_model.rds')
fluorouracil_ccle_least_fit_elnet <- readRDS('GLM_Models/fluorouracil_ccle_least_model.rds')

gemcitabine_ccle_most_fit_elnet <- readRDS('GLM_Models/gemcitabine_ccle_most_model.rds')
gemcitabine_ccle_least_fit_elnet <- readRDS('GLM_Models/gemcitabine_ccle_least_model.rds')

methotrexate_ccle_most_fit_elnet <- readRDS('GLM_Models/methotrexate_ccle_most_model.rds')
methotrexate_ccle_least_fit_elnet <- readRDS('GLM_Models/methotrexate_ccle_least_model.rds')


### get accuracy on testing data --------

## CARBOPLATIN
#new_carboplatin_most_sensitive_min <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#carboplatin_most_test_gdsc_auc_min <- auc(carboplatin_test$most_sensitive, new_carboplatin_most_sensitive_min)
#carboplatin_most_test_gdsc_auc_min <- round(carboplatin_most_test_gdsc_auc_min, digits = 2)

new_carboplatin_most_sensitive_1se <- predict(carboplatin_ccle_most_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

carboplatin_most_test_gdsc_auc_1se <- auc(carboplatin_test$most_sensitive, new_carboplatin_most_sensitive_1se)
carboplatin_most_test_gdsc_auc_1se <- round(carboplatin_most_test_gdsc_auc_1se, digits = 2)


#new_carboplatin_least_sensitive_min <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#carboplatin_least_test_gdsc_auc_min <- auc(carboplatin_test$least_sensitive, new_carboplatin_least_sensitive_min)
#carboplatin_least_test_gdsc_auc_min <- round(carboplatin_least_test_gdsc_auc_min, digits = 2)

new_carboplatin_least_sensitive_1se <- predict(carboplatin_ccle_least_fit_elnet, newx = as.matrix(carboplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

carboplatin_least_test_gdsc_auc_1se <- auc(carboplatin_test$least_sensitive, new_carboplatin_least_sensitive_1se)
carboplatin_least_test_gdsc_auc_1se <- round(carboplatin_least_test_gdsc_auc_1se, digits = 2)

## cyclophosphamide
#new_cyclophosphamide_most_sensitive_min <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#cyclophosphamide_most_test_gdsc_auc_min <- auc(cyclophosphamide_test$most_sensitive, new_cyclophosphamide_most_sensitive_min)
#cyclophosphamide_most_test_gdsc_auc_min <- round(cyclophosphamide_most_test_gdsc_auc_min, digits = 2)

new_cyclophosphamide_most_sensitive_1se <- predict(cyclophosphamide_ccle_most_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cyclophosphamide_most_test_gdsc_auc_1se <- auc(cyclophosphamide_test$most_sensitive, new_cyclophosphamide_most_sensitive_1se)
cyclophosphamide_most_test_gdsc_auc_1se <- round(cyclophosphamide_most_test_gdsc_auc_1se, digits = 2)


#new_cyclophosphamide_least_sensitive_min <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#cyclophosphamide_least_test_gdsc_auc_min <- auc(cyclophosphamide_test$least_sensitive, new_cyclophosphamide_least_sensitive_min)
#cyclophosphamide_least_test_gdsc_auc_min <- round(cyclophosphamide_least_test_gdsc_auc_min, digits = 2)

new_cyclophosphamide_least_sensitive_1se <- predict(cyclophosphamide_ccle_least_fit_elnet, newx = as.matrix(cyclophosphamide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cyclophosphamide_least_test_gdsc_auc_1se <- auc(cyclophosphamide_test$least_sensitive, new_cyclophosphamide_least_sensitive_1se)
cyclophosphamide_least_test_gdsc_auc_1se <- round(cyclophosphamide_least_test_gdsc_auc_1se, digits = 2)

## docetaxel
#new_docetaxel_most_sensitive_min <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#docetaxel_most_test_gdsc_auc_min <- auc(docetaxel_test$most_sensitive, new_docetaxel_most_sensitive_min)
#docetaxel_most_test_gdsc_auc_min <- round(docetaxel_most_test_gdsc_auc_min, digits = 2)

new_docetaxel_most_sensitive_1se <- predict(docetaxel_ccle_most_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

docetaxel_most_test_gdsc_auc_1se <- auc(docetaxel_test$most_sensitive, new_docetaxel_most_sensitive_1se)
docetaxel_most_test_gdsc_auc_1se <- round(docetaxel_most_test_gdsc_auc_1se, digits = 2)


new_docetaxel_least_sensitive_min <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

docetaxel_least_test_gdsc_auc_min <- auc(docetaxel_test$least_sensitive, new_docetaxel_least_sensitive_min)
docetaxel_least_test_gdsc_auc_min <- round(docetaxel_least_test_gdsc_auc_min, digits = 2)

#new_docetaxel_least_sensitive_1se <- predict(docetaxel_ccle_least_fit_elnet, newx = as.matrix(docetaxel_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

#docetaxel_least_test_gdsc_auc_1se <- auc(docetaxel_test$least_sensitive, new_docetaxel_least_sensitive_1se)
#docetaxel_least_test_gdsc_auc_1se <- round(docetaxel_least_test_gdsc_auc_1se, digits = 2)

## etoposide
#new_etoposide_most_sensitive_min <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#etoposide_most_test_gdsc_auc_min <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_min)
#etoposide_most_test_gdsc_auc_min <- round(etoposide_most_test_gdsc_auc_min, digits = 2)

new_etoposide_most_sensitive_1se <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_most_test_gdsc_auc_1se <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_1se)
etoposide_most_test_gdsc_auc_1se <- round(etoposide_most_test_gdsc_auc_1se, digits = 2)


#new_etoposide_least_sensitive_min <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#etoposide_least_test_gdsc_auc_min <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_min)
#etoposide_least_test_gdsc_auc_min <- round(etoposide_least_test_gdsc_auc_min, digits = 2)

new_etoposide_least_sensitive_1se <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_least_test_gdsc_auc_1se <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_1se)
etoposide_least_test_gdsc_auc_1se <- round(etoposide_least_test_gdsc_auc_1se, digits = 2)

## fluorouracil
#new_fluorouracil_most_sensitive_min <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#fluorouracil_most_test_gdsc_auc_min <- auc(fluorouracil_test$most_sensitive, new_fluorouracil_most_sensitive_min)
#fluorouracil_most_test_gdsc_auc_min <- round(fluorouracil_most_test_gdsc_auc_min, digits = 2)

new_fluorouracil_most_sensitive_1se <- predict(fluorouracil_ccle_most_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

fluorouracil_most_test_gdsc_auc_1se <- auc(fluorouracil_test$most_sensitive, new_fluorouracil_most_sensitive_1se)
fluorouracil_most_test_gdsc_auc_1se <- round(fluorouracil_most_test_gdsc_auc_1se, digits = 2)


#new_fluorouracil_least_sensitive_min <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#fluorouracil_least_test_gdsc_auc_min <- auc(fluorouracil_test$least_sensitive, new_fluorouracil_least_sensitive_min)
#fluorouracil_least_test_gdsc_auc_min <- round(fluorouracil_least_test_gdsc_auc_min, digits = 2)

new_fluorouracil_least_sensitive_1se <- predict(fluorouracil_ccle_least_fit_elnet, newx = as.matrix(fluorouracil_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

fluorouracil_least_test_gdsc_auc_1se <- auc(fluorouracil_test$least_sensitive, new_fluorouracil_least_sensitive_1se)
fluorouracil_least_test_gdsc_auc_1se <- round(fluorouracil_least_test_gdsc_auc_1se, digits = 2)

## gemcitabine
#new_gemcitabine_most_sensitive_min <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#gemcitabine_most_test_gdsc_auc_min <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_min)
#gemcitabine_most_test_gdsc_auc_min <- round(gemcitabine_most_test_gdsc_auc_min, digits = 2)

new_gemcitabine_most_sensitive_1se <- predict(gemcitabine_ccle_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_most_test_gdsc_auc_1se <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_1se)
gemcitabine_most_test_gdsc_auc_1se <- round(gemcitabine_most_test_gdsc_auc_1se, digits = 2)


#new_gemcitabine_least_sensitive_min <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#gemcitabine_least_test_gdsc_auc_min <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_min)
#gemcitabine_least_test_gdsc_auc_min <- round(gemcitabine_least_test_gdsc_auc_min, digits = 2)

new_gemcitabine_least_sensitive_1se <- predict(gemcitabine_ccle_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_least_test_gdsc_auc_1se <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_1se)
gemcitabine_least_test_gdsc_auc_1se <- round(gemcitabine_least_test_gdsc_auc_1se, digits = 2)

## methotrexate
#new_methotrexate_most_sensitive_min <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#methotrexate_most_test_gdsc_auc_min <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_min)
#methotrexate_most_test_gdsc_auc_min <- round(methotrexate_most_test_gdsc_auc_min, digits = 2)

new_methotrexate_most_sensitive_1se <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_most_test_gdsc_auc_1se <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_1se)
methotrexate_most_test_gdsc_auc_1se <- round(methotrexate_most_test_gdsc_auc_1se, digits = 2)


#new_methotrexate_least_sensitive_min <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

#methotrexate_least_test_gdsc_auc_min <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_min)
#methotrexate_least_test_gdsc_auc_min <- round(methotrexate_least_test_gdsc_auc_min, digits = 2)

new_methotrexate_least_sensitive_1se <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_least_test_gdsc_auc_1se <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_1se)
methotrexate_least_test_gdsc_auc_1se <- round(methotrexate_least_test_gdsc_auc_1se, digits = 2)

overall_auc <- c(carboplatin_most_test_gdsc_auc_1se, carboplatin_least_test_gdsc_auc_1se, 
                 cyclophosphamide_most_test_gdsc_auc_1se, cyclophosphamide_least_test_gdsc_auc_1se, 
                 docetaxel_most_test_gdsc_auc_1se, docetaxel_least_test_gdsc_auc_min, 
                 etoposide_most_test_gdsc_auc_1se, etoposide_least_test_gdsc_auc_1se, 
                 fluorouracil_most_test_gdsc_auc_1se, fluorouracil_least_test_gdsc_auc_1se, 
                 gemcitabine_most_test_gdsc_auc_1se, gemcitabine_least_test_gdsc_auc_1se, 
                 methotrexate_most_test_gdsc_auc_1se, methotrexate_least_test_gdsc_auc_1se)
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

etoposide_autonomic_ganglia_lines              <- etoposide_test$Tissue == 'autonomic_ganglia'
etoposide_biliary_tract_lines                  <- etoposide_test$Tissue == 'biliary_tract'
etoposide_bone_lines                           <- etoposide_test$Tissue == 'bone'
etoposide_breast_lines                         <- etoposide_test$Tissue == 'breast'
etoposide_central_nervous_system_lines         <- etoposide_test$Tissue == 'central_nervous_system'
etoposide_endometrium_lines                    <- etoposide_test$Tissue == 'endometrium'
etoposide_kidney_lines                         <- etoposide_test$Tissue == 'kidney'
etoposide_large_intestine_lines                <- etoposide_test$Tissue == 'large_intestine'
etoposide_liver_lines                          <- etoposide_test$Tissue == 'liver'
etoposide_lung_lines                           <- etoposide_test$Tissue == 'lung'
etoposide_oesaphagus_lines                     <- etoposide_test$Tissue == 'oesphagus'
etoposide_ovary_lines                          <- etoposide_test$Tissue == 'ovary'
etoposide_pancreas_lines                       <- etoposide_test$Tissue == 'pancreas'
etoposide_pleura_lines                         <- etoposide_test$Tissue == 'pleura'
etoposide_prostate_lines                       <- etoposide_test$Tissue == 'prostate'
etoposide_salivary_gland_lines                 <- etoposide_test$Tissue == 'salivary_gland'
etoposide_skin_lines                           <- etoposide_test$Tissue == 'skin'
etoposide_soft_tissue_lines                    <- etoposide_test$Tissue == 'soft_tissue'
etoposide_stomach_lines                        <- etoposide_test$Tissue == 'stomach'
etoposide_thyroid_lines                        <- etoposide_test$Tissue == 'thyroid'
etoposide_upper_aerodigestive_tract_lines      <- etoposide_test$Tissue == 'upper_aerodigestive_tract'
etoposide_urinary_tract_lines                  <- etoposide_test$Tissue == 'urinary_tract'

etoposide_test_scaled_autonomic_ganglia <- etoposide_rna_seq_test_scaled[etoposide_autonomic_ganglia_lines, ]
etoposide_test_autonomic_ganglia <- etoposide_test[which(etoposide_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
etoposide_autonomic_ganglia_most_auc <- auc(etoposide_test_autonomic_ganglia$most_sensitive, new_ic50)

etoposide_test_scaled_biliary_tract <- etoposide_rna_seq_test_scaled[etoposide_biliary_tract_lines, ]
etoposide_test_biliary_tract <- etoposide_test[which(etoposide_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
etoposide_biliary_tract_most_auc <- auc(etoposide_test_biliary_tract$most_sensitive, new_ic50)

etoposide_test_scaled_bone <- etoposide_rna_seq_test_scaled[etoposide_bone_lines, ]
etoposide_test_bone <- etoposide_test[which(etoposide_test$Tissue == 'bone'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
etoposide_bone_most_auc <- auc(etoposide_test_bone$most_sensitive, new_ic50)

etoposide_test_scaled_breast <- etoposide_rna_seq_test_scaled[etoposide_breast_lines, ]
etoposide_test_breast <- etoposide_test[which(etoposide_test$Tissue == 'breast'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
etoposide_breast_most_auc <- auc(etoposide_test_breast$most_sensitive, new_ic50)

etoposide_test_scaled_central_nervous_system <- etoposide_rna_seq_test_scaled[etoposide_central_nervous_system_lines, ]
etoposide_test_central_nervous_system <- etoposide_test[which(etoposide_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
etoposide_central_nervous_system_most_auc <- auc(etoposide_test_central_nervous_system$most_sensitive, new_ic50)

etoposide_test_scaled_endometrium <- etoposide_rna_seq_test_scaled[etoposide_endometrium_lines, ]
etoposide_test_endometrium <- etoposide_test[which(etoposide_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
etoposide_endometrium_most_auc <- auc(etoposide_test_endometrium$most_sensitive, new_ic50)

etoposide_test_scaled_kidney <- etoposide_rna_seq_test_scaled[etoposide_kidney_lines, ]
etoposide_test_kidney <- etoposide_test[which(etoposide_test$Tissue == 'kidney'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
etoposide_kidney_most_auc <- auc(etoposide_test_kidney$most_sensitive, new_ic50)

etoposide_test_scaled_large_intestine <- etoposide_rna_seq_test_scaled[etoposide_large_intestine_lines, ]
etoposide_test_large_intestine <- etoposide_test[which(etoposide_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
etoposide_large_intestine_most_auc <- auc(etoposide_test_large_intestine$most_sensitive, new_ic50)

etoposide_test_scaled_liver <- etoposide_rna_seq_test_scaled[etoposide_liver_lines, ]
etoposide_test_liver <- etoposide_test[which(etoposide_test$Tissue == 'liver'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
etoposide_liver_most_auc <- auc(etoposide_test_liver$most_sensitive, new_ic50)

etoposide_test_scaled_lung <- etoposide_rna_seq_test_scaled[etoposide_lung_lines, ]
etoposide_test_lung <- etoposide_test[which(etoposide_test$Tissue == 'lung'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
etoposide_lung_most_auc <- auc(etoposide_test_lung$most_sensitive, new_ic50)

etoposide_test_scaled_oesaphagus <- etoposide_rna_seq_test_scaled[etoposide_oesaphagus_lines, ]
etoposide_test_oesaphagus <- etoposide_test[which(etoposide_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
etoposide_oesaphagus_most_auc <- auc(etoposide_test_oesaphagus$most_sensitive, new_ic50)

etoposide_test_scaled_ovary <- etoposide_rna_seq_test_scaled[etoposide_ovary_lines, ]
etoposide_test_ovary <- etoposide_test[which(etoposide_test$Tissue == 'ovary'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
etoposide_ovary_most_auc <- auc(etoposide_test_ovary$most_sensitive, new_ic50)

etoposide_test_scaled_pancreas <- etoposide_rna_seq_test_scaled[etoposide_pancreas_lines, ]
etoposide_test_pancreas <- etoposide_test[which(etoposide_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
etoposide_pancreas_most_auc <- auc(etoposide_test_pancreas$most_sensitive, new_ic50)

etoposide_test_scaled_pleura <- etoposide_rna_seq_test_scaled[etoposide_pleura_lines, ]
etoposide_test_pleura <- etoposide_test[which(etoposide_test$Tissue == 'pleura'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
etoposide_pleura_most_auc <- auc(etoposide_test_pleura$most_sensitive, new_ic50)

etoposide_test_scaled_prostate <- etoposide_rna_seq_test_scaled[etoposide_prostate_lines, ]
etoposide_test_prostate <- etoposide_test[which(etoposide_test$Tissue == 'prostate'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
etoposide_prostate_most_auc <- auc(etoposide_test_prostate$most_sensitive, new_ic50)

etoposide_test_scaled_salivary_gland <- etoposide_rna_seq_test_scaled[etoposide_salivary_gland_lines, ]
etoposide_test_salivary_gland <- etoposide_test[which(etoposide_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
etoposide_salivary_gland_most_auc <- auc(etoposide_test_salivary_gland$most_sensitive, new_ic50)

etoposide_test_scaled_skin <- etoposide_rna_seq_test_scaled[etoposide_skin_lines, ]
etoposide_test_skin <- etoposide_test[which(etoposide_test$Tissue == 'skin'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
etoposide_skin_most_auc <- auc(etoposide_test_skin$most_sensitive, new_ic50)

etoposide_test_scaled_soft_tissue <- etoposide_rna_seq_test_scaled[etoposide_soft_tissue_lines, ]
etoposide_test_soft_tissue <- etoposide_test[which(etoposide_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
etoposide_soft_tissue_most_auc <- auc(etoposide_test_soft_tissue$most_sensitive, new_ic50)

etoposide_test_scaled_stomach <- etoposide_rna_seq_test_scaled[etoposide_stomach_lines, ]
etoposide_test_stomach <- etoposide_test[which(etoposide_test$Tissue == 'stomach'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
etoposide_stomach_most_auc <- auc(etoposide_test_stomach$most_sensitive, new_ic50)

etoposide_test_scaled_thyroid <- etoposide_rna_seq_test_scaled[etoposide_thyroid_lines, ]
etoposide_test_thyroid <- etoposide_test[which(etoposide_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
etoposide_thyroid_most_auc <- auc(etoposide_test_thyroid$most_sensitive, new_ic50)

etoposide_test_scaled_upper_aerodigestive_tract <- etoposide_rna_seq_test_scaled[etoposide_upper_aerodigestive_tract_lines, ]
etoposide_test_upper_aerodigestive_tract <- etoposide_test[which(etoposide_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
etoposide_upper_aerodigestive_tract_most_auc <- auc(etoposide_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

etoposide_test_scaled_urinary_tract <- etoposide_rna_seq_test_scaled[etoposide_urinary_tract_lines, ]
etoposide_test_urinary_tract <- etoposide_test[which(etoposide_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(etoposide_ccle_most_fit_elnet, newx = as.matrix(etoposide_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
etoposide_urinary_tract_most_auc <- auc(etoposide_test_urinary_tract$most_sensitive, new_ic50)

etoposide_most_auc <- c(etoposide_autonomic_ganglia_most_auc, etoposide_biliary_tract_most_auc, 
                          etoposide_bone_most_auc, etoposide_breast_most_auc, etoposide_central_nervous_system_most_auc, 
                          etoposide_endometrium_most_auc, etoposide_kidney_most_auc, 
                          etoposide_large_intestine_most_auc, etoposide_liver_most_auc, 
                          etoposide_lung_most_auc, etoposide_oesaphagus_most_auc, 
                          etoposide_ovary_most_auc, etoposide_pancreas_most_auc, 
                          etoposide_pleura_most_auc, etoposide_prostate_most_auc, 
                          etoposide_salivary_gland_most_auc, etoposide_skin_most_auc, 
                          etoposide_soft_tissue_most_auc, etoposide_stomach_most_auc, 
                          etoposide_thyroid_most_auc, etoposide_upper_aerodigestive_tract_most_auc, 
                          etoposide_urinary_tract_most_auc)

etoposide_test_scaled_autonomic_ganglia <- etoposide_rna_seq_test_scaled[etoposide_autonomic_ganglia_lines, ]
etoposide_test_autonomic_ganglia <- etoposide_test[which(etoposide_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
etoposide_autonomic_ganglia_least_auc <- auc(etoposide_test_autonomic_ganglia$least_sensitive, new_ic50)

etoposide_test_scaled_biliary_tract <- etoposide_rna_seq_test_scaled[etoposide_biliary_tract_lines, ]
etoposide_test_biliary_tract <- etoposide_test[which(etoposide_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
etoposide_biliary_tract_least_auc <- auc(etoposide_test_biliary_tract$least_sensitive, new_ic50)

etoposide_test_scaled_bone <- etoposide_rna_seq_test_scaled[etoposide_bone_lines, ]
etoposide_test_bone <- etoposide_test[which(etoposide_test$Tissue == 'bone'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
etoposide_bone_least_auc <- auc(etoposide_test_bone$least_sensitive, new_ic50)

etoposide_test_scaled_breast <- etoposide_rna_seq_test_scaled[etoposide_breast_lines, ]
etoposide_test_breast <- etoposide_test[which(etoposide_test$Tissue == 'breast'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
etoposide_breast_least_auc <- auc(etoposide_test_breast$least_sensitive, new_ic50)

etoposide_test_scaled_central_nervous_system <- etoposide_rna_seq_test_scaled[etoposide_central_nervous_system_lines, ]
etoposide_test_central_nervous_system <- etoposide_test[which(etoposide_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
etoposide_central_nervous_system_least_auc <- auc(etoposide_test_central_nervous_system$least_sensitive, new_ic50)

etoposide_test_scaled_endometrium <- etoposide_rna_seq_test_scaled[etoposide_endometrium_lines, ]
etoposide_test_endometrium <- etoposide_test[which(etoposide_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
etoposide_endometrium_least_auc <- auc(etoposide_test_endometrium$least_sensitive, new_ic50)

etoposide_test_scaled_kidney <- etoposide_rna_seq_test_scaled[etoposide_kidney_lines, ]
etoposide_test_kidney <- etoposide_test[which(etoposide_test$Tissue == 'kidney'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
etoposide_kidney_least_auc <- auc(etoposide_test_kidney$least_sensitive, new_ic50)

etoposide_test_scaled_large_intestine <- etoposide_rna_seq_test_scaled[etoposide_large_intestine_lines, ]
etoposide_test_large_intestine <- etoposide_test[which(etoposide_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
etoposide_large_intestine_least_auc <- auc(etoposide_test_large_intestine$least_sensitive, new_ic50)

etoposide_test_scaled_liver <- etoposide_rna_seq_test_scaled[etoposide_liver_lines, ]
etoposide_test_liver <- etoposide_test[which(etoposide_test$Tissue == 'liver'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
etoposide_liver_least_auc <- auc(etoposide_test_liver$least_sensitive, new_ic50)

etoposide_test_scaled_lung <- etoposide_rna_seq_test_scaled[etoposide_lung_lines, ]
etoposide_test_lung <- etoposide_test[which(etoposide_test$Tissue == 'lung'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
etoposide_lung_least_auc <- auc(etoposide_test_lung$least_sensitive, new_ic50)

etoposide_test_scaled_oesaphagus <- etoposide_rna_seq_test_scaled[etoposide_oesaphagus_lines, ]
etoposide_test_oesaphagus <- etoposide_test[which(etoposide_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
etoposide_oesaphagus_least_auc <- auc(etoposide_test_oesaphagus$least_sensitive, new_ic50)

etoposide_test_scaled_ovary <- etoposide_rna_seq_test_scaled[etoposide_ovary_lines, ]
etoposide_test_ovary <- etoposide_test[which(etoposide_test$Tissue == 'ovary'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
etoposide_ovary_least_auc <- auc(etoposide_test_ovary$least_sensitive, new_ic50)

etoposide_test_scaled_pancreas <- etoposide_rna_seq_test_scaled[etoposide_pancreas_lines, ]
etoposide_test_pancreas <- etoposide_test[which(etoposide_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
etoposide_pancreas_least_auc <- auc(etoposide_test_pancreas$least_sensitive, new_ic50)

etoposide_test_scaled_pleura <- etoposide_rna_seq_test_scaled[etoposide_pleura_lines, ]
etoposide_test_pleura <- etoposide_test[which(etoposide_test$Tissue == 'pleura'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
etoposide_pleura_least_auc <- auc(etoposide_test_pleura$least_sensitive, new_ic50)

etoposide_test_scaled_prostate <- etoposide_rna_seq_test_scaled[etoposide_prostate_lines, ]
etoposide_test_prostate <- etoposide_test[which(etoposide_test$Tissue == 'prostate'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
etoposide_prostate_least_auc <- auc(etoposide_test_prostate$least_sensitive, new_ic50)

etoposide_test_scaled_salivary_gland <- etoposide_rna_seq_test_scaled[etoposide_salivary_gland_lines, ]
etoposide_test_salivary_gland <- etoposide_test[which(etoposide_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
etoposide_salivary_gland_least_auc <- auc(etoposide_test_salivary_gland$least_sensitive, new_ic50)

etoposide_test_scaled_skin <- etoposide_rna_seq_test_scaled[etoposide_skin_lines, ]
etoposide_test_skin <- etoposide_test[which(etoposide_test$Tissue == 'skin'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
etoposide_skin_least_auc <- auc(etoposide_test_skin$least_sensitive, new_ic50)

etoposide_test_scaled_soft_tissue <- etoposide_rna_seq_test_scaled[etoposide_soft_tissue_lines, ]
etoposide_test_soft_tissue <- etoposide_test[which(etoposide_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
etoposide_soft_tissue_least_auc <- auc(etoposide_test_soft_tissue$least_sensitive, new_ic50)

etoposide_test_scaled_stomach <- etoposide_rna_seq_test_scaled[etoposide_stomach_lines, ]
etoposide_test_stomach <- etoposide_test[which(etoposide_test$Tissue == 'stomach'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
etoposide_stomach_least_auc <- auc(etoposide_test_stomach$least_sensitive, new_ic50)

etoposide_test_scaled_thyroid <- etoposide_rna_seq_test_scaled[etoposide_thyroid_lines, ]
etoposide_test_thyroid <- etoposide_test[which(etoposide_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
etoposide_thyroid_least_auc <- auc(etoposide_test_thyroid$least_sensitive, new_ic50)

etoposide_test_scaled_upper_aerodigestive_tract <- etoposide_rna_seq_test_scaled[etoposide_upper_aerodigestive_tract_lines, ]
etoposide_test_upper_aerodigestive_tract <- etoposide_test[which(etoposide_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
etoposide_upper_aerodigestive_tract_least_auc <- auc(etoposide_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

etoposide_test_scaled_urinary_tract <- etoposide_rna_seq_test_scaled[etoposide_urinary_tract_lines, ]
etoposide_test_urinary_tract <- etoposide_test[which(etoposide_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(etoposide_ccle_least_fit_elnet, newx = as.matrix(etoposide_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
etoposide_urinary_tract_least_auc <- auc(etoposide_test_urinary_tract$least_sensitive, new_ic50)

etoposide_least_auc <- c(etoposide_autonomic_ganglia_least_auc, etoposide_biliary_tract_least_auc, 
                           etoposide_bone_least_auc, etoposide_breast_least_auc, etoposide_central_nervous_system_least_auc, 
                           etoposide_endometrium_least_auc, etoposide_kidney_least_auc, 
                           etoposide_large_intestine_least_auc, etoposide_liver_least_auc, 
                           etoposide_lung_least_auc, etoposide_oesaphagus_least_auc, 
                           etoposide_ovary_least_auc, etoposide_pancreas_least_auc, 
                           etoposide_pleura_least_auc, etoposide_prostate_least_auc, 
                           etoposide_salivary_gland_least_auc, etoposide_skin_least_auc, 
                           etoposide_soft_tissue_least_auc, etoposide_stomach_least_auc, 
                           etoposide_thyroid_least_auc, etoposide_upper_aerodigestive_tract_least_auc, 
                           etoposide_urinary_tract_least_auc)

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

methotrexate_autonomic_ganglia_lines              <- methotrexate_test$Tissue == 'autonomic_ganglia'
methotrexate_biliary_tract_lines                  <- methotrexate_test$Tissue == 'biliary_tract'
methotrexate_bone_lines                           <- methotrexate_test$Tissue == 'bone'
methotrexate_breast_lines                         <- methotrexate_test$Tissue == 'breast'
methotrexate_central_nervous_system_lines         <- methotrexate_test$Tissue == 'central_nervous_system'
methotrexate_endometrium_lines                    <- methotrexate_test$Tissue == 'endometrium'
methotrexate_kidney_lines                         <- methotrexate_test$Tissue == 'kidney'
methotrexate_large_intestine_lines                <- methotrexate_test$Tissue == 'large_intestine'
methotrexate_liver_lines                          <- methotrexate_test$Tissue == 'liver'
methotrexate_lung_lines                           <- methotrexate_test$Tissue == 'lung'
methotrexate_oesaphagus_lines                     <- methotrexate_test$Tissue == 'oesphagus'
methotrexate_ovary_lines                          <- methotrexate_test$Tissue == 'ovary'
methotrexate_pancreas_lines                       <- methotrexate_test$Tissue == 'pancreas'
methotrexate_pleura_lines                         <- methotrexate_test$Tissue == 'pleura'
methotrexate_prostate_lines                       <- methotrexate_test$Tissue == 'prostate'
methotrexate_salivary_gland_lines                 <- methotrexate_test$Tissue == 'salivary_gland'
methotrexate_skin_lines                           <- methotrexate_test$Tissue == 'skin'
methotrexate_soft_tissue_lines                    <- methotrexate_test$Tissue == 'soft_tissue'
methotrexate_stomach_lines                        <- methotrexate_test$Tissue == 'stomach'
methotrexate_thyroid_lines                        <- methotrexate_test$Tissue == 'thyroid'
methotrexate_upper_aerodigestive_tract_lines      <- methotrexate_test$Tissue == 'upper_aerodigestive_tract'
methotrexate_urinary_tract_lines                  <- methotrexate_test$Tissue == 'urinary_tract'

methotrexate_test_scaled_autonomic_ganglia <- methotrexate_rna_seq_test_scaled[methotrexate_autonomic_ganglia_lines, ]
methotrexate_test_autonomic_ganglia <- methotrexate_test[which(methotrexate_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
methotrexate_autonomic_ganglia_most_auc <- auc(methotrexate_test_autonomic_ganglia$most_sensitive, new_ic50)

methotrexate_test_scaled_biliary_tract <- methotrexate_rna_seq_test_scaled[methotrexate_biliary_tract_lines, ]
methotrexate_test_biliary_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_biliary_tract_most_auc <- auc(methotrexate_test_biliary_tract$most_sensitive, new_ic50)

methotrexate_test_scaled_bone <- methotrexate_rna_seq_test_scaled[methotrexate_bone_lines, ]
methotrexate_test_bone <- methotrexate_test[which(methotrexate_test$Tissue == 'bone'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
methotrexate_bone_most_auc <- auc(methotrexate_test_bone$most_sensitive, new_ic50)

methotrexate_test_scaled_breast <- methotrexate_rna_seq_test_scaled[methotrexate_breast_lines, ]
methotrexate_test_breast <- methotrexate_test[which(methotrexate_test$Tissue == 'breast'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
methotrexate_breast_most_auc <- auc(methotrexate_test_breast$most_sensitive, new_ic50)

methotrexate_test_scaled_central_nervous_system <- methotrexate_rna_seq_test_scaled[methotrexate_central_nervous_system_lines, ]
methotrexate_test_central_nervous_system <- methotrexate_test[which(methotrexate_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
methotrexate_central_nervous_system_most_auc <- auc(methotrexate_test_central_nervous_system$most_sensitive, new_ic50)

methotrexate_test_scaled_endometrium <- methotrexate_rna_seq_test_scaled[methotrexate_endometrium_lines, ]
methotrexate_test_endometrium <- methotrexate_test[which(methotrexate_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
methotrexate_endometrium_most_auc <- auc(methotrexate_test_endometrium$most_sensitive, new_ic50)

methotrexate_test_scaled_kidney <- methotrexate_rna_seq_test_scaled[methotrexate_kidney_lines, ]
methotrexate_test_kidney <- methotrexate_test[which(methotrexate_test$Tissue == 'kidney'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
methotrexate_kidney_most_auc <- auc(methotrexate_test_kidney$most_sensitive, new_ic50)

methotrexate_test_scaled_large_intestine <- methotrexate_rna_seq_test_scaled[methotrexate_large_intestine_lines, ]
methotrexate_test_large_intestine <- methotrexate_test[which(methotrexate_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
methotrexate_large_intestine_most_auc <- auc(methotrexate_test_large_intestine$most_sensitive, new_ic50)

methotrexate_test_scaled_liver <- methotrexate_rna_seq_test_scaled[methotrexate_liver_lines, ]
methotrexate_test_liver <- methotrexate_test[which(methotrexate_test$Tissue == 'liver'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
methotrexate_liver_most_auc <- auc(methotrexate_test_liver$most_sensitive, new_ic50)

methotrexate_test_scaled_lung <- methotrexate_rna_seq_test_scaled[methotrexate_lung_lines, ]
methotrexate_test_lung <- methotrexate_test[which(methotrexate_test$Tissue == 'lung'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
methotrexate_lung_most_auc <- auc(methotrexate_test_lung$most_sensitive, new_ic50)

methotrexate_test_scaled_oesaphagus <- methotrexate_rna_seq_test_scaled[methotrexate_oesaphagus_lines, ]
methotrexate_test_oesaphagus <- methotrexate_test[which(methotrexate_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
methotrexate_oesaphagus_most_auc <- auc(methotrexate_test_oesaphagus$most_sensitive, new_ic50)

methotrexate_test_scaled_ovary <- methotrexate_rna_seq_test_scaled[methotrexate_ovary_lines, ]
methotrexate_test_ovary <- methotrexate_test[which(methotrexate_test$Tissue == 'ovary'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
methotrexate_ovary_most_auc <- auc(methotrexate_test_ovary$most_sensitive, new_ic50)

methotrexate_test_scaled_pancreas <- methotrexate_rna_seq_test_scaled[methotrexate_pancreas_lines, ]
methotrexate_test_pancreas <- methotrexate_test[which(methotrexate_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
methotrexate_pancreas_most_auc <- auc(methotrexate_test_pancreas$most_sensitive, new_ic50)

methotrexate_test_scaled_pleura <- methotrexate_rna_seq_test_scaled[methotrexate_pleura_lines, ]
methotrexate_test_pleura <- methotrexate_test[which(methotrexate_test$Tissue == 'pleura'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
methotrexate_pleura_most_auc <- auc(methotrexate_test_pleura$most_sensitive, new_ic50)

methotrexate_test_scaled_prostate <- methotrexate_rna_seq_test_scaled[methotrexate_prostate_lines, ]
methotrexate_test_prostate <- methotrexate_test[which(methotrexate_test$Tissue == 'prostate'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
methotrexate_prostate_most_auc <- auc(methotrexate_test_prostate$most_sensitive, new_ic50)

methotrexate_test_scaled_salivary_gland <- methotrexate_rna_seq_test_scaled[methotrexate_salivary_gland_lines, ]
methotrexate_test_salivary_gland <- methotrexate_test[which(methotrexate_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
methotrexate_salivary_gland_most_auc <- auc(methotrexate_test_salivary_gland$most_sensitive, new_ic50)

methotrexate_test_scaled_skin <- methotrexate_rna_seq_test_scaled[methotrexate_skin_lines, ]
methotrexate_test_skin <- methotrexate_test[which(methotrexate_test$Tissue == 'skin'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
methotrexate_skin_most_auc <- auc(methotrexate_test_skin$most_sensitive, new_ic50)

methotrexate_test_scaled_soft_tissue <- methotrexate_rna_seq_test_scaled[methotrexate_soft_tissue_lines, ]
methotrexate_test_soft_tissue <- methotrexate_test[which(methotrexate_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
methotrexate_soft_tissue_most_auc <- auc(methotrexate_test_soft_tissue$most_sensitive, new_ic50)

methotrexate_test_scaled_stomach <- methotrexate_rna_seq_test_scaled[methotrexate_stomach_lines, ]
methotrexate_test_stomach <- methotrexate_test[which(methotrexate_test$Tissue == 'stomach'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
methotrexate_stomach_most_auc <- auc(methotrexate_test_stomach$most_sensitive, new_ic50)

methotrexate_test_scaled_thyroid <- methotrexate_rna_seq_test_scaled[methotrexate_thyroid_lines, ]
methotrexate_test_thyroid <- methotrexate_test[which(methotrexate_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
methotrexate_thyroid_most_auc <- auc(methotrexate_test_thyroid$most_sensitive, new_ic50)

methotrexate_test_scaled_upper_aerodigestive_tract <- methotrexate_rna_seq_test_scaled[methotrexate_upper_aerodigestive_tract_lines, ]
methotrexate_test_upper_aerodigestive_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_upper_aerodigestive_tract_most_auc <- auc(methotrexate_test_upper_aerodigestive_tract$most_sensitive, new_ic50)

methotrexate_test_scaled_urinary_tract <- methotrexate_rna_seq_test_scaled[methotrexate_urinary_tract_lines, ]
methotrexate_test_urinary_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(methotrexate_ccle_most_fit_elnet, newx = as.matrix(methotrexate_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_urinary_tract_most_auc <- auc(methotrexate_test_urinary_tract$most_sensitive, new_ic50)

methotrexate_most_auc <- c(methotrexate_autonomic_ganglia_most_auc, methotrexate_biliary_tract_most_auc, 
                          methotrexate_bone_most_auc, methotrexate_breast_most_auc, methotrexate_central_nervous_system_most_auc, 
                          methotrexate_endometrium_most_auc, methotrexate_kidney_most_auc, 
                          methotrexate_large_intestine_most_auc, methotrexate_liver_most_auc, 
                          methotrexate_lung_most_auc, methotrexate_oesaphagus_most_auc, 
                          methotrexate_ovary_most_auc, methotrexate_pancreas_most_auc, 
                          methotrexate_pleura_most_auc, methotrexate_prostate_most_auc, 
                          methotrexate_salivary_gland_most_auc, methotrexate_skin_most_auc, 
                          methotrexate_soft_tissue_most_auc, methotrexate_stomach_most_auc, 
                          methotrexate_thyroid_most_auc, methotrexate_upper_aerodigestive_tract_most_auc, 
                          methotrexate_urinary_tract_most_auc)

methotrexate_test_scaled_autonomic_ganglia <- methotrexate_rna_seq_test_scaled[methotrexate_autonomic_ganglia_lines, ]
methotrexate_test_autonomic_ganglia <- methotrexate_test[which(methotrexate_test$Tissue == 'autonomic_ganglia'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_autonomic_ganglia), s = 'lambda.1se', interval = 'conf')
methotrexate_autonomic_ganglia_least_auc <- auc(methotrexate_test_autonomic_ganglia$least_sensitive, new_ic50)

methotrexate_test_scaled_biliary_tract <- methotrexate_rna_seq_test_scaled[methotrexate_biliary_tract_lines, ]
methotrexate_test_biliary_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'biliary_tract'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_biliary_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_biliary_tract_least_auc <- auc(methotrexate_test_biliary_tract$least_sensitive, new_ic50)

methotrexate_test_scaled_bone <- methotrexate_rna_seq_test_scaled[methotrexate_bone_lines, ]
methotrexate_test_bone <- methotrexate_test[which(methotrexate_test$Tissue == 'bone'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_bone), s = 'lambda.1se', interval = 'conf')
methotrexate_bone_least_auc <- auc(methotrexate_test_bone$least_sensitive, new_ic50)

methotrexate_test_scaled_breast <- methotrexate_rna_seq_test_scaled[methotrexate_breast_lines, ]
methotrexate_test_breast <- methotrexate_test[which(methotrexate_test$Tissue == 'breast'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_breast), s = 'lambda.1se', interval = 'conf')
methotrexate_breast_least_auc <- auc(methotrexate_test_breast$least_sensitive, new_ic50)

methotrexate_test_scaled_central_nervous_system <- methotrexate_rna_seq_test_scaled[methotrexate_central_nervous_system_lines, ]
methotrexate_test_central_nervous_system <- methotrexate_test[which(methotrexate_test$Tissue == 'central_nervous_system'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_central_nervous_system), s = 'lambda.1se', interval = 'conf')
methotrexate_central_nervous_system_least_auc <- auc(methotrexate_test_central_nervous_system$least_sensitive, new_ic50)

methotrexate_test_scaled_endometrium <- methotrexate_rna_seq_test_scaled[methotrexate_endometrium_lines, ]
methotrexate_test_endometrium <- methotrexate_test[which(methotrexate_test$Tissue == 'endometrium'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_endometrium), s = 'lambda.1se', interval = 'conf')
methotrexate_endometrium_least_auc <- auc(methotrexate_test_endometrium$least_sensitive, new_ic50)

methotrexate_test_scaled_kidney <- methotrexate_rna_seq_test_scaled[methotrexate_kidney_lines, ]
methotrexate_test_kidney <- methotrexate_test[which(methotrexate_test$Tissue == 'kidney'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_kidney), s = 'lambda.1se', interval = 'conf')
methotrexate_kidney_least_auc <- auc(methotrexate_test_kidney$least_sensitive, new_ic50)

methotrexate_test_scaled_large_intestine <- methotrexate_rna_seq_test_scaled[methotrexate_large_intestine_lines, ]
methotrexate_test_large_intestine <- methotrexate_test[which(methotrexate_test$Tissue == 'large_intestine'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_large_intestine), s = 'lambda.1se', interval = 'conf')
methotrexate_large_intestine_least_auc <- auc(methotrexate_test_large_intestine$least_sensitive, new_ic50)

methotrexate_test_scaled_liver <- methotrexate_rna_seq_test_scaled[methotrexate_liver_lines, ]
methotrexate_test_liver <- methotrexate_test[which(methotrexate_test$Tissue == 'liver'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_liver), s = 'lambda.1se', interval = 'conf')
methotrexate_liver_least_auc <- auc(methotrexate_test_liver$least_sensitive, new_ic50)

methotrexate_test_scaled_lung <- methotrexate_rna_seq_test_scaled[methotrexate_lung_lines, ]
methotrexate_test_lung <- methotrexate_test[which(methotrexate_test$Tissue == 'lung'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_lung), s = 'lambda.1se', interval = 'conf')
methotrexate_lung_least_auc <- auc(methotrexate_test_lung$least_sensitive, new_ic50)

methotrexate_test_scaled_oesaphagus <- methotrexate_rna_seq_test_scaled[methotrexate_oesaphagus_lines, ]
methotrexate_test_oesaphagus <- methotrexate_test[which(methotrexate_test$Tissue == 'oesaphagus'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_oesaphagus), s = 'lambda.1se', interval = 'conf')
methotrexate_oesaphagus_least_auc <- auc(methotrexate_test_oesaphagus$least_sensitive, new_ic50)

methotrexate_test_scaled_ovary <- methotrexate_rna_seq_test_scaled[methotrexate_ovary_lines, ]
methotrexate_test_ovary <- methotrexate_test[which(methotrexate_test$Tissue == 'ovary'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_ovary), s = 'lambda.1se', interval = 'conf')
methotrexate_ovary_least_auc <- auc(methotrexate_test_ovary$least_sensitive, new_ic50)

methotrexate_test_scaled_pancreas <- methotrexate_rna_seq_test_scaled[methotrexate_pancreas_lines, ]
methotrexate_test_pancreas <- methotrexate_test[which(methotrexate_test$Tissue == 'pancreas'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pancreas), s = 'lambda.1se', interval = 'conf')
methotrexate_pancreas_least_auc <- auc(methotrexate_test_pancreas$least_sensitive, new_ic50)

methotrexate_test_scaled_pleura <- methotrexate_rna_seq_test_scaled[methotrexate_pleura_lines, ]
methotrexate_test_pleura <- methotrexate_test[which(methotrexate_test$Tissue == 'pleura'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_pleura), s = 'lambda.1se', interval = 'conf')
methotrexate_pleura_least_auc <- auc(methotrexate_test_pleura$least_sensitive, new_ic50)

methotrexate_test_scaled_prostate <- methotrexate_rna_seq_test_scaled[methotrexate_prostate_lines, ]
methotrexate_test_prostate <- methotrexate_test[which(methotrexate_test$Tissue == 'prostate'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_prostate), s = 'lambda.1se', interval = 'conf')
methotrexate_prostate_least_auc <- auc(methotrexate_test_prostate$least_sensitive, new_ic50)

methotrexate_test_scaled_salivary_gland <- methotrexate_rna_seq_test_scaled[methotrexate_salivary_gland_lines, ]
methotrexate_test_salivary_gland <- methotrexate_test[which(methotrexate_test$Tissue == 'salivary_gland'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_salivary_gland), s = 'lambda.1se', interval = 'conf')
methotrexate_salivary_gland_least_auc <- auc(methotrexate_test_salivary_gland$least_sensitive, new_ic50)

methotrexate_test_scaled_skin <- methotrexate_rna_seq_test_scaled[methotrexate_skin_lines, ]
methotrexate_test_skin <- methotrexate_test[which(methotrexate_test$Tissue == 'skin'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_skin), s = 'lambda.1se', interval = 'conf')
methotrexate_skin_least_auc <- auc(methotrexate_test_skin$least_sensitive, new_ic50)

methotrexate_test_scaled_soft_tissue <- methotrexate_rna_seq_test_scaled[methotrexate_soft_tissue_lines, ]
methotrexate_test_soft_tissue <- methotrexate_test[which(methotrexate_test$Tissue == 'soft_tissue'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_soft_tissue), s = 'lambda.1se', interval = 'conf')
methotrexate_soft_tissue_least_auc <- auc(methotrexate_test_soft_tissue$least_sensitive, new_ic50)

methotrexate_test_scaled_stomach <- methotrexate_rna_seq_test_scaled[methotrexate_stomach_lines, ]
methotrexate_test_stomach <- methotrexate_test[which(methotrexate_test$Tissue == 'stomach'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_stomach), s = 'lambda.1se', interval = 'conf')
methotrexate_stomach_least_auc <- auc(methotrexate_test_stomach$least_sensitive, new_ic50)

methotrexate_test_scaled_thyroid <- methotrexate_rna_seq_test_scaled[methotrexate_thyroid_lines, ]
methotrexate_test_thyroid <- methotrexate_test[which(methotrexate_test$Tissue == 'thyroid'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_thyroid), s = 'lambda.1se', interval = 'conf')
methotrexate_thyroid_least_auc <- auc(methotrexate_test_thyroid$least_sensitive, new_ic50)

methotrexate_test_scaled_upper_aerodigestive_tract <- methotrexate_rna_seq_test_scaled[methotrexate_upper_aerodigestive_tract_lines, ]
methotrexate_test_upper_aerodigestive_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'upper_aerodigestive_tract'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_upper_aerodigestive_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_upper_aerodigestive_tract_least_auc <- auc(methotrexate_test_upper_aerodigestive_tract$least_sensitive, new_ic50)

methotrexate_test_scaled_urinary_tract <- methotrexate_rna_seq_test_scaled[methotrexate_urinary_tract_lines, ]
methotrexate_test_urinary_tract <- methotrexate_test[which(methotrexate_test$Tissue == 'urinary_tract'), ]
new_ic50 <- predict(methotrexate_ccle_least_fit_elnet, newx = as.matrix(methotrexate_test_scaled_urinary_tract), s = 'lambda.1se', interval = 'conf')
methotrexate_urinary_tract_least_auc <- auc(methotrexate_test_urinary_tract$least_sensitive, new_ic50)

methotrexate_least_auc <- c(methotrexate_autonomic_ganglia_least_auc, methotrexate_biliary_tract_least_auc, 
                           methotrexate_bone_least_auc, methotrexate_breast_least_auc, methotrexate_central_nervous_system_least_auc, 
                           methotrexate_endometrium_least_auc, methotrexate_kidney_least_auc, 
                           methotrexate_large_intestine_least_auc, methotrexate_liver_least_auc, 
                           methotrexate_lung_least_auc, methotrexate_oesaphagus_least_auc, 
                           methotrexate_ovary_least_auc, methotrexate_pancreas_least_auc, 
                           methotrexate_pleura_least_auc, methotrexate_prostate_least_auc, 
                           methotrexate_salivary_gland_least_auc, methotrexate_skin_least_auc, 
                           methotrexate_soft_tissue_least_auc, methotrexate_stomach_least_auc, 
                           methotrexate_thyroid_least_auc, methotrexate_upper_aerodigestive_tract_least_auc, 
                           methotrexate_urinary_tract_least_auc)

all_ccle_auc_by_type <- data.frame(carboplatin_most_auc, carboplatin_least_auc, 
                                   cyclophosphamide_most_auc, cyclophosphamide_least_auc, 
                                   docetaxel_most_auc, docetaxel_least_auc, 
                                   etoposide_most_auc, etoposide_least_auc, 
                                   fluorouracil_most_auc, fluorouracil_least_auc, 
                                   gemcitabine_most_auc, gemcitabine_least_auc, 
                                   methotrexate_most_auc, methotrexate_least_auc)

#all_ccle_auc_by_type <- data.frame(t(all_ccle_auc_by_type))
all_ccle_auc_by_type <- rbind(all_ccle_auc_by_type, overall_auc)
#all_ccle_auc_by_type <- all_ccle_auc_by_type[c(15, 1:14), ]
colnames(all_ccle_auc_by_type) <- c('carboplatin_most', 'carboplatin_least', 
                                    'cyclophosphamide_most', 'cyclophosphamide_least', 'docetaxel_most', 
                                    'docetaxel_least', 'etoposide_most', 'etoposide_least', 
                                    'fluorouracil_most', 'fluorouracil_least', 
                                    'gemcitabine_most', 'gemcitabine_least', 'methotrexate_most', 
                                    'methotrexate_least')
rownames(all_ccle_auc_by_type) <- c('autonomic_ganglia', 'biliary_tract', 'bone', 'breast', 'central_nervous_system', 'endometrium', 'kidney', 
                                    'large_intestine', 'liver', 'lung', 'oesaphagus', 'ovary', 
                                    'pancreas', 'pleura', 'prostate', 'salivary_gland', 'skin', 
                                    'soft_tissue', 'stomach', 'thyroid', 'upper_aerodigestive_tract', 'urinary_tract', 'overall')
for (i in 1:nrow(all_ccle_auc_by_type)) {
  for (j in 1:ncol(all_ccle_auc_by_type)) {
    if (all_ccle_auc_by_type[i,j] == 'NaN') {
      all_ccle_auc_by_type[i,j] <- 0.50
    }
  }
}
for (i in 1:nrow(all_ccle_auc_by_type)) {
  for (j in 1:ncol(all_ccle_auc_by_type)) {
    if (all_ccle_auc_by_type[i,j] < 0.50) {
      all_ccle_auc_by_type[i,j] <- 0.50
    }
  }
}
all_ccle_auc_by_type <- round(all_ccle_auc_by_type, digits = 2)
all_ccle_auc_by_type <- data.frame(t(all_ccle_auc_by_type))
all_ccle_auc_by_type <- all_ccle_auc_by_type[, c(23, 1:22)]
colors <- colorRampPalette(c("dodgerblue", "white", "red"))(100)

png(filename = 'Images/CCLE_AUC_heatmap.png')
heatmap.2(as.matrix(all_ccle_auc_by_type), trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, density.info = 'none', key.xlab = 'AUC', key.title = '', cexRow = 0.9, cexCol = 0.8, cellnote = all_ccle_auc_by_type, notecol = 'black', colsep = 1, sepwidth = c(0.1,0.1), srtCol = 45, margins = c(8,8))
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
plot(perf2, col = 'red', lty = 1, main = 'BLCA treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.66)', 'least_sensitive (AUC = 0.62)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
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
plot(perf2, col = 'red', lty = 1, main = 'BRCA treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 1.0)', 'least_sensitive (AUC = 1.0)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'ESCA treated with fluorouracil')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.58)', 'least_sensitive (AUC = 0.92)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'LUSC treated with carboplatin')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.85)', 'least_sensitive (AUC = 0.50)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
plot(perf, col = 'red')
plot(perf2, col = 'red', lty = 2)
plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 2, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_min (0.83)', 'most_1se (0.83)', 'least_min (0.50)', 'least_1se (0.40)'), lty = c(1,2,1,2), col = c('red', 'red', 'blue', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'SARC treated with gemcitabine (CCLE)')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.64)', 'least_sensitive (AUC = 1.0)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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
plot(perf2, col = 'red', lty = 1, main = 'UCS treated with carboplatin')
#plot(perf3, col = 'blue', add = TRUE)
plot(perf4, col = 'blue', lty = 1, add = TRUE)
abline(0,1, lty = 3)
legend(x = 0.7, y = 0.4, legend = c('most_sensitive (AUC = 0.80)', 'least_sensitive (AUC = 0.59)'), lty = c(1,1), col = c('red', 'blue'), bty = 'n')
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