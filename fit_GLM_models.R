### THIS SCRIPT FITS GLM MODELS

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

## GDSC ------
# load clinical data ----
bleomycin     <- read.csv('Processed_Clinical_Data/bleomycin_gdsc_clinical_processed.csv', row.names = 1)
camptothecin     <- read.csv('Processed_Clinical_Data/camptothecin_gdsc_clinical_processed.csv', row.names = 1)
cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)
cytarabine     <- read.csv('Processed_Clinical_Data/cytarabine_gdsc_clinical_processed.csv', row.names = 1)
doxorubicin     <- read.csv('Processed_Clinical_Data/doxorubicin_gdsc_clinical_processed.csv', row.names = 1)
etoposide     <- read.csv('Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv', row.names = 1)
gemcitabine     <- read.csv('Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv', row.names = 1)
methotrexate  <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)
mitomycin     <- read.csv('Processed_Clinical_Data/mitomycin_gdsc_clinical_processed.csv', row.names = 1)
sn38     <- read.csv('Processed_Clinical_Data/sn38_gdsc_clinical_processed.csv', row.names = 1)
temozolomide     <- read.csv('Processed_Clinical_Data/temozolomide_gdsc_clinical_processed.csv', row.names = 1)
#
## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names
#
#
#
#
### set up data for model building ----------
# get names of GDSC cell lines treated with each drug
bleomycin_lines           <- bleomycin$COSMIC_ID #766
camptothecin_lines        <- camptothecin$COSMIC_ID #678
cisplatin_lines           <- cisplatin$COSMIC_ID #680
cytarabine_lines          <- cytarabine$COSMIC_ID #676
doxorubicin_lines         <- doxorubicin$COSMIC_ID #711
etoposide_lines           <- etoposide$COSMIC_ID #718
gemcitabine_lines         <- gemcitabine$COSMIC_ID #707
methotrexate_lines        <- methotrexate$COSMIC_ID # 679
mitomycin_lines           <- mitomycin$COSMIC_ID #712
sn38_lines                <- sn38$COSMIC_ID #761
temozolomide_lines        <- temozolomide$COSMIC_ID #739
#
#
#
# # get split indices
# set.seed(5)
# bleomycin_train_index <- createDataPartition(bleomycin$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# camptothecin_train_index <- createDataPartition(camptothecin$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# cisplatin_train_index <- createDataPartition(cisplatin$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# cytarabine_train_index <- createDataPartition(cytarabine$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# doxorubicin_train_index <- createDataPartition(doxorubicin$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# etoposide_train_index <- createDataPartition(etoposide$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# gemcitabine_train_index <- createDataPartition(gemcitabine$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# methotrexate_train_index <- createDataPartition(methotrexate$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# mitomycin_train_index <- createDataPartition(mitomycin$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# sn38_train_index <- createDataPartition(sn38$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# temozolomide_train_index <- createDataPartition(temozolomide$Cell_line_tissue_type, p = .8,
#                                   list = FALSE,
#                                   times = 1)
#
# # split clinical data
# bleomycin_train <- bleomycin[bleomycin_train_index, ]
# bleomycin_test <- bleomycin[-bleomycin_train_index, ]
#
# camptothecin_train <- camptothecin[camptothecin_train_index, ]
# camptothecin_test <- camptothecin[-camptothecin_train_index, ]
#
# cisplatin_train <- cisplatin[cisplatin_train_index, ]
# cisplatin_test <- cisplatin[-cisplatin_train_index, ]
#
# cytarabine_train <- cytarabine[cytarabine_train_index, ]
# cytarabine_test <- cytarabine[-cytarabine_train_index,]
#
# doxorubicin_train <- doxorubicin[doxorubicin_train_index, ]
# doxorubicin_test <- doxorubicin[-doxorubicin_train_index, ]
#
# etoposide_train <- etoposide[etoposide_train_index, ]
# etoposide_test <- etoposide[-etoposide_train_index, ]
#
# gemcitabine_train <- gemcitabine[gemcitabine_train_index, ]
# gemcitabine_test <- gemcitabine[-gemcitabine_train_index, ]
#
# methotrexate_train <- methotrexate[methotrexate_train_index, ]
# methotrexate_test <- methotrexate[-methotrexate_train_index, ]
#
# mitomycin_train <- mitomycin[mitomycin_train_index, ]
# mitomycin_test <- mitomycin[-mitomycin_train_index, ]
#
# sn38_train <- sn38[sn38_train_index, ]
# sn38_test <- sn38[-sn38_train_index, ]
#
# temozolomide_train <- temozolomide[temozolomide_train_index, ]
# temozolomide_test <- temozolomide[temozolomide_train_index, ]
#
#
# split expression data to match
bleomycin_rna_seq <- gdsc[rownames(gdsc) %in% bleomycin$COSMIC_ID, ]
bleomycin_rna_seq <- as.data.frame(bleomycin_rna_seq)
#bleomycin_rna_seq_test <- gdsc[rownames(gdsc) %in% bleomycin_test$COSMIC_ID, ]
#bleomycin_rna_seq_test <- as.data.frame(bleomycin_rna_seq_test)

bleomycin_rna_seq$res_sens <- bleomycin$res_sens
bleomycin_rose <- ROSE(res_sens ~ ., data = bleomycin_rna_seq)$data
#bleomycin_rose_res_sens <- bleomycin_rose$res_sens
#bleomycin_rose <- bleomycin_rose[, -14210]

camptothecin_rna_seq <- gdsc[rownames(gdsc) %in% camptothecin$COSMIC_ID, ]
camptothecin_rna_seq <- as.data.frame(camptothecin_rna_seq)
#camptothecin_rna_seq_test <- gdsc[rownames(gdsc) %in% camptothecin_test$COSMIC_ID, ]
#camptothecin_rna_seq_test <- as.data.frame(camptothecin_rna_seq_test)

camptothecin_rna_seq$res_sens <- camptothecin$res_sens
camptothecin_rose <- ROSE(res_sens ~ ., data = camptothecin_rna_seq)$data
#camptothecin_rose_res_sens <- camptothecin_rose$res_sens
#camptothecin_rose <- camptothecin_rose[, -14210]

cisplatin_rna_seq <- gdsc[rownames(gdsc) %in% cisplatin$COSMIC_ID, ]
cisplatin_rna_seq <- as.data.frame(cisplatin_rna_seq)
#cisplatin_rna_seq_test <- gdsc[rownames(gdsc) %in% cisplatin_test$COSMIC_ID, ]
#cisplatin_rna_seq_test <- as.data.frame(cisplatin_rna_seq_test)

cisplatin_rna_seq$res_sens <- cisplatin$res_sens
cisplatin_rose <- ROSE(res_sens ~ ., data = cisplatin_rna_seq)$data
#cisplatin_rose_res_sens <- cisplatin_rose$res_sens
#cisplatin_rose <- cisplatin_rose[, -14210]

cytarabine_rna_seq <- gdsc[rownames(gdsc) %in% cytarabine$COSMIC_ID, ]
cytarabine_rna_seq <- as.data.frame(cytarabine_rna_seq)
#cytarabine_rna_seq_test <- gdsc[rownames(gdsc) %in% cytarabine_test$COSMIC_ID, ]
#cytarabine_rna_seq_test <- as.data.frame(cytarabine_rna_seq_test)

cytarabine_rna_seq$res_sens <- cytarabine$res_sens
cytarabine_rose <- ROSE(res_sens ~ ., data = cytarabine_rna_seq)$data
#cytarabine_rose_res_sens <- cytarabine_rose$res_sens
#cytarabine_rose <- cytarabine_rose[, -14210]

doxorubicin_rna_seq <- gdsc[rownames(gdsc) %in% doxorubicin$COSMIC_ID, ]
doxorubicin_rna_seq <- as.data.frame(doxorubicin_rna_seq)
#doxorubicin_rna_seq_test <- gdsc[rownames(gdsc) %in% doxorubicin_test$COSMIC_ID, ]
#doxorubicin_rna_seq_test <- as.data.frame(doxorubicin_rna_seq_test)

doxorubicin_rna_seq$res_sens <- doxorubicin$res_sens
doxorubicin_rose <- ROSE(res_sens ~ ., data = doxorubicin_rna_seq)$data
#doxorubicin_rose_res_sens <- doxorubicin_rose$res_sens
#doxorubicin_rose <- doxorubicin_rose[, -14210]

etoposide_rna_seq <- gdsc[rownames(gdsc) %in% etoposide$COSMIC_ID, ]
etoposide_rna_seq <- as.data.frame(etoposide_rna_seq)
#etoposide_rna_seq_test <- gdsc[rownames(gdsc) %in% etoposide_test$COSMIC_ID, ]
#etoposide_rna_seq_test <- as.data.frame(etoposide_rna_seq_test)

etoposide_rna_seq$res_sens <- etoposide$res_sens
etoposide_rose <- ROSE(res_sens ~ ., data = etoposide_rna_seq)$data
#etoposide_rose_res_sens <- etoposide_rose$res_sens
#etoposide_rose <- etoposide_rose[, -14210]

gemcitabine_rna_seq <- gdsc[rownames(gdsc) %in% gemcitabine$COSMIC_ID, ]
gemcitabine_rna_seq <- as.data.frame(gemcitabine_rna_seq)
#gemcitabine_rna_seq_test <- gdsc[rownames(gdsc) %in% gemcitabine_test$COSMIC_ID, ]
#gemcitabine_rna_seq_test <- as.data.frame(gemcitabine_rna_seq_test)

gemcitabine_rna_seq$res_sens <- gemcitabine$res_sens
gemcitabine_rose <- ROSE(res_sens ~ ., data = gemcitabine_rna_seq)$data
#gemcitabine_rose_res_sens <- gemcitabine_rose$res_sens
#gemcitabine_rose <- gemcitabine_rose[, -14210]

methotrexate_rna_seq <- gdsc[rownames(gdsc) %in% methotrexate$COSMIC_ID, ]
methotrexate_rna_seq <- as.data.frame(methotrexate_rna_seq)
#methotrexate_rna_seq_test <- gdsc[rownames(gdsc) %in% methotrexate_test$COSMIC_ID, ]
#methotrexate_rna_seq_test <- as.data.frame(methotrexate_rna_seq_test)

methotrexate_rna_seq$res_sens <- methotrexate$res_sens
methotrexate_rose <- ROSE(res_sens ~ ., data = methotrexate_rna_seq)$data
#methotrexate_rose_res_sens <- methotrexate_rose$res_sens
#methotrexate_rose <- methotrexate_rose[, -14210]

mitomycin_rna_seq <- gdsc[rownames(gdsc) %in% mitomycin$COSMIC_ID, ]
mitomycin_rna_seq <- as.data.frame(mitomycin_rna_seq)
#mitomycin_rna_seq_test <- gdsc[rownames(gdsc) %in% mitomycin_test$COSMIC_ID, ]
#mitomycin_rna_seq_test <- as.data.frame(mitomycin_rna_seq_test)

mitomycin_rna_seq$res_sens <- mitomycin$res_sens
mitomycin_rose <- ROSE(res_sens ~ ., data = mitomycin_rna_seq)$data
#mitomycin_rose_res_sens <- mitomycin_rose$res_sens
#mitomycin_rose <- mitomycin_rose[, -14210]

sn38_rna_seq <- gdsc[rownames(gdsc) %in% sn38$COSMIC_ID, ]
sn38_rna_seq <- as.data.frame(sn38_rna_seq)
#sn38_rna_seq_test <- gdsc[rownames(gdsc) %in% sn38_test$COSMIC_ID, ]
#sn38_rna_seq_test <- as.data.frame(sn38_rna_seq_test)

sn38_rna_seq$res_sens <- sn38$res_sens
sn38_rose <- ROSE(res_sens ~ ., data = sn38_rna_seq)$data
#sn38_rose_res_sens <- sn38_rose$res_sens
#sn38_rose <- sn38_rose[, -14210]

temozolomide_rna_seq <- gdsc[rownames(gdsc) %in% temozolomide$COSMIC_ID, ]
temozolomide_rna_seq <- as.data.frame(temozolomide_rna_seq)
#temozolomide_rna_seq_test <- gdsc[rownames(gdsc) %in% temozolomide_test$COSMIC_ID, ]
#temozolomide_rna_seq_test <- as.data.frame(temozolomide_rna_seq_test)

temozolomide_rna_seq$res_sens <- temozolomide$res_sens
temozolomide_rose <- ROSE(res_sens ~ ., data = temozolomide_rna_seq)$data
#temozolomide_rose_res_sens <- temozolomide_rose$res_sens
#temozolomide_rose <- temozolomide_rose[, -14210]

#
# temozolomide_rose$res_sens <- temozolomide_rose_res_sens


write.csv(doxorubicin_rose, file = 'Processed_Gene_Expression/doxorubicin_rose_full.csv', row.names = TRUE)

# write.csv(temozolomide_test, file = 'Processed_Clinical_Data/temozolomide_test.csv', row.names = TRUE)
# write.csv(temozolomide_rna_seq_test, file = 'Processed_Gene_Expression/temozolomide_rna_seq_test.csv', row.names = TRUE)




# # split GDSC randomly
# set.seed(5)
# # get random numbers to use for split
# random_sample <- sample(x = rownames(gdsc), size = nrow(gdsc)/5)
#
# # get training and testing sets
# gdsc_train         <- gdsc[which(rownames(gdsc) %ni% random_sample), ] #770 x 14209
#
# gdsc_test          <- gdsc[random_sample, ] #192 x 14209
#
# # make sure zero overlap
# intersect(rownames(gdsc_train), rownames(gdsc_test)) #0

# create training/testing sets for each drug
# bleomycin_rna_seq_train     <- gdsc_train[intersect(bleomycin_lines, rownames(gdsc_train)), ]
# bleomycin_rna_seq_train_df     <- as.data.frame(bleomycin_rna_seq_train)
# # 622 x 14209
# bleomycin_rna_seq_test      <- gdsc_test[intersect(bleomycin_lines, rownames(gdsc_test)), ]
# bleomycin_rna_seq_test      <- as.data.frame(bleomycin_rna_seq_test)
# # 144 x 14209
# camptothecin_rna_seq_train     <- gdsc_train[intersect(camptothecin_lines, rownames(gdsc_train)), ]
# camptothecin_rna_seq_train_df <- as.data.frame(camptothecin_rna_seq_train)
# # 549 x 14209
# camptothecin_rna_seq_test      <- gdsc_test[intersect(camptothecin_lines, rownames(gdsc_test)), ]
# camptothecin_rna_seq_test <- as.data.frame(camptothecin_rna_seq_test)
# # 129 x 14209
# cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# cisplatin_rna_seq_train_df <- as.data.frame(cisplatin_rna_seq_train)
# # 551 x 14209
# cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]
# cisplatin_rna_seq_test <- as.data.frame(cisplatin_rna_seq_test)
# # 129 x 14209
# cytarabine_rna_seq_train     <- gdsc_train[intersect(cytarabine_lines, rownames(gdsc_train)), ]
# cytarabine_rna_seq_train_df <- as.data.frame(cytarabine_rna_seq_train)
# # 547 x 14209
# cytarabine_rna_seq_test      <- gdsc_test[intersect(cytarabine_lines, rownames(gdsc_test)), ]
# cytarabine_rna_seq_test <- as.data.frame(cytarabine_rna_seq_test)
# # 129 x 14209
# doxorubicin_rna_seq_train     <- gdsc_train[intersect(doxorubicin_lines, rownames(gdsc_train)), ]
# doxorubicin_rna_seq_train_df <- as.data.frame(doxorubicin_rna_seq_train)
# # 581 x 14209
# doxorubicin_rna_seq_test      <- gdsc_test[intersect(doxorubicin_lines, rownames(gdsc_test)), ]
# doxorubicin_rna_seq_test <- as.data.frame(doxorubicin_rna_seq_test)
# # 130 x 14209
# etoposide_rna_seq_train     <- gdsc_train[intersect(etoposide_lines, rownames(gdsc_train)), ]
# etoposide_rna_seq_train_df <- as.data.frame(etoposide_rna_seq_train)
# # 584 x 14209
# etoposide_rna_seq_test      <- gdsc_test[intersect(etoposide_lines, rownames(gdsc_test)), ]
# etoposide_rna_seq_test <- as.data.frame(etoposide_rna_seq_test)
# # 134 x 14209
# gemcitabine_rna_seq_train     <- gdsc_train[intersect(gemcitabine_lines, rownames(gdsc_train)), ]
# gemcitabine_rna_seq_train_df <- as.data.frame(gemcitabine_rna_seq_train)
# # 577 x 14209
# gemcitabine_rna_seq_test      <- gdsc_test[intersect(gemcitabine_lines, rownames(gdsc_test)), ]
# gemcitabine_rna_seq_test <- as.data.frame(gemcitabine_rna_seq_test)
# # 130 x 14209
# methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# methotrexate_rna_seq_train_df <- as.data.frame(methotrexate_rna_seq_train)
# # 550 x 14209
# methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# methotrexate_rna_seq_test <- as.data.frame(methotrexate_rna_seq_test)
# # 129 x 14209
# mitomycin_rna_seq_train     <- gdsc_train[intersect(mitomycin_lines, rownames(gdsc_train)), ]
# mitomycin_rna_seq_train_df <- as.data.frame(mitomycin_rna_seq_train)
# # 581 x 14209
# mitomycin_rna_seq_test      <- gdsc_test[intersect(mitomycin_lines, rownames(gdsc_test)), ]
# mitomycin_rna_seq_test <- as.data.frame(mitomycin_rna_seq_test)
# # 131 x 14209
# sn38_rna_seq_train     <- gdsc_train[intersect(sn38_lines, rownames(gdsc_train)), ]
# sn38_rna_seq_train_df <- as.data.frame(sn38_rna_seq_train)
# # 619 x 14209
# sn38_rna_seq_test      <- gdsc_test[intersect(sn38_lines, rownames(gdsc_test)), ]
# sn38_rna_seq_test <- as.data.frame(sn38_rna_seq_test)
# # 142 x 14209
# temozolomide_rna_seq_train     <- gdsc_train[intersect(temozolomide_lines, rownames(gdsc_train)), ]
# temozolomide_rna_seq_train_df <- as.data.frame(temozolomide_rna_seq_train)
# # 601 x 14209
# temozolomide_rna_seq_test      <- gdsc_test[intersect(temozolomide_lines, rownames(gdsc_test)), ]
# temozolomide_rna_seq_test <- as.data.frame(temozolomide_rna_seq_test)
# # 138 x 14209

# split clinical data
# bleomycin_train        <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_train)), ]
# bleomycin_rna_seq_train_df$res_sens <- bleomycin_train$res_sens
# bleomycin_rose <- ROSE(res_sens ~ ., data = bleomycin_rna_seq_train_df)$data
# bleomycin_rose_res_sens <- bleomycin_rose$res_sens
# bleomycin_rose <- bleomycin_rose[, -14210]
#
# bleomycin_test         <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_test)), ]
#
# camptothecin_train        <- camptothecin[which(camptothecin$COSMIC_ID %in% rownames(camptothecin_rna_seq_train)), ]
# camptothecin_rna_seq_train_df$res_sens <- camptothecin_train$res_sens
# camptothecin_rose <- ROSE(res_sens ~ ., data = camptothecin_rna_seq_train_df)$data
# camptothecin_rose_res_sens <- camptothecin_rose$res_sens
# camptothecin_rose <- camptothecin_rose[, -14210]
#
# camptothecin_test         <- camptothecin[which(camptothecin$COSMIC_ID %in% rownames(camptothecin_rna_seq_test)), ]
#
# cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
# cisplatin_rna_seq_train_df$res_sens <- cisplatin_train$res_sens
# cisplatin_rose <- ROSE(res_sens ~ ., data = cisplatin_rna_seq_train_df)$data
# cisplatin_rose_res_sens <- cisplatin_rose$res_sens
# cisplatin_rose <- cisplatin_rose[, -14210]
#
# cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]
#
# cytarabine_train        <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_train)), ]
# cytarabine_rna_seq_train_df$res_sens <- cytarabine_train$res_sens
# cytarabine_rose <- ROSE(res_sens ~ ., data = cytarabine_rna_seq_train_df)$data
# cytarabine_rose_res_sens <- cytarabine_rose$res_sens
# cytarabine_rose <- cytarabine_rose[, -14210]
#
# cytarabine_test         <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_test)), ]
#
# doxorubicin_train        <- doxorubicin[which(doxorubicin$COSMIC_ID %in% rownames(doxorubicin_rna_seq_train)), ]
# doxorubicin_rna_seq_train_df$res_sens <- doxorubicin_train$res_sens
# doxorubicin_rose <- ROSE(res_sens ~ ., data = doxorubicin_rna_seq_train_df)$data
# doxorubicin_rose_res_sens <- doxorubicin_rose$res_sens
# doxorubicin_rose <- doxorubicin_rose[, -14210]
#
# doxorubicin_test         <- doxorubicin[which(doxorubicin$COSMIC_ID %in% rownames(doxorubicin_rna_seq_test)), ]
#
# etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
# etoposide_rna_seq_train_df$res_sens <- etoposide_train$res_sens
# etoposide_rose <- ROSE(res_sens ~ ., data = etoposide_rna_seq_train_df)$data
# etoposide_rose_res_sens <- etoposide_rose$res_sens
# etoposide_rose <- etoposide_rose[, -14210]
#
# etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]
#
# gemcitabine_train        <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_train)), ]
# gemcitabine_rna_seq_train_df$res_sens <- gemcitabine_train$res_sens
# gemcitabine_rose <- ROSE(res_sens ~ ., data = gemcitabine_rna_seq_train_df)$data
# gemcitabine_rose_res_sens <- gemcitabine_rose$res_sens
# gemcitabine_rose <- gemcitabine_rose[, -14210]
#
# gemcitabine_test         <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_test)), ]
#
# methotrexate_train        <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
# methotrexate_rna_seq_train_df$res_sens <- methotrexate_train$res_sens
# methotrexate_rose <- ROSE(res_sens ~ ., data = methotrexate_rna_seq_train_df)$data
# methotrexate_rose_res_sens <- methotrexate_rose$res_sens
# methotrexate_rose <- methotrexate_rose[, -14210]
#
# methotrexate_test         <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]
#
# mitomycin_train        <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_train)), ]
# mitomycin_rna_seq_train_df$res_sens <- mitomycin_train$res_sens
# mitomycin_rose <- ROSE(res_sens ~ ., data = mitomycin_rna_seq_train_df)$data
# mitomycin_rose_res_sens <- mitomycin_rose$res_sens
# mitomycin_rose <- mitomycin_rose[, -14210]
#
# mitomycin_test         <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_test)), ]
#
# sn38_train        <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_train)), ]
# sn38_rna_seq_train_df$res_sens <- sn38_train$res_sens
# sn38_rose <- ROSE(res_sens ~ ., data = sn38_rna_seq_train_df)$data
# sn38_rose_res_sens <- sn38_rose$res_sens
# sn38_rose <- sn38_rose[, -14210]
#
# sn38_test         <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_test)), ]
#
# temozolomide_train        <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_train)), ]
# temozolomide_rna_seq_train_df$res_sens <- temozolomide_train$res_sens
# temozolomide_rose <- ROSE(res_sens ~ ., data = temozolomide_rna_seq_train_df)$data
# temozolomide_rose_res_sens <- temozolomide_rose$res_sens
# temozolomide_rose <- temozolomide_rose[, -14210]
#
# temozolomide_test         <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_test)), ]



### fit models --------
#training sequence data
bleomycin_rose <- read.csv('Processed_Gene_Expression/bleomycin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#bleomycin_test <- read.csv('Processed_Clinical_Data/bleomycin_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#bleomycin_rna_seq_test <- read.csv('Processed_Gene_Expression/bleomycin_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
camptothecin_rose <- read.csv('Processed_Gene_Expression/camptothecin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#camptothecin_test <- read.csv('Processed_Clinical_Data/camptothecin_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#camptothecin_rna_seq_test <- read.csv('Processed_Gene_Expression/camptothecin_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
cisplatin_rose <- read.csv('Processed_Gene_Expression/cisplatin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#cisplatin_test <- read.csv('Processed_Clinical_Data/cisplatin_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#cisplatin_rna_seq_test <- read.csv('Processed_Gene_Expression/cisplatin_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
cytarabine_rose <- read.csv('Processed_Gene_Expression/cytarabine_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#cytarabine_test <- read.csv('Processed_Clinical_Data/cytarabine_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#cytarabine_rna_seq_test <- read.csv('Processed_Gene_Expression/cytarabine_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
doxorubicin_rose <- read.csv('Processed_Gene_Expression/doxorubicin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#doxorubicin_test <- read.csv('Processed_Clinical_Data/doxorubicin_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#doxorubicin_rna_seq_test <- read.csv('Processed_Gene_Expression/doxorubicin_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
etoposide_rose <- read.csv('Processed_Gene_Expression/etoposide_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#etoposide_test <- read.csv('Processed_Clinical_Data/etoposide_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#etoposide_rna_seq_test <- read.csv('Processed_Gene_Expression/etoposide_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
gemcitabine_rose <- read.csv('Processed_Gene_Expression/gemcitabine_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#gemcitabine_test <- read.csv('Processed_Clinical_Data/gemcitabine_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#gemcitabine_rna_seq_test <- read.csv('Processed_Gene_Expression/gemcitabine_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
methotrexate_rose <- read.csv('Processed_Gene_Expression/methotrexate_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#methotrexate_test <- read.csv('Processed_Clinical_Data/methotrexate_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#methotrexate_rna_seq_test <- read.csv('Processed_Gene_Expression/methotrexate_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
mitomycin_rose <- read.csv('Processed_Gene_Expression/mitomycin_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#mitomycin_test <- read.csv('Processed_Clinical_Data/mitomycin_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#mitomycin_rna_seq_test <- read.csv('Processed_Gene_Expression/mitomycin_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
sn38_rose <- read.csv('Processed_Gene_Expression/sn38_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#sn38_test <- read.csv('Processed_Clinical_Data/sn38_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#sn38_rna_seq_test <- read.csv('Processed_Gene_Expression/sn38_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#training sequence data
temozolomide_rose <- read.csv('Processed_Gene_Expression/temozolomide_rose_full.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

#testing clinical data
#temozolomide_test <- read.csv('Processed_Clinical_Data/temozolomide_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)
#testing sequence data
#temozolomide_rna_seq_test <- read.csv('Processed_Gene_Expression/temozolomide_rna_seq_test.csv', stringsAsFactors = FALSE, row.names = 1, header = TRUE)

## fit models ----

## bleomycin ----
## glm model
set.seed(5)
bleomycin_fit_elnet <- cv.glmnet(x = as.matrix(bleomycin_rose[, -14210]), y = bleomycin_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/bleomycin_auc.png')
plot(bleomycin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/bleomycin_model.rds', bleomycin_fit_elnet)
bleomycin_pred <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#bleomycin_preds <- glmnet::auc(bleomycin$res_sens, bleomycin_pred)
#bleomycin_preds <- round(bleomycin_preds, digits = 2) #0.66

bleomycin_overall_acc <- sum(bleomycin$res_sens == bleomycin_pred)/length(bleomycin_pred) #0.916

bleomycin_glm_confusion <- confusionMatrix(factor(bleomycin_pred), factor(bleomycin$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

bleomycin_svm <- train(factor(res_sens) ~ ., data = bleomycin_rose, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)

#ggplot(svm_Linear)



bleomycin_pred <- predict(bleomycin_svm, newdata = bleomycin_rna_seq)


bleomycin_svm_confusion <- confusionMatrix(bleomycin_pred, factor(bleomycin$res_sens))

saveRDS(file = 'GLM_Models/bleomycin_svm_model_1.rds', bleomycin_svm)



## fit glm mars first degree
set.seed(5)
bleomycin_fit_1 <- earth(res_sens ~ ., data = bleomycin_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
bleomycin_fit_1_summary <- summary(bleomycin_fit_1)
bleomycin_fit_1_evimp <- evimp(bleomycin_fit_1) #plot these w ggplot side by side

bleomycin_fit_1_pred <- predict(bleomycin_fit_1, newdata = as.matrix(bleomycin_rna_seq), type = 'class')

bleomycin_mars_confusion <- confusionMatrix(factor(bleomycin_fit_1_pred), factor(bleomycin$res_sens))


saveRDS(file = 'GLM_Models/bleomycin_cv_mars_glm_model_1.rds', bleomycin_fit_1)


## camptothecin ----
## glm model
set.seed(5)
camptothecin_fit_elnet <- cv.glmnet(x = as.matrix(camptothecin_rose[, -14210]), y = camptothecin_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/camptothecin_auc.png')
plot(camptothecin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/camptothecin_model.rds', camptothecin_fit_elnet)
camptothecin_pred <- predict(camptothecin_fit_elnet, newx = as.matrix(camptothecin_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#camptothecin_preds <- glmnet::auc(camptothecin$res_sens, camptothecin_pred)
#camptothecin_preds <- round(camptothecin_preds, digits = 2) #0.66

camptothecin_overall_acc <- sum(camptothecin$res_sens == camptothecin_pred)/length(camptothecin_pred) #0.916

camptothecin_glm_confusion <- confusionMatrix(factor(camptothecin_pred), factor(camptothecin$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

camptothecin_svm <- train(factor(res_sens) ~ ., data = camptothecin_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



camptothecin_pred <- predict(camptothecin_svm, newdata = camptothecin_rna_seq)


camptothecin_svm_confusion <- confusionMatrix(camptothecin_pred, factor(camptothecin$res_sens))

saveRDS(file = 'GLM_Models/camptothecin_svm_model_1.rds', camptothecin_svm)



## fit glm mars first degree
set.seed(5)
camptothecin_fit_1 <- earth(res_sens ~ ., data = camptothecin_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
camptothecin_fit_1_summary <- summary(camptothecin_fit_1)
camptothecin_fit_1_evimp <- evimp(camptothecin_fit_1) #plot these w ggplot side by side

camptothecin_fit_1_pred <- predict(camptothecin_fit_1, newdata = as.matrix(camptothecin_rna_seq), type = 'class')

camptothecin_mars_confusion <- confusionMatrix(factor(camptothecin_fit_1_pred), factor(camptothecin$res_sens))


saveRDS(file = 'GLM_Models/camptothecin_cv_mars_glm_model_1.rds', camptothecin_fit_1)



## cisplatin ----
## glm model
set.seed(5)
cisplatin_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rose[, -14210]), y = cisplatin_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/cisplatin_auc.png')
plot(cisplatin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/cisplatin_model.rds', cisplatin_fit_elnet)
cisplatin_pred <- predict(cisplatin_fit_elnet, newx = as.matrix(cisplatin_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#cisplatin_preds <- glmnet::auc(cisplatin$res_sens, cisplatin_pred)
#cisplatin_preds <- round(cisplatin_preds, digits = 2) #0.66

cisplatin_overall_acc <- sum(cisplatin$res_sens == cisplatin_pred)/length(cisplatin_pred) #0.916

cisplatin_glm_confusion <- confusionMatrix(factor(cisplatin_pred), factor(cisplatin$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

cisplatin_svm <- train(factor(res_sens) ~ ., data = cisplatin_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



cisplatin_pred <- predict(cisplatin_svm, newdata = cisplatin_rna_seq)


cisplatin_svm_confusion <- confusionMatrix(cisplatin_pred, factor(cisplatin$res_sens))

saveRDS(file = 'GLM_Models/cisplatin_svm_model_1.rds', cisplatin_svm)



## fit glm mars first degree
set.seed(5)
cisplatin_fit_1 <- earth(res_sens ~ ., data = cisplatin_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
cisplatin_fit_1_summary <- summary(cisplatin_fit_1)
cisplatin_fit_1_evimp <- evimp(cisplatin_fit_1) #plot these w ggplot side by side

cisplatin_fit_1_pred <- predict(cisplatin_fit_1, newdata = as.matrix(cisplatin_rna_seq), type = 'class')

cisplatin_mars_confusion <- confusionMatrix(factor(cisplatin_fit_1_pred), factor(cisplatin$res_sens))


saveRDS(file = 'GLM_Models/cisplatin_cv_mars_glm_model_1.rds', cisplatin_fit_1)



## cytarabine ----
## glm model
set.seed(5)
cytarabine_fit_elnet <- cv.glmnet(x = as.matrix(cytarabine_rose[, -14210]), y = cytarabine_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/cytarabine_auc.png')
plot(cytarabine_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/cytarabine_model.rds', cytarabine_fit_elnet)
cytarabine_pred <- predict(cytarabine_fit_elnet, newx = as.matrix(cytarabine_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#cytarabine_preds <- glmnet::auc(cytarabine$res_sens, cytarabine_pred)
#cytarabine_preds <- round(cytarabine_preds, digits = 2) #0.66

cytarabine_overall_acc <- sum(cytarabine$res_sens == cytarabine_pred)/length(cytarabine_pred) #0.916

cytarabine_glm_confusion <- confusionMatrix(factor(cytarabine_pred), factor(cytarabine$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

cytarabine_svm <- train(factor(res_sens) ~ ., data = cytarabine_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



cytarabine_pred <- predict(cytarabine_svm, newdata = cytarabine_rna_seq)


cytarabine_svm_confusion <- confusionMatrix(cytarabine_pred, factor(cytarabine$res_sens))

saveRDS(file = 'GLM_Models/cytarabine_svm_model_1.rds', cytarabine_svm)



## fit glm mars first degree
set.seed(5)
cytarabine_fit_1 <- earth(res_sens ~ ., data = cytarabine_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
cytarabine_fit_1_summary <- summary(cytarabine_fit_1)
cytarabine_fit_1_evimp <- evimp(cytarabine_fit_1) #plot these w ggplot side by side

cytarabine_fit_1_pred <- predict(cytarabine_fit_1, newdata = as.matrix(cytarabine_rna_seq), type = 'class')

cytarabine_mars_confusion <- confusionMatrix(factor(cytarabine_fit_1_pred), factor(cytarabine$res_sens))


saveRDS(file = 'GLM_Models/cytarabine_cv_mars_glm_model_1.rds', cytarabine_fit_1)



## doxorubicin ----
## glm model
set.seed(5)
doxorubicin_fit_elnet <- cv.glmnet(x = as.matrix(doxorubicin_rose[, -14210]), y = doxorubicin_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/doxorubicin_auc.png')
plot(doxorubicin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/doxorubicin_model.rds', doxorubicin_fit_elnet)
doxorubicin_pred <- predict(doxorubicin_fit_elnet, newx = as.matrix(doxorubicin_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#doxorubicin_preds <- glmnet::auc(doxorubicin$res_sens, doxorubicin_pred)
#doxorubicin_preds <- round(doxorubicin_preds, digits = 2) #0.66

doxorubicin_overall_acc <- sum(doxorubicin$res_sens == doxorubicin_pred)/length(doxorubicin_pred) #0.916

doxorubicin_glm_confusion <- confusionMatrix(factor(doxorubicin_pred), factor(doxorubicin$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

doxorubicin_svm <- train(factor(res_sens) ~ ., data = doxorubicin_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



doxorubicin_pred <- predict(doxorubicin_svm, newdata = doxorubicin_rna_seq)


doxorubicin_svm_confusion <- confusionMatrix(doxorubicin_pred, factor(doxorubicin$res_sens))

saveRDS(file = 'GLM_Models/doxorubicin_svm_model_1.rds', doxorubicin_svm)



## fit glm mars first degree
set.seed(5)
doxorubicin_fit_1 <- earth(res_sens ~ ., data = doxorubicin_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
doxorubicin_fit_1_summary <- summary(doxorubicin_fit_1)
doxorubicin_fit_1_evimp <- evimp(doxorubicin_fit_1) #plot these w ggplot side by side

doxorubicin_fit_1_pred <- predict(doxorubicin_fit_1, newdata = as.matrix(doxorubicin_rna_seq), type = 'class')

doxorubicin_mars_confusion <- confusionMatrix(factor(doxorubicin_fit_1_pred), factor(doxorubicin$res_sens))


saveRDS(file = 'GLM_Models/doxorubicin_cv_mars_glm_model_1.rds', doxorubicin_fit_1)



## etoposide ----
## glm model
set.seed(5)
etoposide_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rose[, -14210]), y = etoposide_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/etoposide_auc.png')
plot(etoposide_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/etoposide_model.rds', etoposide_fit_elnet)
etoposide_pred <- predict(etoposide_fit_elnet, newx = as.matrix(etoposide_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#etoposide_preds <- glmnet::auc(etoposide$res_sens, etoposide_pred)
#etoposide_preds <- round(etoposide_preds, digits = 2) #0.66

etoposide_overall_acc <- sum(etoposide$res_sens == etoposide_pred)/length(etoposide_pred) #0.916

etoposide_glm_confusion <- confusionMatrix(factor(etoposide_pred), factor(etoposide$res_sens))

## svm model
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
set.seed(5)

etoposide_svm <- train(factor(res_sens) ~ ., data = etoposide_rose, method = "svmLinear",
                       trControl=trctrl,
                       preProcess = c("center", "scale"),
                       tuneLength = 10)

#ggplot(svm_Linear)



etoposide_pred <- predict(etoposide_svm, newdata = etoposide_rna_seq)


etoposide_svm_confusion <- confusionMatrix(etoposide_pred, factor(etoposide$res_sens))

saveRDS(file = 'GLM_Models/etoposide_svm_model_1.rds', etoposide_svm)



## fit glm mars first degree
set.seed(5)
etoposide_fit_1 <- earth(res_sens ~ ., data = etoposide_rose, ncross = 5, degree=1, nfold=5, pmethod = 'cv', keepxy=TRUE, glm=list(family=binomial), trace = .5) 
etoposide_fit_1_summary <- summary(etoposide_fit_1)
etoposide_fit_1_evimp <- evimp(etoposide_fit_1) #plot these w ggplot side by side

etoposide_fit_1_pred <- predict(etoposide_fit_1, newdata = as.matrix(etoposide_rna_seq), type = 'class')

etoposide_mars_confusion <- confusionMatrix(factor(etoposide_fit_1_pred), factor(etoposide$res_sens))


saveRDS(file = 'GLM_Models/etoposide_cv_mars_glm_model_1.rds', etoposide_fit_1)



## gemcitabine ----
## glm model
set.seed(5)
gemcitabine_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rose[, -14210]), y = gemcitabine_rose$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/gemcitabine_auc.png')
plot(gemcitabine_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/gemcitabine_model.rds', gemcitabine_fit_elnet)
gemcitabine_pred <- predict(gemcitabine_fit_elnet, newx = as.matrix(gemcitabine_rna_seq[, -14210]), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
#gemcitabine_preds <- glmnet::auc(gemcitabine$res_sens, gemcitabine_pred)
#gemcitabine_preds <- round(gemcitabine_preds, digits = 2) #0.66

gemcitabine_overall_acc <- sum(gemcitabine$res_sens == gemcitabine_pred)/length(gemcitabine_pred) #0.916

gemcitabine_glm_confusion <- confusionMatrix(factor(gemcitabine_pred), factor(gemcitabine$res_sens))

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



## methotrexate ----
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



## mitomycin ----
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



## sn38 ----
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



## temozolomide ----
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
















## fitting MARS models
# fit a basic MARS model
bleomycin_mars1 <- earth(res_sens ~ ., data = bleomycin_rose)
#print model summary
print(bleomycin_mars1) #23/10
summary(bleomycin_mars1)
coef(bleomycin_mars1)
evimp(bleomycin_mars1)

png(filename = 'Images/bleomycin_first_degree_model_selection.png')
plot(bleomycin_mars1, which = 1)
dev.off()


bleomycin_mars1_predictions <- predict(bleomycin_mars1, bleomycin_rna_seq_test)
mse <- mean((bleomycin_test$res_sens - bleomycin_mars1_predictions) ^ 2) #0.65

# fit a basic MARS model
bleomycin_mars2 <- earth(res_sens ~ ., data = bleomycin_rose, degree = 2)

print(bleomycin_mars2) #53
summary(bleomycin_mars2)

## create a tuning grid
hyper_grid <- expand.grid(degree = 1:3, nprune = seq(2,70,length.out = 10) %>% floor())

## cv mars model
set.seed(5)
bleomycin_cv_mars <- train(x = subset(bleomycin_rose, select = -res_sens), 
                 y = as.factor(bleomycin_rose$res_sens), 
                 method = 'earth', 
                 trControl = trainControl(method = 'cv', number = 10), 
                 tuneGrid = hyper_grid)

#results
bleomycin_cv_mars$bestTune

#plot it
png(filename = 'Images/bleomycin_cv_mars.png')
ggplot(bleomycin_cv_mars) + ggtitle('Bleomycin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/bleomycin_cv_mars_model.rds', bleomycin_cv_mars)

#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(bleomycin_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('Bleomycin GCV')

p2 <- vip(bleomycin_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('Bleomycin RSS')

png(filename = 'Images/bleomycin_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
bleomycin_svm <- ksvm(res_sens ~ ., bleomycin_rose)

summary(bleomycin_svm) # is there better info than this?

bleomycin_svm_predictions <- predict(bleomycin_svm, bleomycin_rna_seq_test)
mse <- mean((bleomycin_test$res_sens - bleomycin_svm_predictions) ^ 2)
mse #0.42


## kNN
bleomycin_knn <- knnreg(bleomycin_rose[ ,1:14209], bleomycin_rose[, 14210], k = 2)

summary(bleomycin_knn)
bleomycin_knn_predictions <- predict(bleomycin_knn, bleomycin_rna_seq_test)
mse <- mean((bleomycin_test$res_sens - bleomycin_knn_predictions) ^ 2)
mse #0.35


# fit a basic MARS model
camptothecin_mars1 <- earth(res_sens ~ ., data = camptothecin_rose)
#print model summary
print(camptothecin_mars1) #22
summary(camptothecin_mars1)
coef(camptothecin_mars1)

png(filename = 'Images/camptothecin_first_degree_model_selection.png')
plot(camptothecin_mars1, which = 1)
dev.off()


camptothecin_mars1_predictions <- predict(camptothecin_mars1, camptothecin_rna_seq_test)
mse <- mean((camptothecin_test$res_sens - camptothecin_mars1_predictions) ^ 2) #0.71

# fit a basic MARS model
camptothecin_mars2 <- earth(res_sens ~ ., data = camptothecin_rose, degree = 2)

print(camptothecin_mars2) #45
summary(camptothecin_mars2)



## cv mars model
set.seed(5)
camptothecin_cv_mars <- train(x = subset(camptothecin_rose, select = -res_sens),
                           y = as.factor(camptothecin_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
camptothecin_cv_mars$bestTune

#plot it
png(filename = 'Images/camptothecin_cv_mars.png')
ggplot(camptothecin_cv_mars) + ggtitle('Camptothecin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/camptothecin_cv_mars_model.rds', camptothecin_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(camptothecin_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('camptothecin GCV')

p2 <- vip(camptothecin_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('camptothecin RSS')

png(filename = 'Images/camptothecin_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
camptothecin_svm <- ksvm(res_sens ~ ., camptothecin_rose)

summary(camptothecin_svm) # is there better info than this?

camptothecin_svm_predictions <- predict(camptothecin_svm, camptothecin_rna_seq_test)
mse <- mean((camptothecin_test$res_sens - camptothecin_svm_predictions) ^ 2)
mse #0.40


## kNN
camptothecin_knn <- knnreg(camptothecin_rose[ ,1:14209], camptothecin_rose[, 14210], k = 2)

summary(camptothecin_knn)
camptothecin_knn_predictions <- predict(camptothecin_knn, camptothecin_rna_seq_test)
mse <- mean((camptothecin_test$res_sens - camptothecin_knn_predictions) ^ 2)
mse #0.33

# fit a basic MARS model
cisplatin_mars1 <- earth(res_sens ~ ., data = cisplatin_rose)
#print model summary
print(cisplatin_mars1) #26
summary(cisplatin_mars1)
coef(cisplatin_mars1)

png(filename = 'Images/cisplatin_first_degree_model_selection.png')
plot(cisplatin_mars1, which = 1)
dev.off()


cisplatin_mars1_predictions <- predict(cisplatin_mars1, cisplatin_rna_seq_test)
mse <- mean((cisplatin_test$res_sens - cisplatin_mars1_predictions) ^ 2) #0.85

# fit a basic MARS model
cisplatin_mars2 <- earth(res_sens ~ ., data = cisplatin_rose, degree = 2)

print(cisplatin_mars2) #44
summary(cisplatin_mars2)



## cv mars model
set.seed(5)
cisplatin_cv_mars <- train(x = subset(cisplatin_rose, select = -res_sens),
                           y = as.factor(cisplatin_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
cisplatin_cv_mars$bestTune

#plot it
png(filename = 'Images/cisplatin_cv_mars.png')
ggplot(cisplatin_cv_mars) + ggtitle('Cisplatin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/cisplatin_cv_mars_model.rds', cisplatin_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(cisplatin_cv_mars, num_features = 30, bar = FALSE, value = 'gcv') + ggtitle('cisplatin GCV')

p2 <- vip(cisplatin_cv_mars, num_features = 30, bar = FALSE, value = 'rss') + ggtitle('cisplatin RSS')

png(filename = 'Images/cisplatin_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
cisplatin_svm <- ksvm(res_sens ~ ., cisplatin_rose)

summary(cisplatin_svm) # is there better info than this?

cisplatin_svm_predictions <- predict(cisplatin_svm, cisplatin_rna_seq_test)
mse <- mean((cisplatin_test$res_sens - cisplatin_svm_predictions) ^ 2)
mse #0.45


## kNN
cisplatin_knn <- knnreg(cisplatin_rose[ ,1:14209], cisplatin_rose[, 14210], k = 2)

summary(cisplatin_knn)
cisplatin_knn_predictions <- predict(cisplatin_knn, cisplatin_rna_seq_test)
mse <- mean((cisplatin_test$res_sens - cisplatin_knn_predictions) ^ 2)
mse #0.40

# fit a basic MARS model
cytarabine_mars1 <- earth(res_sens ~ ., data = cytarabine_rose)
#print model summary
print(cytarabine_mars1) #22
summary(cytarabine_mars1)
coef(cytarabine_mars1)

png(filename = 'Images/cytarabine_first_degree_model_selection.png')
plot(cytarabine_mars1, which = 1)
dev.off()


cytarabine_mars1_predictions <- predict(cytarabine_mars1, cytarabine_rna_seq_test)
mse <- mean((cytarabine_test$res_sens - cytarabine_mars1_predictions) ^ 2) #0.66

# fit a basic MARS model
cytarabine_mars2 <- earth(res_sens ~ ., data = cytarabine_rose, degree = 2)

print(cytarabine_mars2) #47
summary(cytarabine_mars2)

## cv mars model
set.seed(5)
cytarabine_cv_mars <- train(x = subset(cytarabine_rose, select = -res_sens),
                           y = as.factor(cytarabine_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
cytarabine_cv_mars$bestTune

#plot it
png(filename = 'Images/cytarabine_cv_mars.png')
ggplot(cytarabine_cv_mars) + ggtitle('Cytarabine CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/cytarabine_cv_mars_model.rds', cytarabine_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(cytarabine_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('cytarabine GCV')

p2 <- vip(cytarabine_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('cytarabine RSS')

png(filename = 'Images/cytarabine_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


## svm
cytarabine_svm <- ksvm(res_sens ~ ., cytarabine_rose)

summary(cytarabine_svm) # is there better info than this?

cytarabine_svm_predictions <- predict(cytarabine_svm, cytarabine_rna_seq_test)
mse <- mean((cytarabine_test$res_sens - cytarabine_svm_predictions) ^ 2)
mse #0.46


## kNN
cytarabine_knn <- knnreg(cytarabine_rose[ ,1:14209], cytarabine_rose[, 14210], k = 2)

summary(cytarabine_knn)
cytarabine_knn_predictions <- predict(cytarabine_knn, cytarabine_rna_seq_test)
mse <- mean((cytarabine_test$res_sens - cytarabine_knn_predictions) ^ 2)
mse #0.37

# fit a basic MARS model
doxorubicin_mars1 <- earth(res_sens ~ ., data = doxorubicin_rose)
#print model summary
print(doxorubicin_mars1) #21
summary(doxorubicin_mars1)
coef(doxorubicin_mars1)

png(filename = 'Images/doxorubicin_first_degree_model_selection.png')
plot(doxorubicin_mars1, which = 1)
dev.off()


doxorubicin_mars1_predictions <- predict(doxorubicin_mars1, doxorubicin_rna_seq_test)
mse <- mean((doxorubicin_test$res_sens - doxorubicin_mars1_predictions) ^ 2) #0.73

# fit a basic MARS model
doxorubicin_mars2 <- earth(res_sens ~ ., data = doxorubicin_rose, degree = 2)

print(doxorubicin_mars2) #47
summary(doxorubicin_mars2)

## cv mars model
set.seed(5)
doxorubicin_cv_mars <- train(x = subset(doxorubicin_rose, select = -res_sens),
                           y = as.factor(doxorubicin_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
doxorubicin_cv_mars$bestTune

#plot it
png(filename = 'Images/doxorubicin_cv_mars.png')
ggplot(doxorubicin_cv_mars) + ggtitle('Doxorubicin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/doxorubicin_cv_mars_model.rds', doxorubicin_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(doxorubicin_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('doxorubicin GCV')

p2 <- vip(doxorubicin_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('doxorubicin RSS')

png(filename = 'Images/doxorubicin_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
doxorubicin_svm <- ksvm(res_sens ~ ., doxorubicin_rose)

summary(doxorubicin_svm) # is there better info than this?

doxorubicin_svm_predictions <- predict(doxorubicin_svm, doxorubicin_rna_seq_test)
mse <- mean((doxorubicin_test$res_sens - doxorubicin_svm_predictions) ^ 2)
mse #0.33


## kNN
doxorubicin_knn <- knnreg(doxorubicin_rose[ ,1:14209], doxorubicin_rose[, 14210], k = 2)

summary(doxorubicin_knn)
doxorubicin_knn_predictions <- predict(doxorubicin_knn, doxorubicin_rna_seq_test)
mse <- mean((doxorubicin_test$res_sens - doxorubicin_knn_predictions) ^ 2)
mse #0.27

# fit a basic MARS model
etoposide_mars1 <- earth(res_sens ~ ., data = etoposide_rose)
#print model summary
print(etoposide_mars1) #26
summary(etoposide_mars1)
coef(etoposide_mars1)

png(filename = 'Images/etoposide_first_degree_model_selection.png')
plot(etoposide_mars1, which = 1)
dev.off()


etoposide_mars1_predictions <- predict(etoposide_mars1, etoposide_rna_seq_test)
mse <- mean((etoposide_test$res_sens - etoposide_mars1_predictions) ^ 2) #0.82

# fit a basic MARS model
etoposide_mars2 <- earth(res_sens ~ ., data = etoposide_rose, degree = 2)

print(etoposide_mars2) #44
summary(etoposide_mars2)

## cv mars model
set.seed(5)
etoposide_cv_mars <- train(x = subset(etoposide_rose, select = -res_sens),
                           y = as.factor(etoposide_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
etoposide_cv_mars$bestTune

#plot it
png(filename = 'Images/etoposide_cv_mars.png')
ggplot(etoposide_cv_mars) + ggtitle('Etoposide CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/etoposide_cv_mars_model.rds', etoposide_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(etoposide_cv_mars, num_features = 30, bar = FALSE, value = 'gcv') + ggtitle('etoposide GCV')

p2 <- vip(etoposide_cv_mars, num_features = 30, bar = FALSE, value = 'rss') + ggtitle('etoposide RSS')

png(filename = 'Images/etoposide_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()



## svm
etoposide_svm <- ksvm(res_sens ~ ., etoposide_rose)

summary(etoposide_svm) # is there better info than this?

etoposide_svm_predictions <- predict(etoposide_svm, etoposide_rna_seq_test)
mse <- mean((etoposide_test$res_sens - etoposide_svm_predictions) ^ 2)
mse #0.33


## kNN
etoposide_knn <- knnreg(etoposide_rose[ ,1:14209], etoposide_rose[, 14210], k = 2)

summary(etoposide_knn)
etoposide_knn_predictions <- predict(etoposide_knn, etoposide_rna_seq_test)
mse <- mean((etoposide_test$res_sens - etoposide_knn_predictions) ^ 2)
mse #0.22

# fit a basic MARS model
gemcitabine_mars1 <- earth(res_sens ~ ., data = gemcitabine_rose)
#print model summary
print(gemcitabine_mars1) #22
summary(gemcitabine_mars1)
coef(gemcitabine_mars1)

png(filename = 'Images/gemcitabine_first_degree_model_selection.png')
plot(gemcitabine_mars1, which = 1)
dev.off()


gemcitabine_mars1_predictions <- predict(gemcitabine_mars1, gemcitabine_rna_seq_test)
mse <- mean((gemcitabine_test$res_sens - gemcitabine_mars1_predictions) ^ 2) #0.71

# fit a basic MARS model
gemcitabine_mars2 <- earth(res_sens ~ ., data = gemcitabine_rose, degree = 2)

print(gemcitabine_mars2) #45
summary(gemcitabine_mars2)

## cv mars model
set.seed(5)
gemcitabine_cv_mars <- train(x = subset(gemcitabine_rose, select = -res_sens),
                           y = as.factor(gemcitabine_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
gemcitabine_cv_mars$bestTune

#plot it
png(filename = 'Images/gemcitabine_cv_mars.png')
ggplot(gemcitabine_cv_mars) + ggtitle('Gemcitabine CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/gemcitabine_cv_mars_model.rds', gemcitabine_cv_mars)

#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(gemcitabine_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('gemcitabine GCV')

p2 <- vip(gemcitabine_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('gemcitabine RSS')

png(filename = 'Images/gemcitabine_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


## svm
gemcitabine_svm <- ksvm(res_sens ~ ., gemcitabine_rose)

summary(gemcitabine_svm) # is there better info than this?

gemcitabine_svm_predictions <- predict(gemcitabine_svm, gemcitabine_rna_seq_test)
mse <- mean((gemcitabine_test$res_sens - gemcitabine_svm_predictions) ^ 2)
mse #0.46


## kNN
gemcitabine_knn <- knnreg(gemcitabine_rose[ ,1:14209], gemcitabine_rose[, 14210], k = 2)

summary(gemcitabine_knn)
gemcitabine_knn_predictions <- predict(gemcitabine_knn, gemcitabine_rna_seq_test)
mse <- mean((gemcitabine_test$res_sens - gemcitabine_knn_predictions) ^ 2)
mse #0.46

# fit a basic MARS model
methotrexate_mars1 <- earth(res_sens ~ ., data = methotrexate_rose)
#print model summary
print(methotrexate_mars1) #22
summary(methotrexate_mars1)
coef(methotrexate_mars1)

png(filename = 'Images/methotrexate_first_degree_model_selection.png')
plot(methotrexate_mars1, which = 1)
dev.off()


methotrexate_mars1_predictions <- predict(methotrexate_mars1, methotrexate_rna_seq_test)
mse <- mean((methotrexate_test$res_sens - methotrexate_mars1_predictions) ^ 2) #0.67

# fit a basic MARS model
methotrexate_mars2 <- earth(res_sens ~ ., data = methotrexate_rose, degree = 2)

print(methotrexate_mars2) #42
summary(methotrexate_mars2)

## cv mars model
set.seed(5)
methotrexate_cv_mars <- train(x = subset(methotrexate_rose, select = -res_sens),
                           y = as.factor(methotrexate_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
methotrexate_cv_mars$bestTune

#plot it
png(filename = 'Images/methotrexate_cv_mars.png')
ggplot(methotrexate_cv_mars) + ggtitle('Methotrexate CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/methotrexate_cv_mars_model.rds', methotrexate_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(methotrexate_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('methotrexate GCV')

p2 <- vip(methotrexate_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('methotrexate RSS')

png(filename = 'Images/methotrexate_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


## svm
methotrexate_svm <- ksvm(res_sens ~ ., methotrexate_rose)

summary(methotrexate_svm) # is there better info than this?

methotrexate_svm_predictions <- predict(methotrexate_svm, methotrexate_rna_seq_test)
mse <- mean((methotrexate_test$res_sens - methotrexate_svm_predictions) ^ 2)
mse #0.31


## kNN
methotrexate_knn <- knnreg(methotrexate_rose[ ,1:14209], methotrexate_rose[, 14210], k = 2)

summary(methotrexate_knn)
methotrexate_knn_predictions <- predict(methotrexate_knn, methotrexate_rna_seq_test)
mse <- mean((methotrexate_test$res_sens - methotrexate_knn_predictions) ^ 2)
mse #0.15

# fit a basic MARS model
mitomycin_mars1 <- earth(res_sens ~ ., data = mitomycin_rose)
#print model summary
print(mitomycin_mars1) #21
summary(mitomycin_mars1)
coef(mitomycin_mars1)

png(filename = 'Images/mitomycin_first_degree_model_selection.png')
plot(mitomycin_mars1, which = 1)
dev.off()


mitomycin_mars1_predictions <- predict(mitomycin_mars1, mitomycin_rna_seq_test)
mse <- mean((mitomycin_test$res_sens - mitomycin_mars1_predictions) ^ 2) #0.76

# fit a basic MARS model
mitomycin_mars2 <- earth(res_sens ~ ., data = mitomycin_rose, degree = 2)

print(mitomycin_mars2) #43
summary(mitomycin_mars2)

## cv mars model
set.seed(5)
mitomycin_cv_mars <- train(x = subset(mitomycin_rose, select = -res_sens),
                           y = as.factor(mitomycin_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
mitomycin_cv_mars$bestTune

#plot it
png(filename = 'Images/mitomycin_cv_mars.png')
ggplot(mitomycin_cv_mars) + ggtitle('Mitomycin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/mitomycin_cv_mars_model.rds', mitomycin_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(mitomycin_cv_mars, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('mitomycin GCV')

p2 <- vip(mitomycin_cv_mars, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('mitomycin RSS')

png(filename = 'Images/mitomycin_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


## svm
mitomycin_svm <- ksvm(res_sens ~ ., mitomycin_rose)

summary(mitomycin_svm) # is there better info than this?

mitomycin_svm_predictions <- predict(mitomycin_svm, mitomycin_rna_seq_test)
mse <- mean((mitomycin_test$res_sens - mitomycin_svm_predictions) ^ 2)
mse #0.45


## kNN
mitomycin_knn <- knnreg(mitomycin_rose[ ,1:14209], mitomycin_rose[, 14210], k = 2)

summary(mitomycin_knn)
mitomycin_knn_predictions <- predict(mitomycin_knn, mitomycin_rna_seq_test)
mse <- mean((mitomycin_test$res_sens - mitomycin_knn_predictions) ^ 2)
mse #0.30

# fit a basic MARS model
sn38_mars1 <- earth(res_sens ~ ., data = sn38_rose)
#print model summary
print(sn38_mars1) #24
summary(sn38_mars1)
coef(sn38_mars1)

png(filename = 'Images/sn38_first_degree_model_selection.png')
plot(sn38_mars1, which = 1)
dev.off()


sn38_mars1_predictions <- predict(sn38_mars1, sn38_rna_seq_test)
mse <- mean((sn38_test$res_sens - sn38_mars1_predictions) ^ 2) #0.63

# fit a basic MARS model
sn38_mars2 <- earth(res_sens ~ ., data = sn38_rose, degree = 2)

print(sn38_mars2) #49
summary(sn38_mars2)

## cv mars model
set.seed(5)
sn38_cv_mars <- train(x = subset(sn38_rose, select = -res_sens),
                           y = as.factor(sn38_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
sn38_cv_mars$bestTune

#plot it
png(filename = 'Images/sn38_cv_mars.png')
ggplot(sn38_cv_mars) + ggtitle('Sn38 CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/sn38_cv_mars_model.rds', sn38_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(sn38_cv_mars, num_features = 30, bar = FALSE, value = 'gcv') + ggtitle('sn38 GCV')

p2 <- vip(sn38_cv_mars, num_features = 30, bar = FALSE, value = 'rss') + ggtitle('sn38 RSS')

png(filename = 'Images/sn38_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
sn38_svm <- ksvm(res_sens ~ ., sn38_rose)

summary(sn38_svm) # is there better info than this?

sn38_svm_predictions <- predict(sn38_svm, sn38_rna_seq_test)
mse <- mean((sn38_test$res_sens - sn38_svm_predictions) ^ 2)
mse #0.35


## kNN
sn38_knn <- knnreg(sn38_rose[ ,1:14209], sn38_rose[, 14210], k = 2)

summary(sn38_knn)
sn38_knn_predictions <- predict(sn38_knn, sn38_rna_seq_test)
mse <- mean((sn38_test$res_sens - sn38_knn_predictions) ^ 2)
mse #0.36

# fit a basic MARS model
temozolomide_mars1 <- earth(res_sens ~ ., data = temozolomide_rose)
#print model summary
print(temozolomide_mars1) #25
summary(temozolomide_mars1)
coef(temozolomide_mars1)

png(filename = 'Images/temozolomide_first_degree_model_selection.png')
plot(temozolomide_mars1, which = 1)
dev.off()


temozolomide_mars1_predictions <- predict(temozolomide_mars1, temozolomide_rna_seq_test)
mse <- mean((temozolomide_test$res_sens - temozolomide_mars1_predictions) ^ 2) #0.74

# fit a basic MARS model
temozolomide_mars2 <- earth(res_sens ~ ., data = temozolomide_rose, degree = 2)

print(temozolomide_mars2) #43
summary(temozolomide_mars2)

## cv mars model
set.seed(5)
temozolomide_cv_mars <- train(x = subset(temozolomide_rose, select = -res_sens),
                           y = as.factor(temozolomide_rose$res_sens),
                           method = 'earth',
                           trControl = trainControl(method = 'cv', number = 10),
                           tuneGrid = hyper_grid)

#results
temozolomide_cv_mars$bestTune

#plot it
png(filename = 'Images/temozolomide_cv_mars.png')
ggplot(temozolomide_cv_mars) + ggtitle('Temozolomide CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/temozolomide_cv_mars_model.rds', temozolomide_cv_mars)


#refine grid search (nprune)

## feature importance
# this should be done on cv_mars
p1 <- vip(temozolomide_cv_mars, num_features = 30, bar = FALSE, value = 'gcv') + ggtitle('temozolomide GCV')

p2 <- vip(temozolomide_cv_mars, num_features = 30, bar = FALSE, value = 'rss') + ggtitle('temozolomide RSS')

png(filename = 'Images/temozolomide_feat_imp.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## svm
temozolomide_svm <- ksvm(res_sens ~ ., temozolomide_rose)

summary(temozolomide_svm) # is there better info than this?

temozolomide_svm_predictions <- predict(temozolomide_svm, temozolomide_rna_seq_test)
mse <- mean((temozolomide_test$res_sens - temozolomide_svm_predictions) ^ 2)
mse #0.35


## kNN
temozolomide_knn <- knnreg(temozolomide_rose[ ,1:14209], temozolomide_rose[, 14210], k = 2)

summary(temozolomide_knn)
temozolomide_knn_predictions <- predict(temozolomide_knn, temozolomide_rna_seq_test)
mse <- mean((temozolomide_test$res_sens - temozolomide_knn_predictions) ^ 2)
mse #0.12





##refined grid searches
## bleomycin
hyper_grid <- expand.grid(degree = 1, nprune = seq(29,35,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
bleomycin_cv_mars_refined <- train(x = subset(bleomycin_rose, select = -res_sens), 
                           y = as.factor(bleomycin_rose$res_sens), 
                           method = 'earth', 
                           trControl = trainControl(method = 'cv', number = 10), 
                           tuneGrid = hyper_grid)

#results
bleomycin_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/bleomycin_cv_mars_refined.png')
ggplot(bleomycin_cv_mars_refined) + ggtitle('Bleomycin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/bleomycin_cv_mars_refined_model.rds', bleomycin_cv_mars_refined)

p1 <- vip(bleomycin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('Bleomycin GCV')
# 
p2 <- vip(bleomycin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('Bleomycin RSS')
# 
png(filename = 'Images/bleomycin_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## camptothecin
hyper_grid <- expand.grid(degree = 1:3, nprune = seq(14,20,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
camptothecin_cv_mars_refined <- train(x = subset(camptothecin_rose, select = -res_sens), 
                                   y = as.factor(camptothecin_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
camptothecin_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/camptothecin_cv_mars_refined.png')
ggplot(camptothecin_cv_mars_refined) + ggtitle('Camptothecin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/camptothecin_cv_mars_refined_model.rds', camptothecin_cv_mars_refined)

p1 <- vip(camptothecin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('Camptothecin GCV')
# 
p2 <- vip(camptothecin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('Camptothecin RSS')
# 
png(filename = 'Images/camptothecin_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## cisplatin
hyper_grid <- expand.grid(degree = 1, nprune = seq(21,27,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
cisplatin_cv_mars_refined <- train(x = subset(cisplatin_rose, select = -res_sens), 
                                   y = as.factor(cisplatin_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
cisplatin_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/cisplatin_cv_mars_refined.png')
ggplot(cisplatin_cv_mars_refined) + ggtitle('cisplatin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/cisplatin_cv_mars_refined_model.rds', cisplatin_cv_mars_refined)

p1 <- vip(cisplatin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('cisplatin GCV')
# 
p2 <- vip(cisplatin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('cisplatin RSS')
# 
png(filename = 'Images/cisplatin_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## cytarabine
hyper_grid <- expand.grid(degree = 1, nprune = seq(14,20,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
cytarabine_cv_mars_refined <- train(x = subset(cytarabine_rose, select = -res_sens), 
                                   y = as.factor(cytarabine_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
cytarabine_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/cytarabine_cv_mars_refined.png')
ggplot(cytarabine_cv_mars_refined) + ggtitle('cytarabine CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/cytarabine_cv_mars_refined_model.rds', cytarabine_cv_mars_refined)

p1 <- vip(cytarabine_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('cytarabine GCV')
# 
p2 <- vip(cytarabine_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('cytarabine RSS')
# 
png(filename = 'Images/cytarabine_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## doxorubicin
hyper_grid <- expand.grid(degree = 1, nprune = seq(21,27,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
doxorubicin_cv_mars_refined <- train(x = subset(doxorubicin_rose, select = -res_sens), 
                                   y = as.factor(doxorubicin_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
doxorubicin_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/doxorubicin_cv_mars_refined.png')
ggplot(doxorubicin_cv_mars_refined) + ggtitle('doxorubicin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/doxorubicin_cv_mars_refined_model.rds', doxorubicin_cv_mars_refined)

p1 <- vip(doxorubicin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('doxorubicin GCV')
# 
p2 <- vip(doxorubicin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('doxorubicin RSS')
# 
png(filename = 'Images/doxorubicin_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## etoposide
hyper_grid <- expand.grid(degree = 1, nprune = seq(21,27,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
etoposide_cv_mars_refined <- train(x = subset(etoposide_rose, select = -res_sens), 
                                   y = as.factor(etoposide_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
etoposide_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/etoposide_cv_mars_refined.png')
ggplot(etoposide_cv_mars_refined) + ggtitle('etoposide CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/etoposide_cv_mars_refined_model.rds', etoposide_cv_mars_refined)

p1 <- vip(etoposide_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('etoposide GCV')
# 
p2 <- vip(etoposide_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('etoposide RSS')
# 
png(filename = 'Images/etoposide_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## gemcitabine
hyper_grid <- expand.grid(degree = 1, nprune = seq(29,35,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
gemcitabine_cv_mars_refined <- train(x = subset(gemcitabine_rose, select = -res_sens), 
                                   y = as.factor(gemcitabine_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
gemcitabine_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/gemcitabine_cv_mars_refined.png')
ggplot(gemcitabine_cv_mars_refined) + ggtitle('gemcitabine CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/gemcitabine_cv_mars_refined_model.rds', gemcitabine_cv_mars_refined)

p1 <- vip(gemcitabine_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('gemcitabine GCV')
# 
p2 <- vip(gemcitabine_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('gemcitabine RSS')
# 
png(filename = 'Images/gemcitabine_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## methotrexate
hyper_grid <- expand.grid(degree = 1:3, nprune = seq(14,20,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
methotrexate_cv_mars_refined <- train(x = subset(methotrexate_rose, select = -res_sens), 
                                   y = as.factor(methotrexate_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
methotrexate_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/methotrexate_cv_mars_refined.png')
ggplot(methotrexate_cv_mars_refined) + ggtitle('methotrexate CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/methotrexate_cv_mars_refined_model.rds', methotrexate_cv_mars_refined)

p1 <- vip(methotrexate_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('methotrexate GCV')
# 
p2 <- vip(methotrexate_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('methotrexate RSS')
# 
png(filename = 'Images/methotrexate_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## mitomycin
hyper_grid <- expand.grid(degree = 1, nprune = seq(21,27,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
mitomycin_cv_mars_refined <- train(x = subset(mitomycin_rose, select = -res_sens), 
                                   y = as.factor(mitomycin_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
mitomycin_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/mitomycin_cv_mars_refined.png')
ggplot(mitomycin_cv_mars_refined) + ggtitle('mitomycin CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/mitomycin_cv_mars_refined_model.rds', mitomycin_cv_mars_refined)

p1 <- vip(mitomycin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('mitomycin GCV')
# 
p2 <- vip(mitomycin_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('mitomycin RSS')
# 
png(filename = 'Images/mitomycin_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## sn38
hyper_grid <- expand.grid(degree = 1, nprune = seq(21,27,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
sn38_cv_mars_refined <- train(x = subset(sn38_rose, select = -res_sens), 
                                   y = as.factor(sn38_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
sn38_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/sn38_cv_mars_refined.png')
ggplot(sn38_cv_mars_refined) + ggtitle('sn38 CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/sn38_cv_mars_refined_model.rds', sn38_cv_mars_refined)

p1 <- vip(sn38_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('sn38 GCV')
# 
p2 <- vip(sn38_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('sn38 RSS')
# 
png(filename = 'Images/sn38_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

## temozolomide
hyper_grid <- expand.grid(degree = 1, nprune = seq(29,35,length.out = 7) %>% floor())
hyper_grid

set.seed(5)
temozolomide_cv_mars_refined <- train(x = subset(temozolomide_rose, select = -res_sens), 
                                   y = as.factor(temozolomide_rose$res_sens), 
                                   method = 'earth', 
                                   trControl = trainControl(method = 'cv', number = 10), 
                                   tuneGrid = hyper_grid)

#results
temozolomide_cv_mars_refined$bestTune

#plot it
png(filename = 'Images/temozolomide_cv_mars_refined.png')
ggplot(temozolomide_cv_mars_refined) + ggtitle('temozolomide CV Mars')
dev.off()

saveRDS(file = 'GLM_Models/temozolomide_cv_mars_refined_model.rds', temozolomide_cv_mars_refined)

p1 <- vip(temozolomide_cv_mars_refined, num_features = 23, bar = FALSE, value = 'gcv') + ggtitle('temozolomide GCV')
# 
p2 <- vip(temozolomide_cv_mars_refined, num_features = 23, bar = FALSE, value = 'rss') + ggtitle('temozolomide RSS')
# 
png(filename = 'Images/temozolomide_feat_imp_refined.png')
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()



## testing for accuracy
bleomycin_test_cv_mars <- predict(bleomycin_mars1, newdata = as.matrix(bleomycin_rna_seq), type = 'class') #class for everything else

bleomycin_test_cv_mars_auc <- auc(bleomycin_test$res_sens, bleomycin_test_cv_mars)
bleomycin_test_cv_mars_auc <- round(bleomycin_test_cv_mars_auc, digits = 2)
bleomycin_cv_mars_acc <- sum(bleomycin$res_sens == bleomycin_test_cv_mars)/length(bleomycin_test_cv_mars)

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

new_blca_tcga_cisplatin <- predict(cisplatin_mars1, newdata = as.matrix(blca_gene_short_scaled), type = 'class', na.action = 'na.pass')

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
legend(x = 450, y = 0.4, legend = c('predicted sensitive', 'predicted resistant'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)

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

new_blca_tcga_gemcitabine <- predict(gemcitabine_mars1, newdata = as.matrix(blca_gene_short_scaled), type = 'class')

blca_surv_times <- blca_clinical_gemcitabine_short$PFS
blca_status <- ifelse(blca_clinical_gemcitabine_short$PFS == blca_clinical_gemcitabine_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_gemcitabine)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
                data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('darkviolet', 'limegreen'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2)
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted sensitive', 'predicted resistant'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)






# ### CCLE -----
# ## load clinical data ----
# carboplatin_ccle <- read.csv('Processed_Clinical_Data/carboplatin_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
# cyclophosphamide_ccle <- read.csv('Processed_Clinical_Data/cyclophosphamide_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
# docetaxel_ccle <- read.csv('Processed_Clinical_Data/docetaxel_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
# fluorouracil_ccle <- read.csv('Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
# gemcitabine_ccle <- read.csv('Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
# 
# ## load gene expression data -----
# ccle_microarray <- read.csv('Processed_Gene_Expression/ccle_microarray_processed.csv', row.names = 1)
# colnames(ccle_microarray) <- gsub('_.*$', '', colnames(ccle_microarray))
# colnames(ccle_microarray) <- gsub('^X', '', colnames(ccle_microarray))
# 
# ## set up data for model building ----
# carboplatin_lines       <- carboplatin_ccle$Cell.line.name #650
# carboplatin_lines       <- as.character(carboplatin_lines)
# cyclophosphamide_lines  <- cyclophosphamide_ccle$Cell.line.name #335
# cyclophosphamide_lines  <- as.character(cyclophosphamide_lines)
# docetaxel_lines         <- docetaxel_ccle$Cell.line.name #302
# docetaxel_lines         <- as.character(docetaxel_lines)
# fluorouracil_lines      <- fluorouracil_ccle$Cell.line.name #652
# fluorouracil_lines      <- as.character(fluorouracil_lines)
# gemcitabine_lines       <- gemcitabine_ccle$Cell.line.name #589
# gemcitabine_lines       <- as.character(gemcitabine_lines)
# 
# # set ccle dataframe correct orientation
# ccle <- data.frame(t(ccle_microarray))
# #855 x 14209
# 
# #split into training/testing
# set.seed(5)
# 
# random_sample <- sample(x = rownames(ccle), size = nrow(ccle)/5)
# 
# 
# ccle_train         <- ccle[which(rownames(ccle) %ni% random_sample), ] #427 x 14209
# 
# ccle_test          <- ccle[random_sample, ] #481 x 14209
# 
# intersect(rownames(ccle_train), rownames(ccle_test)) #0
# 
# #split expression data
# carboplatin_rna_seq_train           <- ccle_train[intersect(carboplatin_lines, rownames(ccle_train)), ]
# # 312 x 14209
# carboplatin_rna_seq_test            <- ccle_test[intersect(carboplatin_lines, rownames(ccle_test)), ]
# # 306 x 14209
# 
# cyclophosphamide_rna_seq_train      <- ccle_train[intersect(cyclophosphamide_lines, rownames(ccle_train)), ]
# # 165 x 14209
# cyclophosphamide_rna_seq_test       <- ccle_test[intersect(cyclophosphamide_lines, rownames(ccle_test)), ]
# # 151 x 14209
# 
# docetaxel_rna_seq_train             <- ccle_train[intersect(docetaxel_lines, rownames(ccle_train)), ]
# # 142 x 14209
# docetaxel_rna_seq_test              <- ccle_test[intersect(docetaxel_lines, rownames(ccle_test)), ]
# # 155 x 14209
# 
# fluorouracil_rna_seq_train          <- ccle_train[intersect(fluorouracil_lines, rownames(ccle_train)), ]
# # 313 x 14209
# fluorouracil_rna_seq_test           <- ccle_test[intersect(fluorouracil_lines, rownames(ccle_test)), ]
# # 307 x 14209
# 
# gemcitabine_rna_seq_train           <- ccle_train[intersect(gemcitabine_lines, rownames(ccle_train)), ]
# # 276 x 14209
# gemcitabine_rna_seq_test            <- ccle_test[intersect(gemcitabine_lines, rownames(ccle_test)), ]
# # 286 x 14209
# 
# #split clinical data
# carboplatin_train         <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_train)), ]
# carboplatin_test          <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_test)), ]
# 
# cyclophosphamide_train    <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_train)), ]
# cyclophosphamide_test     <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_test)), ]
# 
# docetaxel_train           <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_train)), ]
# docetaxel_test            <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_test)), ]
# 
# fluorouracil_train        <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_train)), ]
# fluorouracil_test         <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_test)), ]
# 
# gemcitabine_train         <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_train)), ]
# gemcitabine_test          <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_test)), ]
# 
# #scale data
# carboplatin_rna_seq_train_scaled          <- apply(carboplatin_rna_seq_train, 2, scale)
# carboplatin_rna_seq_test_scaled           <- as.data.frame(apply(carboplatin_rna_seq_test, 2, scale))
# 
# cyclophosphamide_rna_seq_train_scaled     <- apply(cyclophosphamide_rna_seq_train, 2, scale)
# cyclophosphamide_rna_seq_test_scaled      <- as.data.frame(apply(cyclophosphamide_rna_seq_test, 2, scale))
# 
# docetaxel_rna_seq_train_scaled            <- apply(docetaxel_rna_seq_train, 2, scale)
# docetaxel_rna_seq_test_scaled             <- as.data.frame(apply(docetaxel_rna_seq_test, 2, scale))
# 
# fluorouracil_rna_seq_train_scaled         <- apply(fluorouracil_rna_seq_train, 2, scale)
# fluorouracil_rna_seq_test_scaled          <- as.data.frame(apply(fluorouracil_rna_seq_test, 2, scale))
# 
# gemcitabine_rna_seq_train_scaled          <- apply(gemcitabine_rna_seq_train, 2, scale)
# gemcitabine_rna_seq_test_scaled           <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))
# 
# 
# ### build models ----
# ## carboplatin
# carboplatin_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(carboplatin_rna_seq_train_scaled), y = carboplatin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# #save plot
# png(filename = 'Images/carboplatin_ccle_most_auc.png')
# plot(carboplatin_ccle_most_fit_elnet)
# dev.off()
# #save model
# saveRDS(file = 'GLM_Models/carboplatin_ccle_most_model.rds', carboplatin_ccle_most_fit_elnet)
# 
# carboplatin_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(carboplatin_rna_seq_train_scaled), y = carboplatin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/carboplatin_ccle_least_auc.png')
# plot(carboplatin_ccle_least_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/carboplatin_ccle_least_model.rds', carboplatin_ccle_least_fit_elnet)
# 
# ## cyclophosphamide
# cyclophosphamide_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(cyclophosphamide_rna_seq_train_scaled), y = cyclophosphamide_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/cyclophosphamide_ccle_most_auc.png')
# plot(cyclophosphamide_ccle_most_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/cyclophosphamide_ccle_most_model.rds', cyclophosphamide_ccle_most_fit_elnet)
# 
# cyclophosphamide_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(cyclophosphamide_rna_seq_train_scaled), y = cyclophosphamide_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/cyclophosphamide_ccle_least_auc.png')
# plot(cyclophosphamide_ccle_least_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/cyclophosphamide_ccle_least_model.rds', cyclophosphamide_ccle_least_fit_elnet)
# 
# ## docetaxel
# docetaxel_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(docetaxel_rna_seq_train_scaled), y = docetaxel_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/docetaxel_ccle_most_auc.png')
# plot(docetaxel_ccle_most_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/docetaxel_ccle_most_model.rds', docetaxel_ccle_most_fit_elnet)
# 
# docetaxel_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(docetaxel_rna_seq_train_scaled), y = docetaxel_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/docetaxel_ccle_least_auc.png')
# plot(docetaxel_ccle_least_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/docetaxel_ccle_least_model.rds', docetaxel_ccle_least_fit_elnet)
# 
# ## fluorouracil
# fluorouracil_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(fluorouracil_rna_seq_train_scaled), y = fluorouracil_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/fluorouracil_ccle_most_auc.png')
# plot(fluorouracil_ccle_most_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/fluorouracil_ccle_most_model.rds', fluorouracil_ccle_most_fit_elnet)
# 
# fluorouracil_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(fluorouracil_rna_seq_train_scaled), y = fluorouracil_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/fluorouracil_ccle_least_auc.png')
# plot(fluorouracil_ccle_least_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/fluorouracil_ccle_least_model.rds', fluorouracil_ccle_least_fit_elnet)
# 
# ## gemcitabine
# gemcitabine_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/gemcitabine_ccle_most_auc.png')
# plot(gemcitabine_ccle_most_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/gemcitabine_ccle_most_model.rds', gemcitabine_ccle_most_fit_elnet)
# 
# gemcitabine_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# 
# png(filename = 'Images/gemcitabine_ccle_least_auc.png')
# plot(gemcitabine_ccle_least_fit_elnet)
# dev.off()
# 
# saveRDS(file = 'GLM_Models/gemcitabine_ccle_least_model.rds', gemcitabine_ccle_least_fit_elnet)