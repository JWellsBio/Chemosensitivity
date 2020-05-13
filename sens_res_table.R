# create table showing resistant/sensitive numbers, thus the need for oversampling

## load necessary packages ----
if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) # for model building

if (!require ('ROSE')) install.packages('ROSE')
library(ROSE)

if (!require('vip')) install.packages('vip')
library(vip)

if (!require('pdp')) install.pacakges('pdp') 
library(pdp)

if(!require('patchwork')) install.packages('patchwork')
library(patchwork)

if(!require('pROC')) install.packages('pROC')
library(pROC)


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

bleomycin_rna_seq <- gdsc[rownames(gdsc) %in% bleomycin$COSMIC_ID, ]
bleomycin_rna_seq <- as.data.frame(bleomycin_rna_seq)

bleomycin_rna_seq$res_sens <- bleomycin$res_sens
bleomycin_table <- table(bleomycin_rna_seq$res_sens)

camptothecin_rna_seq <- gdsc[rownames(gdsc) %in% camptothecin$COSMIC_ID, ]
camptothecin_rna_seq <- as.data.frame(camptothecin_rna_seq)

camptothecin_rna_seq$res_sens <- camptothecin$res_sens
camptothecin_table <- table(camptothecin_rna_seq$res_sens)

cisplatin_rna_seq <- gdsc[rownames(gdsc) %in% cisplatin$COSMIC_ID, ]
cisplatin_rna_seq <- as.data.frame(cisplatin_rna_seq)

cisplatin_rna_seq$res_sens <- cisplatin$res_sens
cisplatin_table <- table(cisplatin_rna_seq$res_sens)

cytarabine_rna_seq <- gdsc[rownames(gdsc) %in% cytarabine$COSMIC_ID, ]
cytarabine_rna_seq <- as.data.frame(cytarabine_rna_seq)

cytarabine_rna_seq$res_sens <- cytarabine$res_sens
cytarabine_table <- table(cytarabine_rna_seq$res_sens)

doxorubicin_rna_seq <- gdsc[rownames(gdsc) %in% doxorubicin$COSMIC_ID, ]
doxorubicin_rna_seq <- as.data.frame(doxorubicin_rna_seq)

doxorubicin_rna_seq$res_sens <- doxorubicin$res_sens
doxorubicin_table <- table(doxorubicin_rna_seq$res_sens)

etoposide_rna_seq <- gdsc[rownames(gdsc) %in% etoposide$COSMIC_ID, ]
etoposide_rna_seq <- as.data.frame(etoposide_rna_seq)

etoposide_rna_seq$res_sens <- etoposide$res_sens
etoposide_table <- table(etoposide_rna_seq$res_sens)

gemcitabine_rna_seq <- gdsc[rownames(gdsc) %in% gemcitabine$COSMIC_ID, ]
gemcitabine_rna_seq <- as.data.frame(gemcitabine_rna_seq)

gemcitabine_rna_seq$res_sens <- gemcitabine$res_sens
gemcitabine_table <- table(gemcitabine_rna_seq$res_sens)

methotrexate_rna_seq <- gdsc[rownames(gdsc) %in% methotrexate$COSMIC_ID, ]
methotrexate_rna_seq <- as.data.frame(methotrexate_rna_seq)

methotrexate_rna_seq$res_sens <- methotrexate$res_sens
methotrexate_table <- table(methotrexate_rna_seq$res_sens)

mitomycin_rna_seq <- gdsc[rownames(gdsc) %in% mitomycin$COSMIC_ID, ]
mitomycin_rna_seq <- as.data.frame(mitomycin_rna_seq)

mitomycin_rna_seq$res_sens <- mitomycin$res_sens
mitomycin_table <- table(mitomycin_rna_seq$res_sens)

sn38_rna_seq <- gdsc[rownames(gdsc) %in% sn38$COSMIC_ID, ]
sn38_rna_seq <- as.data.frame(sn38_rna_seq)

sn38_rna_seq$res_sens <- sn38$res_sens
sn38_table <- table(sn38_rna_seq$res_sens)

temozolomide_rna_seq <- gdsc[rownames(gdsc) %in% temozolomide$COSMIC_ID, ]
temozolomide_rna_seq <- as.data.frame(temozolomide_rna_seq)

temozolomide_rna_seq$res_sens <- temozolomide$res_sens
temozolomide_table <- table(temozolomide_rna_seq$res_sens)

sensitivity_thresholds <- c(-1.4805, -6.584, 1.3801, -1.9516, -3.9565, -1.2198, -5.9903, -2.4743, -2.9647, -6.559, 4.6032)

sensitive_samples <- c(bleomycin_table[1], camptothecin_table[1], cisplatin_table[1], cytarabine_table[1], doxorubicin_table[1], 
                       etoposide_table[1], gemcitabine_table[1], methotrexate_table[1], mitomycin_table[1], sn38_table[1], 
                       temozolomide_table[1])

resistant_sampels <- c(bleomycin_table[2], camptothecin_table[2], cisplatin_table[2], cytarabine_table[2], doxorubicin_table[2], 
                       etoposide_table[2], gemcitabine_table[2], methotrexate_table[2], mitomycin_table[2], sn38_table[2], 
                       temozolomide_table[2])

sample_df <- data.frame(sensitivity_thresholds, sensitive_samples, resistant_sampels)
colnames(sample_df) <- c('Sensitivity\nThreshold', 'Sensitive\nSamples', 'Resistant\nSamples')
rownames(sample_df) <- c('Bleomycin', 'Camptothecin', 'Cisplatin', 'Cytarabine', 'Doxorubicin', 
                         'Etoposide', 'Gemcitabine', 'Methotrexate', 'Mitomycin', 'SN38', 
                         'Temozolomide')

