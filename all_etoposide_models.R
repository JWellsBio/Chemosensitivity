#all etoposide models

## make this a function

fix_scale <- function(gene_scaled) {
  var_check <- rep(NA, ncol(gene_scaled))
  var_seq <- seq(-1, 1, 2/nrow(gene_scaled))
  var_seq <- var_seq[-1]
  var_seq <- replicate(nrow(gene_scaled), sample(var_seq))
  for (i in 1:length(var_check)) {
    if (is.na(var(gene_scaled[, i]))) {
      gene_scaled[, i] <- sample(var_seq)[1:nrow(gene_scaled)]
    }
  }
  return(gene_scaled)
}
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
library(caret)



etoposide     <- read.csv('Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

etoposide_lines           <- etoposide$COSMIC_ID #680


etoposide_rna_seq     <- gdsc[intersect(etoposide_lines, rownames(gdsc)), ]


etoposide_rna_seq <- as.data.frame(etoposide_rna_seq)
etoposide$med_sens <- ifelse(etoposide$LN_IC50 <= median(etoposide$LN_IC50), 1, 0)
etoposide_rna_seq$med_sens <- etoposide$med_sens

set.seed(5)
train_index <- createDataPartition(etoposide_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

etoposide_rna_seq_train <- etoposide_rna_seq[ train_index,]
etoposide_rna_seq_test  <- etoposide_rna_seq[-train_index,]

etoposide_elastic      <- readRDS('GLM_Models/etoposide_glm_alpha_search_model.rds')
etoposide_lasso <- readRDS('GLM_Models/etoposide_glm_lasso_model.rds')
etoposide_med_lasso <- readRDS('GLM_Models/etoposide_glm_med_lasso_model.rds')
etoposide_med_elastic <- readRDS('GLM_Models/etoposide_med_glm_alpha_search_model.rds')


# LUAD W ETOPOSIDE (7)
luad_clinical <- read.csv('Processed_Clinical_Data/luad_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(luad_clinical$most_sensitive)
luad_clinical <- luad_clinical[!na_idx, ]
table(luad_clinical$drug_name)
luad_clinical_etoposide <- luad_clinical[which(luad_clinical$drug_name == 'Etoposide'), ]


luad_gene <- read.csv('Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv', row.names = 1)
colnames(luad_gene) <- gsub('\\.', '-', colnames(luad_gene))
luad_matching_idx <- luad_clinical_etoposide$submitter_id.samples %in% colnames(luad_gene)
luad_clinical_etoposide_short <- luad_clinical_etoposide[luad_matching_idx, ]
luad_matching_idx <- colnames(luad_gene) %in% luad_clinical_etoposide_short$submitter_id.samples
luad_gene_short <- luad_gene[, luad_matching_idx]
luad_gene_short <- t(luad_gene_short)
luad_gene_short_scaled <- apply(luad_gene_short, 2, scale)
luad_gene_short_scaled <- fix_scale(luad_gene_short_scaled)


## elastic
new_luad_tcga_etoposide <- predict(etoposide_elastic, newx = as.matrix(luad_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

luad_surv_times <- luad_clinical_etoposide_short$PFS
luad_status <- ifelse(luad_clinical_etoposide_short$PFS == luad_clinical_etoposide_short$OS, 0, 1)

luad_surv_df <- data.frame(luad_surv_times, luad_status, new_luad_tcga_etoposide)
fit <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
               data = luad_surv_df)
fit2 <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
                data = luad_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(luad_surv_df, file = 'Survival_Data/luad_etoposide_surv_df.txt', row.names = FALSE)



# LUSC W ETOPOSIDE (6)
lusc_clinical <- read.csv('Processed_Clinical_Data/lusc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(lusc_clinical$most_sensitive)
lusc_clinical <- lusc_clinical[!na_idx, ]
table(lusc_clinical$drug_name)
lusc_clinical_etoposide <- lusc_clinical[which(lusc_clinical$drug_name == 'etoposide' | lusc_clinical$drug_name == 'Etoposide'), ]

lusc_gene <- read.csv('Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(lusc_gene) <- gsub('\\.', '-', colnames(lusc_gene))
lusc_matching_idx <- lusc_clinical_etoposide$submitter_id.samples %in% colnames(lusc_gene)
lusc_clinical_etoposide_short <- lusc_clinical_etoposide[lusc_matching_idx, ]
lusc_matching_idx <- colnames(lusc_gene) %in% lusc_clinical_etoposide_short$submitter_id.samples
lusc_gene_short <- lusc_gene[, lusc_matching_idx]
lusc_gene_short <- t(lusc_gene_short)
lusc_gene_short_scaled <- apply(lusc_gene_short, 2, scale)
lusc_gene_short_scaled <- fix_scale(lusc_gene_short_scaled)


## elastic
new_lusc_tcga_etoposide <- predict(etoposide_elastic, newx = as.matrix(lusc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

lusc_surv_times <- lusc_clinical_etoposide_short$PFS
lusc_status <- ifelse(lusc_clinical_etoposide_short$PFS == lusc_clinical_etoposide_short$OS, 0, 1)

lusc_surv_df <- data.frame(lusc_surv_times, lusc_status, new_lusc_tcga_etoposide)
fit <- survfit(Surv(lusc_surv_times, lusc_status) ~ X1,
               data = lusc_surv_df)
fit2 <- survfit(Surv(lusc_surv_times, lusc_status) ~ X1,
                data = lusc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(lusc_surv_df, file = 'Survival_Data/lusc_etoposide_surv_df.txt', row.names = FALSE)
