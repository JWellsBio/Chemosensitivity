# all methotrexate models

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



methotrexate     <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

methotrexate_lines           <- methotrexate$COSMIC_ID #680


methotrexate_rna_seq     <- gdsc[intersect(methotrexate_lines, rownames(gdsc)), ]


methotrexate_rna_seq <- as.data.frame(methotrexate_rna_seq)
methotrexate$med_sens <- ifelse(methotrexate$LN_IC50 <= median(methotrexate$LN_IC50), 1, 0)
methotrexate_rna_seq$med_sens <- methotrexate$med_sens

set.seed(5)
train_index <- createDataPartition(methotrexate_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

methotrexate_rna_seq_train <- methotrexate_rna_seq[ train_index,]
methotrexate_rna_seq_test  <- methotrexate_rna_seq[-train_index,]

methotrexate_elastic      <- readRDS('GLM_Models/methotrexate_glm_alpha_search_model.rds')
methotrexate_lasso <- readRDS('GLM_Models/methotrexate_glm_lasso_model.rds')
methotrexate_med_lasso <- readRDS('GLM_Models/methotrexate_glm_med_lasso_model.rds')
methotrexate_med_elastic <- readRDS('GLM_Models/methotrexate_med_glm_alpha_search_model.rds')


# BRCA W METHOTREXATE (9)
brca_clinical <- read.csv('Processed_Clinical_Data/brca_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(brca_clinical$most_sensitive)
brca_clinical <- brca_clinical[!na_idx, ]
table(brca_clinical$drug_name)
brca_clinical_methotrexate <- brca_clinical[which(brca_clinical$drug_name == 'Methotrexate' | brca_clinical$drug_name == 'METHOTREXATE'), ]


brca_gene <- read.csv('Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv', row.names = 1)
colnames(brca_gene) <- gsub('\\.', '-', colnames(brca_gene))
brca_matching_idx <- brca_clinical_methotrexate$submitter_id.samples %in% colnames(brca_gene)
brca_clinical_methotrexate_short <- brca_clinical_methotrexate[brca_matching_idx, ]
brca_matching_idx <- colnames(brca_gene) %in% brca_clinical_methotrexate_short$submitter_id.samples
brca_gene_short <- brca_gene[, brca_matching_idx]
brca_gene_short <- t(brca_gene_short)
brca_gene_short_scaled <- apply(brca_gene_short, 2, scale)
brca_gene_short_scaled <- fix_scale(brca_gene_short_scaled)


## elastic
new_brca_tcga_methotrexate <- predict(methotrexate_elastic, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

brca_surv_times <- brca_clinical_methotrexate_short$PFS
brca_status <- ifelse(brca_clinical_methotrexate_short$PFS == brca_clinical_methotrexate_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_methotrexate)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_methotrexate,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_methotrexate,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
#all called resistant
