#all temozolomide models

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



temozolomide     <- read.csv('Processed_Clinical_Data/temozolomide_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

temozolomide_lines           <- temozolomide$COSMIC_ID #680


temozolomide_rna_seq     <- gdsc[intersect(temozolomide_lines, rownames(gdsc)), ]


temozolomide_rna_seq <- as.data.frame(temozolomide_rna_seq)
temozolomide$med_sens <- ifelse(temozolomide$LN_IC50 <= median(temozolomide$LN_IC50), 1, 0)
temozolomide_rna_seq$med_sens <- temozolomide$med_sens

set.seed(5)
train_index <- createDataPartition(temozolomide_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

temozolomide_rna_seq_train <- temozolomide_rna_seq[ train_index,]
temozolomide_rna_seq_test  <- temozolomide_rna_seq[-train_index,]

temozolomide_elastic      <- readRDS('GLM_Models/temozolomide_glm_alpha_search_model.rds')
temozolomide_lasso <- readRDS('GLM_Models/temozolomide_glm_lasso_model.rds')
temozolomide_med_lasso <- readRDS('GLM_Models/temozolomide_glm_med_lasso_model.rds')
temozolomide_med_elastic <- readRDS('GLM_Models/temozolomide_med_glm_alpha_search_model.rds')


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
lihc_gene_short_scaled <- fix_scale(lihc_gene_short_scaled)


## elastic
new_lihc_tcga_temozolomide <- predict(temozolomide_elastic, newx = as.matrix(lihc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

lihc_surv_times <- lihc_clinical_temozolomide_short$PFS
lihc_status <- ifelse(lihc_clinical_temozolomide_short$PFS == lihc_clinical_temozolomide_short$OS, 0, 1)

lihc_surv_df <- data.frame(lihc_surv_times, lihc_status, new_lihc_tcga_temozolomide)
fit <- survfit(Surv(lihc_surv_times, lihc_status) ~ new_lihc_tcga_temozolomide,
               data = lihc_surv_df)
fit2 <- survfit(Surv(lihc_surv_times, lihc_status) ~ new_lihc_tcga_temozolomide,
                data = lihc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
# all called resistant



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
skcm_gene_short_scaled <- fix_scale(skcm_gene_short_scaled)


## elastic
new_skcm_tcga_temozolomide <- predict(temozolomide_elastic, newx = as.matrix(skcm_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

skcm_surv_times <- skcm_clinical_temozolomide_short$PFS
skcm_status <- ifelse(skcm_clinical_temozolomide_short$PFS == skcm_clinical_temozolomide_short$OS, 0, 1)

skcm_surv_df <- data.frame(skcm_surv_times, skcm_status, new_skcm_tcga_temozolomide)
fit <- survfit(Surv(skcm_surv_times, skcm_status) ~ X1,
               data = skcm_surv_df)
fit2 <- survfit(Surv(skcm_surv_times, skcm_status) ~ X1,
                data = skcm_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(skcm_surv_df, file = 'Survival_Data/skcm_temozolomide_surv_df.txt', row.names = FALSE)
