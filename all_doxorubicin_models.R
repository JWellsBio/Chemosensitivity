# all doxorubicin models

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



doxorubicin     <- read.csv('Processed_Clinical_Data/doxorubicin_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

doxorubicin_lines           <- doxorubicin$COSMIC_ID #680


doxorubicin_rna_seq     <- gdsc[intersect(doxorubicin_lines, rownames(gdsc)), ]


doxorubicin_rna_seq <- as.data.frame(doxorubicin_rna_seq)
doxorubicin$med_sens <- ifelse(doxorubicin$LN_IC50 <= median(doxorubicin$LN_IC50), 1, 0)
doxorubicin_rna_seq$med_sens <- doxorubicin$med_sens

set.seed(5)
train_index <- createDataPartition(doxorubicin_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

doxorubicin_rna_seq_train <- doxorubicin_rna_seq[ train_index,]
doxorubicin_rna_seq_test  <- doxorubicin_rna_seq[-train_index,]

doxorubicin_elastic      <- readRDS('GLM_Models/doxorubicin_glm_alpha_search_model.rds')
doxorubicin_lasso <- readRDS('GLM_Models/doxorubicin_glm_lasso_model.rds')
doxorubicin_med_lasso <- readRDS('GLM_Models/doxorubicin_glm_med_lasso_model.rds')
doxorubicin_med_elastic <- readRDS('GLM_Models/doxorubicin_med_glm_alpha_search_model.rds')


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
brca_gene_short_scaled <- fix_scale(brca_gene_short_scaled)



## elastic
new_brca_tcga_doxorubicin <- predict(doxorubicin_elastic, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

brca_surv_times <- brca_clinical_doxorubicin_short$PFS
brca_status <- ifelse(brca_clinical_doxorubicin_short$PFS == brca_clinical_doxorubicin_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_doxorubicin)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_doxorubicin,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_doxorubicin,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(brca_surv_df, file = 'Survival_Data/brca_doxorubicin_surv_df.txt', row.names = FALSE)



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
ov_gene_short_scaled <- fix_scale(ov_gene_short_scaled)


## elastic
new_ov_tcga_doxorubicin <- predict(doxorubicin_elastic, newx = as.matrix(ov_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

ov_surv_times <- ov_clinical_doxorubicin_short$PFS
ov_status <- ifelse(ov_clinical_doxorubicin_short$PFS == ov_clinical_doxorubicin_short$OS, 0, 1)

ov_surv_df <- data.frame(ov_surv_times, ov_status, new_ov_tcga_doxorubicin)
fit <- survfit(Surv(ov_surv_times, ov_status) ~ new_ov_tcga_doxorubicin,
               data = ov_surv_df)
fit2 <- survfit(Surv(ov_surv_times, ov_status) ~ new_ov_tcga_doxorubicin,
                data = ov_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2, main = 'elastic')
legend(x = 2000, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)
# no events here



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
ucec_gene_short_scaled <- fix_scale(ucec_gene_short_scaled)


## elastic
new_ucec_tcga_doxorubicin <- predict(doxorubicin_elastic, newx = as.matrix(ucec_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

ucec_surv_times <- ucec_clinical_doxorubicin_short$PFS
ucec_status <- ifelse(ucec_clinical_doxorubicin_short$PFS == ucec_clinical_doxorubicin_short$OS, 0, 1)

ucec_surv_df <- data.frame(ucec_surv_times, ucec_status, new_ucec_tcga_doxorubicin)
fit <- survfit(Surv(ucec_surv_times, ucec_status) ~ new_ucec_tcga_doxorubicin,
               data = ucec_surv_df)
fit2 <- survfit(Surv(ucec_surv_times, ucec_status) ~ new_ucec_tcga_doxorubicin,
                data = ucec_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(ucec_surv_df, file = 'Survival_Data/ucec_doxorubicin_surv_df.txt', row.names = FALSE)

