#all gemcitabine models

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



gemcitabine     <- read.csv('Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

gemcitabine_lines           <- gemcitabine$COSMIC_ID #680


gemcitabine_rna_seq     <- gdsc[intersect(gemcitabine_lines, rownames(gdsc)), ]


gemcitabine_rna_seq <- as.data.frame(gemcitabine_rna_seq)
gemcitabine$med_sens <- ifelse(gemcitabine$LN_IC50 <= median(gemcitabine$LN_IC50), 1, 0)
gemcitabine_rna_seq$med_sens <- gemcitabine$med_sens

set.seed(5)
train_index <- createDataPartition(gemcitabine_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

gemcitabine_rna_seq_train <- gemcitabine_rna_seq[ train_index,]
gemcitabine_rna_seq_test  <- gemcitabine_rna_seq[-train_index,]

gemcitabine_elastic      <- readRDS('GLM_Models/gemcitabine_glm_alpha_search_model.rds')
gemcitabine_lasso <- readRDS('GLM_Models/gemcitabine_glm_lasso_model.rds')
gemcitabine_med_lasso <- readRDS('GLM_Models/gemcitabine_glm_med_lasso_model.rds')
gemcitabine_med_elastic <- readRDS('GLM_Models/gemcitabine_med_glm_alpha_search_model.rds')


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


## elastic
new_blca_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(blca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

blca_surv_times <- blca_clinical_gemcitabine_short$PFS
blca_status <- ifelse(blca_clinical_gemcitabine_short$PFS == blca_clinical_gemcitabine_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_gemcitabine)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ new_blca_tcga_gemcitabine,
                data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(blca_surv_df, file = 'Survival_Data/blca_gemcitabine_surv_df.txt', row.names = FALSE)



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
brca_gene_short_scaled <- fix_scale(brca_gene_short_scaled)


## elastic
new_brca_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(brca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

brca_surv_times <- brca_clinical_gemcitabine_short$PFS
brca_status <- ifelse(brca_clinical_gemcitabine_short$PFS == brca_clinical_gemcitabine_short$OS, 0, 1)

brca_surv_df <- data.frame(brca_surv_times, brca_status, new_brca_tcga_gemcitabine)
fit <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_gemcitabine,
               data = brca_surv_df)
fit2 <- survfit(Surv(brca_surv_times, brca_status) ~ new_brca_tcga_gemcitabine,
                data = brca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(brca_surv_df, file = 'Survival_Data/brca_gemcitabine_surv_df.txt', row.names = FALSE)



# CHOL W GEMCITABINE (8)
chol_clinical <- read.csv('Processed_Clinical_Data/chol_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(chol_clinical$most_sensitive)
chol_clinical <- chol_clinical[!na_idx, ]
table(chol_clinical$drug_name)
chol_clinical_gemcitabine <- chol_clinical[which(chol_clinical$drug_name == 'gemcitabine' | chol_clinical$drug_name == 'Gemcitabine' | 
                                                   chol_clinical$drug_name == 'GEMCITABINE' | chol_clinical$drug_name == 'Gemzar'), ]


chol_gene <- read.csv('Processed_Gene_Expression/chol_tcga_rna_seq_processed.csv', row.names = 1)
colnames(chol_gene) <- gsub('\\.', '-', colnames(chol_gene))
chol_matching_idx <- chol_clinical_gemcitabine$submitter_id.samples %in% colnames(chol_gene)
chol_clinical_gemcitabine_short <- chol_clinical_gemcitabine[chol_matching_idx, ]
chol_matching_idx <- colnames(chol_gene) %in% chol_clinical_gemcitabine_short$submitter_id.samples
chol_gene_short <- chol_gene[, chol_matching_idx]
chol_gene_short <- t(chol_gene_short)
chol_gene_short_scaled <- apply(chol_gene_short, 2, scale)
chol_gene_short_scaled <- fix_scale(chol_gene_short_scaled)


## elastic
new_chol_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(chol_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

chol_surv_times <- chol_clinical_gemcitabine_short$PFS
chol_status <- ifelse(chol_clinical_gemcitabine_short$PFS == chol_clinical_gemcitabine_short$OS, 0, 1)

chol_surv_df <- data.frame(chol_surv_times, chol_status, new_chol_tcga_gemcitabine)
fit <- survfit(Surv(chol_surv_times, chol_status) ~ X1,
               data = chol_surv_df)
fit2 <- survfit(Surv(chol_surv_times, chol_status) ~ X1,
                data = chol_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(chol_surv_df, file = 'Survival_Data/chol_gemcitabine_surv_df.txt', row.names = FALSE)



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
kirc_gene_short_scaled <- fix_scale(kirc_gene_short_scaled)


## elastic
new_kirc_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(kirc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

kirc_surv_times <- kirc_clinical_gemcitabine_short$PFS
kirc_status <- ifelse(kirc_clinical_gemcitabine_short$PFS == kirc_clinical_gemcitabine_short$OS, 0, 1)

kirc_surv_df <- data.frame(kirc_surv_times, kirc_status, new_kirc_tcga_gemcitabine)
fit <- survfit(Surv(kirc_surv_times, kirc_status) ~ X1,
               data = kirc_surv_df)
fit2 <- survfit(Surv(kirc_surv_times, kirc_status) ~ X1,
                data = kirc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2, main = 'elastic')
legend(x = 175, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)
# all predicted resistant


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
lihc_gene_short_scaled <- fix_scale(lihc_gene_short_scaled)


## elastic
new_lihc_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(lihc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

lihc_surv_times <- lihc_clinical_gemcitabine_short$PFS
lihc_status <- ifelse(lihc_clinical_gemcitabine_short$PFS == lihc_clinical_gemcitabine_short$OS, 0, 1)

lihc_surv_df <- data.frame(lihc_surv_times, lihc_status, new_lihc_tcga_gemcitabine)
fit <- survfit(Surv(lihc_surv_times, lihc_status) ~ X1,
               data = lihc_surv_df)
fit2 <- survfit(Surv(lihc_surv_times, lihc_status) ~ X1,
                data = lihc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(lihc_surv_df, file = 'Survival_Data/lihc_gemcitabine_surv_df.txt', row.names = FALSE)


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
luad_gene_short_scaled <- fix_scale(luad_gene_short_scaled)


## elastic
new_luad_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(luad_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

luad_surv_times <- luad_clinical_gemcitabine_short$PFS
luad_status <- ifelse(luad_clinical_gemcitabine_short$PFS == luad_clinical_gemcitabine_short$OS, 0, 1)

luad_surv_df <- data.frame(luad_surv_times, luad_status, new_luad_tcga_gemcitabine)
fit <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
               data = luad_surv_df)
fit2 <- survfit(Surv(luad_surv_times, luad_status) ~ X1,
                data = luad_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(luad_surv_df, file = 'Survival_Data/luad_gemcitabine_surv_df.txt', row.names = FALSE)


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
lusc_gene_short_scaled <- fix_scale(lusc_gene_short_scaled)


## elastic
new_lusc_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(lusc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

lusc_surv_times <- lusc_clinical_gemcitabine_short$PFS
lusc_status <- ifelse(lusc_clinical_gemcitabine_short$PFS == lusc_clinical_gemcitabine_short$OS, 0, 1)

lusc_surv_df <- data.frame(lusc_surv_times, lusc_status, new_lusc_tcga_gemcitabine)
fit <- survfit(Surv(lusc_surv_times, lusc_status) ~ X1,
               data = lusc_surv_df)
fit2 <- survfit(Surv(lusc_surv_times, lusc_status) ~ X1,
                data = lusc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(lusc_surv_df, file = 'Survival_Data/lusc_gemcitabine_surv_df.txt', row.names = FALSE)


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
ov_gene_short_scaled <- fix_scale(ov_gene_short_scaled)


## elastic
new_ov_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(ov_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

ov_surv_times <- ov_clinical_gemcitabine_short$PFS
ov_status <- ifelse(ov_clinical_gemcitabine_short$PFS == ov_clinical_gemcitabine_short$OS, 0, 1)

ov_surv_df <- data.frame(ov_surv_times, ov_status, new_ov_tcga_gemcitabine)
fit <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
               data = ov_surv_df)
fit2 <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
                data = ov_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(ov_surv_df, file = 'Survival_Data/ov_gemcitabine_surv_df.txt', row.names = FALSE)


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
paad_gene_short_scaled <- fix_scale(paad_gene_short_scaled)


## elastic
new_paad_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(paad_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

paad_surv_times <- paad_clinical_gemcitabine_short$PFS
paad_status <- ifelse(paad_clinical_gemcitabine_short$PFS == paad_clinical_gemcitabine_short$OS, 0, 1)

paad_surv_df <- data.frame(paad_surv_times, paad_status, new_paad_tcga_gemcitabine)
fit <- survfit(Surv(paad_surv_times, paad_status) ~ X1,
               data = paad_surv_df)
fit2 <- survfit(Surv(paad_surv_times, paad_status) ~ X1,
                data = paad_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(paad_surv_df, file = 'Survival_Data/paad_gemcitabine_surv_df.txt', row.names = FALSE)


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
sarc_gene_short_scaled <- fix_scale(sarc_gene_short_scaled)


## elastic
new_sarc_tcga_gemcitabine <- predict(gemcitabine_elastic, newx = as.matrix(sarc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

sarc_surv_times <- sarc_clinical_gemcitabine_short$PFS
sarc_status <- ifelse(sarc_clinical_gemcitabine_short$PFS == sarc_clinical_gemcitabine_short$OS, 0, 1)

sarc_surv_df <- data.frame(sarc_surv_times, sarc_status, new_sarc_tcga_gemcitabine)
fit <- survfit(Surv(sarc_surv_times, sarc_status) ~ X1,
               data = sarc_surv_df)
fit2 <- survfit(Surv(sarc_surv_times, sarc_status) ~ X1,
                data = sarc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(sarc_surv_df, file = 'Survival_Data/sarc_gemcitabine_surv_df.txt', row.names = FALSE)

