#all cisplatin models

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



cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)


## load gene expression data ----
gdsc <- read.csv('gdsc_rna_seq_names.csv', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names

cisplatin_lines           <- cisplatin$COSMIC_ID #680


cisplatin_rna_seq     <- gdsc[intersect(cisplatin_lines, rownames(gdsc)), ]


cisplatin_rna_seq <- as.data.frame(cisplatin_rna_seq)
cisplatin$med_sens <- ifelse(cisplatin$LN_IC50 <= median(cisplatin$LN_IC50), 1, 0)
cisplatin_rna_seq$med_sens <- cisplatin$med_sens

set.seed(5)
train_index <- createDataPartition(cisplatin_rna_seq$med_sens, p = .8, 
                                   list = FALSE, 
                                   times = 1)

cisplatin_rna_seq_train <- cisplatin_rna_seq[ train_index,]
cisplatin_rna_seq_test  <- cisplatin_rna_seq[-train_index,]

cisplatin_elastic      <- readRDS('GLM_Models/cisplatin_glm_alpha_search_model.rds')
cisplatin_lasso <- readRDS('GLM_Models/cisplatin_glm_lasso_model.rds')
cisplatin_med_lasso <- readRDS('GLM_Models/cisplatin_glm_med_lasso_model.rds')
cisplatin_med_elastic <- readRDS('GLM_Models/cisplatin_med_glm_alpha_search_model.rds')

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

colnames(blca_gene_short_scaled) <- colnames(gdsc)


## elastic
new_blca_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(blca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class', na.action = 'na.pass')

blca_surv_times <- blca_clinical_cisplatin_short$PFS
blca_status <- ifelse(blca_clinical_cisplatin_short$PFS == blca_clinical_cisplatin_short$OS, 0, 1)

blca_surv_df <- data.frame(blca_surv_times, blca_status, new_blca_tcga_cisplatin)
fit <- survfit(Surv(blca_surv_times, blca_status) ~ X1,
               data = blca_surv_df)
fit2 <- survfit(Surv(blca_surv_times, blca_status) ~ X1,
                data = blca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(blca_surv_df, file = 'Survival_Data/blca_cisplatin_surv.df.txt', row.names = FALSE)


par(mfrow = c(2,2))
# CESC W CISPLATIN (104)
cesc_clinical <- read.csv('Processed_Clinical_Data/cesc_tcga_clinical_processed.csv', row.names = 1)
na_idx <- is.na(cesc_clinical$most_sensitive)
cesc_clinical <- cesc_clinical[!na_idx, ]
table(cesc_clinical$drug_name)
cesc_clinical_cisplatin <- cesc_clinical[which(cesc_clinical$drug_name == 'cisplatin' | cesc_clinical$drug_name == 'Cisplatin' | 
                                                 cesc_clinical$drug_name == 'Cisplatin-xrt' | cesc_clinical$drug_name == 'Cisplatinum'), ]


cesc_gene <- read.csv('Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv', row.names = 1)
colnames(cesc_gene) <- gsub('\\.', '-', colnames(cesc_gene))
cesc_matching_idx <- cesc_clinical_cisplatin$submitter_id.samples %in% colnames(cesc_gene)
cesc_clinical_cisplatin_short <- cesc_clinical_cisplatin[cesc_matching_idx, ]
cesc_matching_idx <- colnames(cesc_gene) %in% cesc_clinical_cisplatin_short$submitter_id.samples
cesc_gene_short <- cesc_gene[, cesc_matching_idx]
cesc_gene_short <- t(cesc_gene_short)

cesc_gene_short_scaled <- apply(cesc_gene_short, 2, scale)
cesc_gene_short_scaled <- scale(cesc_gene_short, center = TRUE, scale = TRUE)
colnames(cesc_gene_short_scaled) <- colnames(gdsc)

cesc_gene_short_scaled <- fix_scale(cesc_gene_short_scaled)

new_cesc_tcga_cisplatin <- predict(cisplatin_elastic, newx = cesc_gene_short_scaled, s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

cesc_surv_times <- cesc_clinical_cisplatin_short$PFS
cesc_status <- ifelse(cesc_clinical_cisplatin_short$PFS == cesc_clinical_cisplatin_short$OS, 0, 1)

cesc_surv_df <- data.frame(cesc_surv_times, cesc_status, new_cesc_tcga_cisplatin)
fit <- survfit(Surv(cesc_surv_times, cesc_status) ~ X1,
               data = cesc_surv_df)
fit2 <- survfit(Surv(cesc_surv_times, cesc_status) ~ X1,
                data = cesc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(cesc_surv_df, file = 'Survival_Data/cesc_cisplatin_surv_df.txt', row.names = FALSE)



# ESCA W CISPLATIN (15)
par(mfrow = c(2,2))

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

colnames(esca_gene_short_scaled) <- colnames(gdsc)

new_esca_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(esca_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

esca_surv_times <- esca_clinical_cisplatin_short$PFS
esca_status <- ifelse(esca_clinical_cisplatin_short$PFS == esca_clinical_cisplatin_short$OS, 0, 1)

esca_surv_df <- data.frame(esca_surv_times, esca_status, new_esca_tcga_cisplatin)
fit <- survfit(Surv(esca_surv_times, esca_status) ~ new_esca_tcga_cisplatin,
               data = esca_surv_df)
fit2 <- survfit(Surv(esca_surv_times, esca_status) ~ new_esca_tcga_cisplatin,
                data = esca_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(esca_surv_df, file = 'Survival_Data/esca_cisplatin_surv_df.txt', row.names = FALSE)


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

par(mfrow = c(2,2))

## elastic
new_hnsc_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(hnsc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

hnsc_surv_times <- hnsc_clinical_cisplatin_short$PFS
hnsc_status <- ifelse(hnsc_clinical_cisplatin_short$PFS == hnsc_clinical_cisplatin_short$OS, 0, 1)

hnsc_surv_df <- data.frame(hnsc_surv_times, hnsc_status, new_hnsc_tcga_cisplatin)
fit <- survfit(Surv(hnsc_surv_times, hnsc_status) ~ X1,
               data = hnsc_surv_df)
fit2 <- survfit(Surv(hnsc_surv_times, hnsc_status) ~ X1,
                data = hnsc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt

write.table(hnsc_surv_df, file = 'Survival_Data/hnsc_cisplatin_surv_df.txt', row.names = FALSE)




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
luad_gene_short_scaled <- fix_scale(luad_gene_short_scaled)


## elastic
new_luad_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(luad_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

luad_surv_times <- luad_clinical_cisplatin_short$PFS
luad_status <- ifelse(luad_clinical_cisplatin_short$PFS == luad_clinical_cisplatin_short$OS, 0, 1)

luad_surv_df <- data.frame(luad_surv_times, luad_status, new_luad_tcga_cisplatin)
fit <- survfit(Surv(luad_surv_times, luad_status) ~ new_luad_tcga_cisplatin,
               data = luad_surv_df)
fit2 <- survfit(Surv(luad_surv_times, luad_status) ~ new_luad_tcga_cisplatin,
                data = luad_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(luad_surv_df, file = 'Survival_Data/luad_cisplatin_surv_df.txt', row.names = FALSE)




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
lusc_gene_short_scaled <- fix_scale(lusc_gene_short_scaled)


## elastic
new_lusc_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(lusc_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

lusc_surv_times <- lusc_clinical_cisplatin_short$PFS
lusc_status <- ifelse(lusc_clinical_cisplatin_short$PFS == lusc_clinical_cisplatin_short$OS, 0, 1)

lusc_surv_df <- data.frame(lusc_surv_times, lusc_status, new_lusc_tcga_cisplatin)
fit <- survfit(Surv(lusc_surv_times, lusc_status) ~ new_lusc_tcga_cisplatin,
               data = lusc_surv_df)
fit2 <- survfit(Surv(lusc_surv_times, lusc_status) ~ new_lusc_tcga_cisplatin,
                data = lusc_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(lusc_surv_df, file = 'Survival_Data/lusc_cisplatin_surv_df.txt', row.names = FALSE)




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
meso_gene_short_scaled <- fix_scale(meso_gene_short_scaled)


## elastic
new_meso_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(meso_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

meso_surv_times <- meso_clinical_cisplatin_short$PFS
meso_status <- ifelse(meso_clinical_cisplatin_short$PFS == meso_clinical_cisplatin_short$OS, 0, 1)

meso_surv_df <- data.frame(meso_surv_times, meso_status, new_meso_tcga_cisplatin)
fit <- survfit(Surv(meso_surv_times, meso_status) ~ new_meso_tcga_cisplatin,
               data = meso_surv_df)
fit2 <- survfit(Surv(meso_surv_times, meso_status) ~ new_meso_tcga_cisplatin,
                data = meso_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(meso_surv_df, file = 'Survival_Data/meso_cisplatin_surv_df.txt', row.names = FALSE)
plot(fit, col = c('limegreen', 'darkviolet'), xlab = 'Time (d)', ylab = 'Percent Recurrence-Free', lwd = 2, main = 'elastic')
legend(x = 100, y = 0.8, legend = paste0('log-rank\n', fit_pvalue), bty = 'n', cex = 0.8)
legend(x = 450, y = 0.4, legend = c('predicted resistant', 'predicted sensitive'), lty = c(1,1), lwd = 2, col = c('darkviolet', 'limegreen'), bty = 'n', cex = 0.8)



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
ov_gene_short_scaled <- fix_scale(ov_gene_short_scaled)


## elastic
new_ov_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(ov_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

ov_surv_times <- ov_clinical_cisplatin_short$PFS
ov_status <- ifelse(ov_clinical_cisplatin_short$PFS == ov_clinical_cisplatin_short$OS, 0, 1)

ov_surv_df <- data.frame(ov_surv_times, ov_status, new_ov_tcga_cisplatin)
fit <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
               data = ov_surv_df)
fit2 <- survfit(Surv(ov_surv_times, ov_status) ~ X1,
                data = ov_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(ov_surv_df, file = 'Survival_Data/ov_cisplatin_surv_df.txt', row.names = FALSE)


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
stad_gene_short_scaled <- fix_scale(stad_gene_short_scaled)


## elastic
new_stad_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(stad_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

stad_surv_times <- stad_clinical_cisplatin_short$PFS
stad_status <- ifelse(stad_clinical_cisplatin_short$PFS == stad_clinical_cisplatin_short$OS, 0, 1)

stad_surv_df <- data.frame(stad_surv_times, stad_status, new_stad_tcga_cisplatin)
fit <- survfit(Surv(stad_surv_times, stad_status) ~ new_stad_tcga_cisplatin,
               data = stad_surv_df)
fit2 <- survfit(Surv(stad_surv_times, stad_status) ~ new_stad_tcga_cisplatin,
                data = stad_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(stad_surv_df, file = 'Survival_Data/stad_cisplatin_surv_df.txt', row.names = FALSE)


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


## elastic
new_tgct_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(tgct_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

tgct_surv_times <- tgct_clinical_cisplatin_short$PFS
tgct_status <- ifelse(tgct_clinical_cisplatin_short$PFS == tgct_clinical_cisplatin_short$OS, 0, 1)

tgct_surv_df <- data.frame(tgct_surv_times, tgct_status, new_tgct_tcga_cisplatin)
fit <- survfit(Surv(tgct_surv_times, tgct_status) ~ X1,
               data = tgct_surv_df)
fit2 <- survfit(Surv(tgct_surv_times, tgct_status) ~ X1,
                data = tgct_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(tgct_surv_df, file = 'Survival_Data/tgct_cisplatin_surv_df.txt', row.names = FALSE)



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
ucec_gene_short_scaled <- fix_scale(ucec_gene_short_scaled)


## elastic
new_ucec_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(ucec_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

ucec_surv_times <- ucec_clinical_cisplatin_short$PFS
ucec_status <- ifelse(ucec_clinical_cisplatin_short$PFS == ucec_clinical_cisplatin_short$OS, 0, 1)

ucec_surv_df <- data.frame(ucec_surv_times, ucec_status, new_ucec_tcga_cisplatin)
fit <- survfit(Surv(ucec_surv_times, ucec_status) ~ new_ucec_tcga_cisplatin,
               data = ucec_surv_df)
fit2 <- survfit(Surv(ucec_surv_times, ucec_status) ~ new_ucec_tcga_cisplatin,
                data = ucec_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
# no censoring or recurrence events here


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
ucs_gene_short_scaled <- fix_scale(ucs_gene_short_scaled)


## elastic
new_ucs_tcga_cisplatin <- predict(cisplatin_elastic, newx = as.matrix(ucs_gene_short_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')

ucs_surv_times <- ucs_clinical_cisplatin_short$PFS
ucs_status <- ifelse(ucs_clinical_cisplatin_short$PFS == ucs_clinical_cisplatin_short$OS, 0, 1)

ucs_surv_df <- data.frame(ucs_surv_times, ucs_status, new_ucs_tcga_cisplatin)
fit <- survfit(Surv(ucs_surv_times, ucs_status) ~ new_ucs_tcga_cisplatin,
               data = ucs_surv_df)
fit2 <- survfit(Surv(ucs_surv_times, ucs_status) ~ new_ucs_tcga_cisplatin,
                data = ucs_surv_df)
fit_pvalue <- surv_pvalue(fit)$pval.txt
write.table(ucs_surv_df, file = 'Survival_Data/ucs_cisplatin_surv_df.txt', row.names = FALSE)

