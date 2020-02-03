### THIS SCRIPT FITS GLM MODELS

## load necessary packages ----
if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) # for model building

if(!require ('ROSE')) install.packages('ROSE')
library(ROSE)

### functions needed ----
# create function opposite of %in%
'%ni%' <- Negate('%in%')

### GDSC ------
## load clinical data ----
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

## load gene expression data ----
gdsc_rna_seq <- read.csv('Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE)
gdsc_rna_seq <- gdsc_rna_seq[, -1]
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))

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

# set GDSC to usable format and scale
gdsc <- data.frame(t(gdsc_rna_seq)) # puts predictors in columns
gdsc_names <- rownames(gdsc)
gdsc <- apply(gdsc, 2, scale)
rownames(gdsc) <- gdsc_names
# dim: 962 x 14209 

# split GDSC randomly
set.seed(5)
# get random numbers to use for split
random_sample <- sample(x = rownames(gdsc), size = nrow(gdsc)/5)

# get training and testing sets
gdsc_train         <- gdsc[which(rownames(gdsc) %ni% random_sample), ] #770 x 14209

gdsc_test          <- gdsc[random_sample, ] #192 x 14209

# make sure zero overlap
intersect(rownames(gdsc_train), rownames(gdsc_test)) #0

# create training/testing sets for each drug
bleomycin_rna_seq_train     <- gdsc_train[intersect(bleomycin_lines, rownames(gdsc_train)), ]
bleomycin_rna_seq_train_df     <- as.data.frame(bleomycin_rna_seq_train)
# 622 x 14209
bleomycin_rna_seq_test      <- gdsc_test[intersect(bleomycin_lines, rownames(gdsc_test)), ]
bleomycin_rna_seq_test      <- as.data.frame(bleomycin_rna_seq_test)
# 144 x 14209
camptothecin_rna_seq_train     <- gdsc_train[intersect(camptothecin_lines, rownames(gdsc_train)), ]
# 549 x 14209
camptothecin_rna_seq_test      <- gdsc_test[intersect(camptothecin_lines, rownames(gdsc_test)), ]
# 129 x 14209
cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# 551 x 14209
cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]
# 129 x 14209
cytarabine_rna_seq_train     <- gdsc_train[intersect(cytarabine_lines, rownames(gdsc_train)), ]
# 547 x 14209
cytarabine_rna_seq_test      <- gdsc_test[intersect(cytarabine_lines, rownames(gdsc_test)), ]
# 129 x 14209
doxorubicin_rna_seq_train     <- gdsc_train[intersect(doxorubicin_lines, rownames(gdsc_train)), ]
# 581 x 14209
doxorubicin_rna_seq_test      <- gdsc_test[intersect(doxorubicin_lines, rownames(gdsc_test)), ]
# 130 x 14209
etoposide_rna_seq_train     <- gdsc_train[intersect(etoposide_lines, rownames(gdsc_train)), ]
# 584 x 14209
etoposide_rna_seq_test      <- gdsc_test[intersect(etoposide_lines, rownames(gdsc_test)), ]
# 134 x 14209
gemcitabine_rna_seq_train     <- gdsc_train[intersect(gemcitabine_lines, rownames(gdsc_train)), ]
# 577 x 14209
gemcitabine_rna_seq_test      <- gdsc_test[intersect(gemcitabine_lines, rownames(gdsc_test)), ]
# 130 x 14209
methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# 550 x 14209
methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# 129 x 14209
mitomycin_rna_seq_train     <- gdsc_train[intersect(mitomycin_lines, rownames(gdsc_train)), ]
# 581 x 14209
mitomycin_rna_seq_test      <- gdsc_test[intersect(mitomycin_lines, rownames(gdsc_test)), ]
# 131 x 14209
sn38_rna_seq_train     <- gdsc_train[intersect(sn38_lines, rownames(gdsc_train)), ]
# 619 x 14209
sn38_rna_seq_test      <- gdsc_test[intersect(sn38_lines, rownames(gdsc_test)), ]
# 142 x 14209
temozolomide_rna_seq_train     <- gdsc_train[intersect(temozolomide_lines, rownames(gdsc_train)), ]
# 601 x 14209
temozolomide_rna_seq_test      <- gdsc_test[intersect(temozolomide_lines, rownames(gdsc_test)), ]
# 138 x 14209

# split clinical data
bleomycin_train        <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_train)), ]
bleomycin_rna_seq_train_df$res_sens <- bleomycin_train$res_sens
bleo_rose <- ROSE(res_sens ~ ., data = bleomycin_rna_seq_train_df)$data
bleo_rose_res_sens <- bleo_rose$res_sens
bleo_rose <- bleo_rose[, -14210]

bleomycin_test         <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_test)), ]
bleomycin_rna_seq_test$res_sens <- bleomycin_test$res_sens

camptothecin_train        <- camptothecin[which(camptothecin$COSMIC_ID %in% rownames(camptothecin_rna_seq_train)), ]
camptothecin_test         <- camptothecin[which(camptothecin$COSMIC_ID %in% rownames(camptothecin_rna_seq_test)), ]

cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

cytarabine_train        <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_train)), ]
cytarabine_test         <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_test)), ]

doxorubicin_train        <- doxorubicin[which(doxorubicin$COSMIC_ID %in% rownames(doxorubicin_rna_seq_train)), ]
doxorubicin_test         <- doxorubicin[which(doxorubicin$COSMIC_ID %in% rownames(doxorubicin_rna_seq_test)), ]

etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]

gemcitabine_train        <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test         <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_test)), ]

methotrexate_train     <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test      <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]

mitomycin_train        <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_train)), ]
mitomycin_test         <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_test)), ]

sn38_train        <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_train)), ]
sn38_test         <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_test)), ]

temozolomide_train        <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_train)), ]
temozolomide_test         <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_test)), ]


### fit models --------
## bleomycin
bleomycin_fit_elnet <- cv.glmnet(x = as.matrix(bleo_rose), y = bleo_rose_res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/bleomycin_auc.png')
plot(bleomycin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/bleomycin_model.rds', bleomycin_fit_elnet)
bleo_pred <- predict(bleomycin_fit_elnet, newx = as.matrix(bleomycin_rna_seq_test), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'class')
bleo_preds <- auc(bleomycin_test$res_sens, bleo_pred)
bleo_preds <- round(bleo_preds, digits = 2)
## camptothecin
camptothecin_fit_elnet <- cv.glmnet(x = as.matrix(camptothecin_rna_seq_train), y = camptothecin_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/camptothecin_auc.png')
plot(camptothecin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/camptothecin_model.rds', camptothecin_fit_elnet)

## cisplatin
cisplatin_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rna_seq_train), y = cisplatin_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/cisplatin_auc.png')
plot(cisplatin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/cisplatin_model.rds', cisplatin_fit_elnet)

## cytarabine
cytarabine_fit_elnet <- cv.glmnet(x = as.matrix(cytarabine_rna_seq_train), y = cytarabine_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/cytarabine_auc.png')
plot(cytarabine_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/cytarabine_model.rds', cytarabine_fit_elnet)

## doxorubicin
doxorubicin_fit_elnet <- cv.glmnet(x = as.matrix(doxorubicin_rna_seq_train), y = doxorubicin_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/doxorubicin_auc.png')
plot(doxorubicin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/doxorubicin_model.rds', doxorubicin_fit_elnet)

## etoposide
etoposide_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rna_seq_train), y = etoposide_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/etoposide_auc.png')
plot(etoposide_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/etoposide_model.rds', etoposide_fit_elnet)

## gemcitabine
gemcitabine_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train), y = gemcitabine_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/gemcitabine_auc.png')
plot(gemcitabine_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/gemcitabine_model.rds', gemcitabine_fit_elnet)

## methotrexate
methotrexate_fit_elnet <- cv.glmnet(x = as.matrix(meth_df_rose$data), y = meth_df_rose$data$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/methotrexate_auc.png')
plot(methotrexate_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/methotrexate_model.rds', methotrexate_fit_elnet)

## mitomycin
mitomycin_fit_elnet <- cv.glmnet(x = as.matrix(mitomycin_rna_seq_train), y = mitomycin_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/mitomycin_auc.png')
plot(mitomycin_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/mitomycin_model.rds', mitomycin_fit_elnet)

## sn38
sn38_fit_elnet <- cv.glmnet(x = as.matrix(sn38_rna_seq_train), y = sn38_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/sn38_auc.png')
plot(sn38_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/sn38_model.rds', sn38_fit_elnet)

## temozolomide
temozolomide_fit_elnet <- cv.glmnet(x = as.matrix(temozolomide_rna_seq_train), y = temozolomide_train$res_sens, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/temozolomide_auc.png')
plot(temozolomide_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/temozolomide_model.rds', temozolomide_fit_elnet)

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