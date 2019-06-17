### THIS SCRIPT FITS GLM MODELS

## load necessary packages ----
if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) # for model building

### functions needed ----
# create function opposite of %in%
'%ni%' <- Negate('%in%')

### GDSC ------
## load clinical data ----
cisplatin     <- read.csv('Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv', row.names = 1)
etoposide     <- read.csv('Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv', row.names = 1)
methotrexate  <- read.csv('Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv', row.names = 1)

## load gene expression data ----
gdsc_rna_seq <- read.csv('Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE)
gdsc_rna_seq <- gdsc_rna_seq[, -1]
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))

### set up data for model building ----------
# get names of GDSC cell lines treated with each drug
cisplatin_lines           <- cisplatin$COSMIC_ID #680
etoposide_lines           <- etoposide$COSMIC_ID #718
methotrexate_lines        <- methotrexate$COSMIC_ID # 679

# set GDSC to usable format
gdsc <- data.frame(t(gdsc_rna_seq)) # puts predictors in columns
rownames(gdsc) <- gsub('X', '', rownames(gdsc))
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
cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# 338 x 14209
cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]
# 342 x 14209
etoposide_rna_seq_train     <- gdsc_train[intersect(etoposide_lines, rownames(gdsc_train)), ]
# 352 x 14209
etoposide_rna_seq_test      <- gdsc_test[intersect(etoposide_lines, rownames(gdsc_test)), ]
# 366 x 14209
methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# 337 x 14209
methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# 342 x 14209

# split clinical data
cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]

methotrexate_train     <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test      <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]

# scale data
cisplatin_rna_seq_train_scaled          <- apply(cisplatin_rna_seq_train, 2, scale)
cisplatin_rna_seq_test_scaled           <- as.data.frame(apply(cisplatin_rna_seq_test, 2, scale))

etoposide_rna_seq_train_scaled          <- apply(etoposide_rna_seq_train, 2, scale)
etoposide_rna_seq_test_scaled           <- as.data.frame(apply(etoposide_rna_seq_test, 2, scale))

methotrexate_rna_seq_train_scaled       <- apply(methotrexate_rna_seq_train, 2, scale)
methotrexate_rna_seq_test_scaled        <- as.data.frame(apply(methotrexate_rna_seq_test, 2, scale))

### fit models --------
## CISPLATIN
cisplatin_most_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rna_seq_train_scaled), y = cisplatin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/cisplatin_most_auc.png')
plot(cisplatin_most_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/cisplatin_most_model.rds', cisplatin_most_fit_elnet)

cisplatin_least_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rna_seq_train_scaled), y = cisplatin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cisplatin_least_auc.png')
plot(cisplatin_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/cisplatin_least_model.rds', cisplatin_least_fit_elnet)

## ETOPOSIDE
etoposide_most_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rna_seq_train_scaled), y = etoposide_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/etoposide_most_auc.png')
plot(etoposide_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/etoposide_most_model.rds', etoposide_most_fit_elnet)

etoposide_least_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rna_seq_train_scaled), y = etoposide_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/etoposide_least_auc.png')
plot(etoposide_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/etoposide_least_model.rds', etoposide_least_fit_elnet)

## METHOTREXATE
methotrexate_most_fit_elnet <- cv.glmnet(x = as.matrix(methotrexate_rna_seq_train_scaled), y = methotrexate_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/methotrexate_most_auc.png')
plot(methotrexate_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/methotrexate_most_model.rds', methotrexate_most_fit_elnet)

methotrexate_least_fit_elnet <- cv.glmnet(x = as.matrix(methotrexate_rna_seq_train_scaled), y = methotrexate_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/methotrexate_least_auc.png')
plot(methotrexate_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/methotrexate_least_model.rds', methotrexate_least_fit_elnet)

### CCLE -----
## load clinical data ----
carboplatin_ccle <- read.csv('Processed_Clinical_Data/carboplatin_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
cyclophosphamide_ccle <- read.csv('Processed_Clinical_Data/cyclophosphamide_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
docetaxel_ccle <- read.csv('Processed_Clinical_Data/docetaxel_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
fluorouracil_ccle <- read.csv('Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)
gemcitabine_ccle <- read.csv('Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv', row.names = 1, stringsAsFactors = FALSE)

## load gene expression data -----
ccle_microarray <- read.csv('Processed_Gene_Expression/ccle_microarray_processed.csv', row.names = 1)
colnames(ccle_microarray) <- gsub('_.*$', '', colnames(ccle_microarray))
colnames(ccle_microarray) <- gsub('^X', '', colnames(ccle_microarray))

## set up data for model building ----
carboplatin_lines       <- carboplatin_ccle$Cell.line.name #650
carboplatin_lines       <- as.character(carboplatin_lines)
cyclophosphamide_lines  <- cyclophosphamide_ccle$Cell.line.name #335
cyclophosphamide_lines  <- as.character(cyclophosphamide_lines)
docetaxel_lines         <- docetaxel_ccle$Cell.line.name #302
docetaxel_lines         <- as.character(docetaxel_lines)
fluorouracil_lines      <- fluorouracil_ccle$Cell.line.name #652
fluorouracil_lines      <- as.character(fluorouracil_lines)
gemcitabine_lines       <- gemcitabine_ccle$Cell.line.name #589
gemcitabine_lines       <- as.character(gemcitabine_lines)

# set ccle dataframe correct orientation
ccle <- data.frame(t(ccle_microarray))
#855 x 14209

#split into training/testing
set.seed(5)

random_sample <- sample(x = rownames(ccle), size = nrow(ccle)/5)


ccle_train         <- ccle[which(rownames(ccle) %ni% random_sample), ] #427 x 14209

ccle_test          <- ccle[random_sample, ] #481 x 14209

intersect(rownames(ccle_train), rownames(ccle_test)) #0

#split expression data
carboplatin_rna_seq_train           <- ccle_train[intersect(carboplatin_lines, rownames(ccle_train)), ]
# 312 x 14209
carboplatin_rna_seq_test            <- ccle_test[intersect(carboplatin_lines, rownames(ccle_test)), ]
# 306 x 14209

cyclophosphamide_rna_seq_train      <- ccle_train[intersect(cyclophosphamide_lines, rownames(ccle_train)), ]
# 165 x 14209
cyclophosphamide_rna_seq_test       <- ccle_test[intersect(cyclophosphamide_lines, rownames(ccle_test)), ]
# 151 x 14209

docetaxel_rna_seq_train             <- ccle_train[intersect(docetaxel_lines, rownames(ccle_train)), ]
# 142 x 14209
docetaxel_rna_seq_test              <- ccle_test[intersect(docetaxel_lines, rownames(ccle_test)), ]
# 155 x 14209

fluorouracil_rna_seq_train          <- ccle_train[intersect(fluorouracil_lines, rownames(ccle_train)), ]
# 313 x 14209
fluorouracil_rna_seq_test           <- ccle_test[intersect(fluorouracil_lines, rownames(ccle_test)), ]
# 307 x 14209

gemcitabine_rna_seq_train           <- ccle_train[intersect(gemcitabine_lines, rownames(ccle_train)), ]
# 276 x 14209
gemcitabine_rna_seq_test            <- ccle_test[intersect(gemcitabine_lines, rownames(ccle_test)), ]
# 286 x 14209

#split clinical data
carboplatin_train         <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_train)), ]
carboplatin_test          <- carboplatin_ccle[which(carboplatin_ccle$Cell.line.name %in% rownames(carboplatin_rna_seq_test)), ]

cyclophosphamide_train    <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_train)), ]
cyclophosphamide_test     <- cyclophosphamide_ccle[which(cyclophosphamide_ccle$Cell.line.name %in% rownames(cyclophosphamide_rna_seq_test)), ]

docetaxel_train           <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_train)), ]
docetaxel_test            <- docetaxel_ccle[which(docetaxel_ccle$Cell.line.name %in% rownames(docetaxel_rna_seq_test)), ]

fluorouracil_train        <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_train)), ]
fluorouracil_test         <- fluorouracil_ccle[which(fluorouracil_ccle$Cell.line.name %in% rownames(fluorouracil_rna_seq_test)), ]

gemcitabine_train         <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test          <- gemcitabine_ccle[which(gemcitabine_ccle$Cell.line.name %in% rownames(gemcitabine_rna_seq_test)), ]

#scale data
carboplatin_rna_seq_train_scaled          <- apply(carboplatin_rna_seq_train, 2, scale)
carboplatin_rna_seq_test_scaled           <- as.data.frame(apply(carboplatin_rna_seq_test, 2, scale))

cyclophosphamide_rna_seq_train_scaled     <- apply(cyclophosphamide_rna_seq_train, 2, scale)
cyclophosphamide_rna_seq_test_scaled      <- as.data.frame(apply(cyclophosphamide_rna_seq_test, 2, scale))

docetaxel_rna_seq_train_scaled            <- apply(docetaxel_rna_seq_train, 2, scale)
docetaxel_rna_seq_test_scaled             <- as.data.frame(apply(docetaxel_rna_seq_test, 2, scale))

fluorouracil_rna_seq_train_scaled         <- apply(fluorouracil_rna_seq_train, 2, scale)
fluorouracil_rna_seq_test_scaled          <- as.data.frame(apply(fluorouracil_rna_seq_test, 2, scale))

gemcitabine_rna_seq_train_scaled          <- apply(gemcitabine_rna_seq_train, 2, scale)
gemcitabine_rna_seq_test_scaled           <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))


### build models ----
## carboplatin
carboplatin_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(carboplatin_rna_seq_train_scaled), y = carboplatin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
#save plot
png(filename = 'Images/carboplatin_ccle_most_auc.png')
plot(carboplatin_ccle_most_fit_elnet)
dev.off()
#save model
saveRDS(file = 'GLM_Models/carboplatin_ccle_most_model.rds', carboplatin_ccle_most_fit_elnet)

carboplatin_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(carboplatin_rna_seq_train_scaled), y = carboplatin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/carboplatin_ccle_least_auc.png')
plot(carboplatin_ccle_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/carboplatin_ccle_least_model.rds', carboplatin_ccle_least_fit_elnet)

## cyclophosphamide
cyclophosphamide_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(cyclophosphamide_rna_seq_train_scaled), y = cyclophosphamide_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cyclophosphamide_ccle_most_auc.png')
plot(cyclophosphamide_ccle_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/cyclophosphamide_ccle_most_model.rds', cyclophosphamide_ccle_most_fit_elnet)

cyclophosphamide_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(cyclophosphamide_rna_seq_train_scaled), y = cyclophosphamide_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cyclophosphamide_ccle_least_auc.png')
plot(cyclophosphamide_ccle_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/cyclophosphamide_ccle_least_model.rds', cyclophosphamide_ccle_least_fit_elnet)

## docetaxel
docetaxel_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(docetaxel_rna_seq_train_scaled), y = docetaxel_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/docetaxel_ccle_most_auc.png')
plot(docetaxel_ccle_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/docetaxel_ccle_most_model.rds', docetaxel_ccle_most_fit_elnet)

docetaxel_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(docetaxel_rna_seq_train_scaled), y = docetaxel_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/docetaxel_ccle_least_auc.png')
plot(docetaxel_ccle_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/docetaxel_ccle_least_model.rds', docetaxel_ccle_least_fit_elnet)

## fluorouracil
fluorouracil_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(fluorouracil_rna_seq_train_scaled), y = fluorouracil_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/fluorouracil_ccle_most_auc.png')
plot(fluorouracil_ccle_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/fluorouracil_ccle_most_model.rds', fluorouracil_ccle_most_fit_elnet)

fluorouracil_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(fluorouracil_rna_seq_train_scaled), y = fluorouracil_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/fluorouracil_ccle_least_auc.png')
plot(fluorouracil_ccle_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/fluorouracil_ccle_least_model.rds', fluorouracil_ccle_least_fit_elnet)

## gemcitabine
gemcitabine_ccle_most_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/gemcitabine_ccle_most_auc.png')
plot(gemcitabine_ccle_most_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/gemcitabine_ccle_most_model.rds', gemcitabine_ccle_most_fit_elnet)

gemcitabine_ccle_least_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/gemcitabine_ccle_least_auc.png')
plot(gemcitabine_ccle_least_fit_elnet)
dev.off()

saveRDS(file = 'GLM_Models/gemcitabine_ccle_least_model.rds', gemcitabine_ccle_least_fit_elnet)