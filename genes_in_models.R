## THIS SCRIPT SIMPLY CAPTURES THE GENES USED IN EACH MODEL AND SAVES THEM

### load necessary packages ----
if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) #for handling glm models already built

if (!require ('treemap')) install.packages('treemap')
library(treemap) #for total gene set figure

if (!require ('flextable')) install.packages('flextable')
library(flextable) #for nice tables


## GDSC -----
### load models ----
cisplatin_most_fit_elnet      <- readRDS('GLM_Models/cisplatin_most_model.rds')
cisplatin_least_fit_elnet     <- readRDS('GLM_Models/cisplatin_least_model.rds')

etoposide_most_fit_elnet      <- readRDS('GLM_Models/etoposide_most_model.rds')
etoposide_least_fit_elnet     <- readRDS('GLM_Models/etoposide_least_model.rds')

methotrexate_most_fit_elnet   <- readRDS('GLM_Models/methotrexate_most_model.rds')
methotrexate_least_fit_elnet  <- readRDS('GLM_Models/methotrexate_least_model.rds')

#### capture genes used in models --------
## CISPLATIN MOST SENSITIVE
cisplatin_most_min_tmp_coefs <- coef(cisplatin_most_fit_elnet, s = 'lambda.min')
cisplatin_most_model_min <- data.frame(name = cisplatin_most_min_tmp_coefs@Dimnames[[1]][cisplatin_most_min_tmp_coefs@i + 1], coefficient = cisplatin_most_min_tmp_coefs@x)
write.csv(cisplatin_most_model_min, file = 'GLM_Models/cisplatin_most_model_min.csv', row.names = FALSE)

cisplatin_most_1se_tmp_coefs <- coef(cisplatin_most_fit_elnet, s = 'lambda.1se')
cisplatin_most_model_1se <- data.frame(name = cisplatin_most_1se_tmp_coefs@Dimnames[[1]][cisplatin_most_1se_tmp_coefs@i + 1], coefficient = cisplatin_most_1se_tmp_coefs@x)
write.csv(cisplatin_most_model_1se, file = 'GLM_Models/cisplatin_most_model_1se.csv', row.names = FALSE)

## CISPLATIN LEAST SENSITIVE
cisplatin_least_min_tmp_coefs <- coef(cisplatin_least_fit_elnet, s = 'lambda.min')
cisplatin_least_model_min <- data.frame(name = cisplatin_least_min_tmp_coefs@Dimnames[[1]][cisplatin_least_min_tmp_coefs@i + 1], coefficient = cisplatin_least_min_tmp_coefs@x)
write.csv(cisplatin_least_model_min, file = 'GLM_Models/cisplatin_least_model_min.csv', row.names = FALSE)

cisplatin_least_1se_tmp_coefs <- coef(cisplatin_least_fit_elnet, s = 'lambda.1se')
cisplatin_least_model_1se <- data.frame(name = cisplatin_least_1se_tmp_coefs@Dimnames[[1]][cisplatin_least_1se_tmp_coefs@i + 1], coefficient = cisplatin_least_1se_tmp_coefs@x)
write.csv(cisplatin_least_model_1se, file = 'GLM_Models/cisplatin_least_model_1se.csv', row.names = FALSE)

## ETOPOSIDE MOST SENSITIVE
etoposide_most_min_tmp_coefs <- coef(etoposide_most_fit_elnet, s = 'lambda.min')
etoposide_most_model_min <- data.frame(name = etoposide_most_min_tmp_coefs@Dimnames[[1]][etoposide_most_min_tmp_coefs@i + 1], coefficient = etoposide_most_min_tmp_coefs@x)
write.csv(etoposide_most_model_min, file = 'GLM_Models/etoposide_most_model_min.csv', row.names = FALSE)

etoposide_most_1se_tmp_coefs <- coef(etoposide_most_fit_elnet, s = 'lambda.1se')
etoposide_most_model_1se <- data.frame(name = etoposide_most_1se_tmp_coefs@Dimnames[[1]][etoposide_most_1se_tmp_coefs@i + 1], coefficient = etoposide_most_1se_tmp_coefs@x)
write.csv(etoposide_most_model_1se, file = 'GLM_Models/etoposide_most_model_1se.csv', row.names = FALSE)

## ETOPOSIDE LEAST SENSITIVE
etoposide_least_min_tmp_coefs <- coef(etoposide_least_fit_elnet, s = 'lambda.min')
etoposide_least_model_min <- data.frame(name = etoposide_least_min_tmp_coefs@Dimnames[[1]][etoposide_least_min_tmp_coefs@i + 1], coefficient = etoposide_least_min_tmp_coefs@x)
write.csv(etoposide_least_model_min, file = 'GLM_Models/etoposide_least_model_min.csv', row.names = FALSE)

etoposide_least_1se_tmp_coefs <- coef(etoposide_least_fit_elnet, s = 'lambda.1se')
etoposide_least_model_1se <- data.frame(name = etoposide_least_1se_tmp_coefs@Dimnames[[1]][etoposide_least_1se_tmp_coefs@i + 1], coefficient = etoposide_least_1se_tmp_coefs@x)
write.csv(etoposide_least_model_1se, file = 'GLM_Models/etoposide_least_model_1se.csv', row.names = FALSE)

## METHOTREXATE MOST SENSITIVE
methotrexate_most_min_tmp_coefs <- coef(methotrexate_most_fit_elnet, s = 'lambda.min')
methotrexate_most_model_min <- data.frame(name = methotrexate_most_min_tmp_coefs@Dimnames[[1]][methotrexate_most_min_tmp_coefs@i + 1], coefficient = methotrexate_most_min_tmp_coefs@x)
write.csv(methotrexate_most_model_min, file = 'GLM_Models/methotrexate_most_model_min.csv', row.names = FALSE)

methotrexate_most_1se_tmp_coefs <- coef(methotrexate_most_fit_elnet, s = 'lambda.1se')
methotrexate_most_model_1se <- data.frame(name = methotrexate_most_1se_tmp_coefs@Dimnames[[1]][methotrexate_most_1se_tmp_coefs@i + 1], coefficient = methotrexate_most_1se_tmp_coefs@x)
write.csv(methotrexate_most_model_1se, file = 'GLM_Models/methotrexate_most_model_1se.csv', row.names = FALSE)

## METHOTREXATE LEAST SENSITIVE
methotrexate_least_min_tmp_coefs <- coef(methotrexate_least_fit_elnet, s = 'lambda.min')
methotrexate_least_model_min <- data.frame(name = methotrexate_least_min_tmp_coefs@Dimnames[[1]][methotrexate_least_min_tmp_coefs@i + 1], coefficient = methotrexate_least_min_tmp_coefs@x)
write.csv(methotrexate_least_model_min, file = 'GLM_Models/methotrexate_least_model_min.csv', row.names = FALSE)

methotrexate_least_1se_tmp_coefs <- coef(methotrexate_least_fit_elnet, s = 'lambda.1se')
methotrexate_least_model_1se <- data.frame(name = methotrexate_least_1se_tmp_coefs@Dimnames[[1]][methotrexate_least_1se_tmp_coefs@i + 1], coefficient = methotrexate_least_1se_tmp_coefs@x)
write.csv(methotrexate_least_model_1se, file = 'GLM_Models/methotrexate_least_model_1se.csv', row.names = FALSE)


## PUT TOGETHER IN TABLE ----
cisplatin_most_min_genes <- as.character(cisplatin_most_model_min$name)
cisplatin_most_min_genes <- cisplatin_most_min_genes[-1]

cisplatin_most_1se_genes <- as.character(cisplatin_most_model_1se$name)
cisplatin_most_1se_genes <- cisplatin_most_1se_genes[-1]

cisplatin_least_min_genes <- as.character(cisplatin_least_model_min$name)
cisplatin_least_min_genes <- cisplatin_least_min_genes[-1]

cisplatin_least_1se_genes <- as.character(cisplatin_least_model_1se$name)
cisplatin_least_1se_genes <- cisplatin_least_1se_genes[-1]

etoposide_most_min_genes <- as.character(etoposide_most_model_min$name)
etoposide_most_min_genes <- etoposide_most_min_genes[-1]

etoposide_most_1se_genes <- as.character(etoposide_most_model_1se$name)
etoposide_most_1se_genes <- etoposide_most_1se_genes[-1]

etoposide_least_min_genes <- as.character(etoposide_least_model_min$name)
etoposide_least_min_genes <- etoposide_least_min_genes[-1]

etoposide_least_1se_genes <- as.character(etoposide_least_model_1se$name)
etoposide_least_1se_genes <- etoposide_least_1se_genes[-1]

methotrexate_most_min_genes <- as.character(methotrexate_most_model_min$name)
methotrexate_most_min_genes <- methotrexate_most_min_genes[-1]

methotrexate_most_1se_genes <- as.character(methotrexate_most_model_1se$name)
methotrexate_most_1se_genes <- methotrexate_most_1se_genes[-1]

methotrexate_least_min_genes <- as.character(methotrexate_least_model_min$name)
methotrexate_least_min_genes <- methotrexate_least_min_genes[-1]

methotrexate_least_1se_genes <- as.character(methotrexate_least_model_1se$name)
methotrexate_least_1se_genes <- methotrexate_least_1se_genes[-1]

# using: cisplatin most 1se, cisplatin least min, etoposide most min, etoposide least min,
# gemcitabine most 1se, gemcitabine least 1se, methotrexate most min, methotrexate least 1se

length(intersect(cisplatin_most_1se_genes, cisplatin_least_min_genes))
#4, gives total of 103 for cisplatin
length(intersect(etoposide_most_min_genes, etoposide_least_min_genes))
#5, gives total of 65 for etoposide
length(intersect(methotrexate_most_min_genes, methotrexate_least_1se_genes))
#1, gives total of 6 genes for methotrexate

cisplatin_genes <- c(cisplatin_most_1se_genes, cisplatin_least_min_genes)
cisplatin_genes <- unique(cisplatin_genes) #103

etoposide_genes <- c(etoposide_most_min_genes, etoposide_least_min_genes)
etoposide_genes <- unique(etoposide_genes) #65

methotrexate_genes <- c(methotrexate_most_min_genes, methotrexate_least_1se_genes)
methotrexate_genes <- unique(methotrexate_genes) #6

intersect(cisplatin_genes, etoposide_genes)
#2 "ENSG00000102882" "ENSG00000172716"
intersect(cisplatin_genes, methotrexate_genes)
#0
intersect(etoposide_genes, methotrexate_genes)
#0

# gene names converted in g:profiler
cisplatin_gene_names <- read.csv('cisplatin_gene_names.csv', stringsAsFactors = FALSE, header = TRUE)
cisplatin_gene_names <- cisplatin_gene_names$name
etoposide_gene_names <- read.csv('etoposide_gene_names.csv', stringsAsFactors = FALSE, header = TRUE)
etoposide_gene_names <- etoposide_gene_names$name
methotrexate_gene_names <- read.csv('methotrexate_gene_names.csv', stringsAsFactors = FALSE, header = TRUE)
methotrexate_gene_names <- methotrexate_gene_names$name

cisp_etop <- intersect(cisplatin_gene_names, etoposide_gene_names)
#2 "MAPK3" "SLFN11"
cisp_meth <- intersect(cisplatin_gene_names, methotrexate_gene_names)
#0
etop_meth <- intersect(etoposide_gene_names, methotrexate_gene_names)
#0

common_genes <- c(cisp_etop, cisp_meth, etop_meth)
common_genes <- unique(common_genes)

gene_list <- list(cisplatin_gene_names, etoposide_gene_names, methotrexate_gene_names)
gene_table1 <- sapply(gene_list, 'length<-', max(lengths(gene_list)))
colnames(gene_table) <- c('Cisplatin_Genes', 'Etoposide_Genes', 'Methotrexate_Genes')

ft <- flextable(as.data.frame(gene_table), cheight = 0.01)
ft <- color(ft, i = ~ Cisplatin_Genes %in% common_genes, j = 1, color = 'red')
ft <- color(ft, i = ~ Etoposide_Genes %in% common_genes, j = 2, color = 'red')
ft <- color(ft, i = ~ Methotrexate_Genes %in% common_genes, j = 3, color = 'red')
ft


## CCLE ----
## load models ----
carboplatin_ccle_most_fit_elnet           <- readRDS('GLM_Models/carboplatin_ccle_most_model.rds')
carboplatin_ccle_least_fit_elnet          <- readRDS('GLM_Models/carboplatin_ccle_least_model.rds')

cyclophosphamide_ccle_most_fit_elnet      <- readRDS('GLM_Models/cyclophosphamide_ccle_most_model.rds')
cyclophosphamide_ccle_least_fit_elnet     <- readRDS('GLM_Models/cyclophosphamide_ccle_least_model.rds')

docetaxel_ccle_most_fit_elnet             <- readRDS('GLM_Models/docetaxel_ccle_most_model.rds')
docetaxel_ccle_least_fit_elnet            <- readRDS('GLM_Models/docetaxel_ccle_least_model.rds')

fluorouracil_ccle_most_fit_elnet          <- readRDS('GLM_Models/fluorouracil_ccle_most_model.rds')
fluorouracil_ccle_least_fit_elnet         <- readRDS('GLM_Models/fluorouracil_ccle_least_model.rds')

gemcitabine_ccle_most_fit_elnet           <- readRDS('GLM_Models/gemcitabine_ccle_most_model.rds')
gemcitabine_ccle_least_fit_elnet          <- readRDS('GLM_Models/gemcitabine_ccle_least_model.rds')

### capture genes used in models ----
## CARBOPLATIN MOST SENSITIVE
carboplatin_ccle_most_min_tmp_coefs <- coef(carboplatin_ccle_most_fit_elnet, s = 'lambda.min')
carboplatin_ccle_most_model_min <- data.frame(name = carboplatin_ccle_most_min_tmp_coefs@Dimnames[[1]][carboplatin_ccle_most_min_tmp_coefs@i + 1], coefficient = carboplatin_ccle_most_min_tmp_coefs@x)
write.csv(carboplatin_ccle_most_model_min, file = 'GLM_Models/carboplatin_ccle_most_model_min.csv', row.names = FALSE)

carboplatin_ccle_most_1se_tmp_coefs <- coef(carboplatin_ccle_most_fit_elnet, s = 'lambda.1se')
carboplatin_ccle_most_model_1se <- data.frame(name = carboplatin_ccle_most_1se_tmp_coefs@Dimnames[[1]][carboplatin_ccle_most_1se_tmp_coefs@i + 1], coefficient = carboplatin_ccle_most_1se_tmp_coefs@x)
write.csv(carboplatin_ccle_most_model_1se, file = 'GLM_Models/carboplatin_ccle_most_model_1se.csv', row.names = FALSE)

## CARBOPLATIN LEAST SENSITIVE
carboplatin_ccle_least_min_tmp_coefs <- coef(carboplatin_ccle_least_fit_elnet, s = 'lambda.min')
carboplatin_ccle_least_model_min <- data.frame(name = carboplatin_ccle_least_min_tmp_coefs@Dimnames[[1]][carboplatin_ccle_least_min_tmp_coefs@i + 1], coefficient = carboplatin_ccle_least_min_tmp_coefs@x)
write.csv(carboplatin_ccle_least_model_min, file = 'GLM_Models/carboplatin_ccle_least_model_min.csv', row.names = FALSE)

carboplatin_ccle_least_1se_tmp_coefs <- coef(carboplatin_ccle_least_fit_elnet, s = 'lambda.1se')
carboplatin_ccle_least_model_1se <- data.frame(name = carboplatin_ccle_least_1se_tmp_coefs@Dimnames[[1]][carboplatin_ccle_least_1se_tmp_coefs@i + 1], coefficient = carboplatin_ccle_least_1se_tmp_coefs@x)
write.csv(carboplatin_ccle_least_model_1se, file = 'GLM_Models/carboplatin_ccle_least_model_1se.csv', row.names = FALSE)

## CYCLOPHOSPHAMIDE MOST SENSITIVE
cyclophosphamide_ccle_most_min_tmp_coefs <- coef(cyclophosphamide_ccle_most_fit_elnet, s = 'lambda.min')
cyclophosphamide_ccle_most_model_min <- data.frame(name = cyclophosphamide_ccle_most_min_tmp_coefs@Dimnames[[1]][cyclophosphamide_ccle_most_min_tmp_coefs@i + 1], coefficient = cyclophosphamide_ccle_most_min_tmp_coefs@x)
write.csv(cyclophosphamide_ccle_most_model_min, file = 'GLM_Models/cyclophosphamide_ccle_most_model_min.csv', row.names = FALSE)

cyclophosphamide_ccle_most_1se_tmp_coefs <- coef(cyclophosphamide_ccle_most_fit_elnet, s = 'lambda.1se')
cyclophosphamide_ccle_most_model_1se <- data.frame(name = cyclophosphamide_ccle_most_1se_tmp_coefs@Dimnames[[1]][cyclophosphamide_ccle_most_1se_tmp_coefs@i + 1], coefficient = cyclophosphamide_ccle_most_1se_tmp_coefs@x)
write.csv(cyclophosphamide_ccle_most_model_1se, file = 'GLM_Models/cyclophosphamide_ccle_most_model_1se.csv', row.names = FALSE)

## CYCLOPHOSPHAMIDE LEAST SENSITIVE
cyclophosphamide_ccle_least_min_tmp_coefs <- coef(cyclophosphamide_ccle_least_fit_elnet, s = 'lambda.min')
cyclophosphamide_ccle_least_model_min <- data.frame(name = cyclophosphamide_ccle_least_min_tmp_coefs@Dimnames[[1]][cyclophosphamide_ccle_least_min_tmp_coefs@i + 1], coefficient = cyclophosphamide_ccle_least_min_tmp_coefs@x)
write.csv(cyclophosphamide_ccle_least_model_min, file = 'GLM_Models/cyclophosphamide_ccle_least_model_min.csv', row.names = FALSE)

cyclophosphamide_ccle_least_1se_tmp_coefs <- coef(cyclophosphamide_ccle_least_fit_elnet, s = 'lambda.1se')
cyclophosphamide_ccle_least_model_1se <- data.frame(name = cyclophosphamide_ccle_least_1se_tmp_coefs@Dimnames[[1]][cyclophosphamide_ccle_least_1se_tmp_coefs@i + 1], coefficient = cyclophosphamide_ccle_least_1se_tmp_coefs@x)
write.csv(cyclophosphamide_ccle_least_model_1se, file = 'GLM_Models/cyclophosphamide_ccle_least_model_1se.csv', row.names = FALSE)

## DOCETAXEL MOST SENSITIVE
docetaxel_ccle_most_min_tmp_coefs <- coef(docetaxel_ccle_most_fit_elnet, s = 'lambda.min')
docetaxel_ccle_most_model_min <- data.frame(name = docetaxel_ccle_most_min_tmp_coefs@Dimnames[[1]][docetaxel_ccle_most_min_tmp_coefs@i + 1], coefficient = docetaxel_ccle_most_min_tmp_coefs@x)
write.csv(docetaxel_ccle_most_model_min, file = 'GLM_Models/docetaxel_ccle_most_model_min.csv', row.names = FALSE)

docetaxel_ccle_most_1se_tmp_coefs <- coef(docetaxel_ccle_most_fit_elnet, s = 'lambda.1se')
docetaxel_ccle_most_model_1se <- data.frame(name = docetaxel_ccle_most_1se_tmp_coefs@Dimnames[[1]][docetaxel_ccle_most_1se_tmp_coefs@i + 1], coefficient = docetaxel_ccle_most_1se_tmp_coefs@x)
write.csv(docetaxel_ccle_most_model_1se, file = 'GLM_Models/docetaxel_ccle_most_model_1se.csv', row.names = FALSE)

## DOCETAXEL LEAST SENSITIVE
docetaxel_ccle_least_min_tmp_coefs <- coef(docetaxel_ccle_least_fit_elnet, s = 'lambda.min')
docetaxel_ccle_least_model_min <- data.frame(name = docetaxel_ccle_least_min_tmp_coefs@Dimnames[[1]][docetaxel_ccle_least_min_tmp_coefs@i + 1], coefficient = docetaxel_ccle_least_min_tmp_coefs@x)
write.csv(docetaxel_ccle_least_model_min, file = 'GLM_Models/docetaxel_ccle_least_model_min.csv', row.names = FALSE)

docetaxel_ccle_least_1se_tmp_coefs <- coef(docetaxel_ccle_least_fit_elnet, s = 'lambda.1se')
docetaxel_ccle_least_model_1se <- data.frame(name = docetaxel_ccle_least_1se_tmp_coefs@Dimnames[[1]][docetaxel_ccle_least_1se_tmp_coefs@i + 1], coefficient = docetaxel_ccle_least_1se_tmp_coefs@x)
write.csv(docetaxel_ccle_least_model_1se, file = 'GLM_Models/docetaxel_ccle_least_model_1se.csv', row.names = FALSE)

## FLUOROURACIL MOST SENSITIVE
fluorouracil_ccle_most_min_tmp_coefs <- coef(fluorouracil_ccle_most_fit_elnet, s = 'lambda.min')
fluorouracil_ccle_most_model_min <- data.frame(name = fluorouracil_ccle_most_min_tmp_coefs@Dimnames[[1]][fluorouracil_ccle_most_min_tmp_coefs@i + 1], coefficient = fluorouracil_ccle_most_min_tmp_coefs@x)
write.csv(fluorouracil_ccle_most_model_min, file = 'GLM_Models/fluorouracil_ccle_most_model_min.csv', row.names = FALSE)

fluorouracil_ccle_most_1se_tmp_coefs <- coef(fluorouracil_ccle_most_fit_elnet, s = 'lambda.1se')
fluorouracil_ccle_most_model_1se <- data.frame(name = fluorouracil_ccle_most_1se_tmp_coefs@Dimnames[[1]][fluorouracil_ccle_most_1se_tmp_coefs@i + 1], coefficient = fluorouracil_ccle_most_1se_tmp_coefs@x)
write.csv(fluorouracil_ccle_most_model_1se, file = 'GLM_Models/fluorouracil_ccle_most_model_1se.csv', row.names = FALSE)

## FLUOROURACIL LEAST SENSITIVE
fluorouracil_ccle_least_min_tmp_coefs <- coef(fluorouracil_ccle_least_fit_elnet, s = 'lambda.min')
fluorouracil_ccle_least_model_min <- data.frame(name = fluorouracil_ccle_least_min_tmp_coefs@Dimnames[[1]][fluorouracil_ccle_least_min_tmp_coefs@i + 1], coefficient = fluorouracil_ccle_least_min_tmp_coefs@x)
write.csv(fluorouracil_ccle_least_model_min, file = 'GLM_Models/fluorouracil_ccle_least_model_min.csv', row.names = FALSE)

fluorouracil_ccle_least_1se_tmp_coefs <- coef(fluorouracil_ccle_least_fit_elnet, s = 'lambda.1se')
fluorouracil_ccle_least_model_1se <- data.frame(name = fluorouracil_ccle_least_1se_tmp_coefs@Dimnames[[1]][fluorouracil_ccle_least_1se_tmp_coefs@i + 1], coefficient = fluorouracil_ccle_least_1se_tmp_coefs@x)
write.csv(fluorouracil_ccle_least_model_1se, file = 'GLM_Models/fluorouracil_ccle_least_model_1se.csv', row.names = FALSE)

## GEMCITABINE MOST SENSITIVE
gemcitabine_ccle_most_min_tmp_coefs <- coef(gemcitabine_ccle_most_fit_elnet, s = 'lambda.min')
gemcitabine_ccle_most_model_min <- data.frame(name = gemcitabine_ccle_most_min_tmp_coefs@Dimnames[[1]][gemcitabine_ccle_most_min_tmp_coefs@i + 1], coefficient = gemcitabine_ccle_most_min_tmp_coefs@x)
write.csv(gemcitabine_ccle_most_model_min, file = 'GLM_Models/gemcitabine_ccle_most_model_min.csv', row.names = FALSE)

gemcitabine_ccle_most_1se_tmp_coefs <- coef(gemcitabine_ccle_most_fit_elnet, s = 'lambda.1se')
gemcitabine_ccle_most_model_1se <- data.frame(name = gemcitabine_ccle_most_1se_tmp_coefs@Dimnames[[1]][gemcitabine_ccle_most_1se_tmp_coefs@i + 1], coefficient = gemcitabine_ccle_most_1se_tmp_coefs@x)
write.csv(gemcitabine_ccle_most_model_1se, file = 'GLM_Models/gemcitabine_ccle_most_model_1se.csv', row.names = FALSE)

## GEMCITABINE LEAST SENSITIVE
gemcitabine_ccle_least_min_tmp_coefs <- coef(gemcitabine_ccle_least_fit_elnet, s = 'lambda.min')
gemcitabine_ccle_least_model_min <- data.frame(name = gemcitabine_ccle_least_min_tmp_coefs@Dimnames[[1]][gemcitabine_ccle_least_min_tmp_coefs@i + 1], coefficient = gemcitabine_ccle_least_min_tmp_coefs@x)
write.csv(gemcitabine_ccle_least_model_min, file = 'GLM_Models/gemcitabine_ccle_least_model_min.csv', row.names = FALSE)

gemcitabine_ccle_least_1se_tmp_coefs <- coef(gemcitabine_ccle_least_fit_elnet, s = 'lambda.1se')
gemcitabine_ccle_least_model_1se <- data.frame(name = gemcitabine_ccle_least_1se_tmp_coefs@Dimnames[[1]][gemcitabine_ccle_least_1se_tmp_coefs@i + 1], coefficient = gemcitabine_ccle_least_1se_tmp_coefs@x)
write.csv(gemcitabine_ccle_least_model_1se, file = 'GLM_Models/gemcitabine_ccle_least_model_1se.csv', row.names = FALSE)

## PUT TOGETHER IN TABLE ----
carboplatin_ccle_most_min_genes <- as.character(carboplatin_ccle_most_model_min$name)
carboplatin_ccle_most_min_genes <- carboplatin_ccle_most_min_genes[-1]

carboplatin_ccle_most_1se_genes <- as.character(carboplatin_ccle_most_model_1se$name)
carboplatin_ccle_most_1se_genes <- carboplatin_ccle_most_1se_genes[-1]

carboplatin_ccle_least_min_genes <- as.character(carboplatin_ccle_least_model_min$name)
carboplatin_ccle_least_min_genes <- carboplatin_ccle_least_min_genes[-1]

carboplatin_ccle_least_1se_genes <- as.character(carboplatin_ccle_least_model_1se$name)
carboplatin_ccle_least_1se_genes <- carboplatin_ccle_least_1se_genes[-1]

cyclophosphamide_ccle_most_min_genes <- as.character(cyclophosphamide_ccle_most_model_min$name)
cyclophosphamide_ccle_most_min_genes <- cyclophosphamide_ccle_most_min_genes[-1]

cyclophosphamide_ccle_most_1se_genes <- as.character(cyclophosphamide_ccle_most_model_1se$name)
cyclophosphamide_ccle_most_1se_genes <- cyclophosphamide_ccle_most_1se_genes[-1]

cyclophosphamide_ccle_least_min_genes <- as.character(cyclophosphamide_ccle_least_model_min$name)
cyclophosphamide_ccle_least_min_genes <- cyclophosphamide_ccle_least_min_genes[-1]

cyclophosphamide_ccle_least_1se_genes <- as.character(cyclophosphamide_ccle_least_model_1se$name)
cyclophosphamide_ccle_least_1se_genes <- cyclophosphamide_ccle_least_1se_genes[-1]

docetaxel_ccle_most_min_genes <- as.character(docetaxel_ccle_most_model_min$name)
docetaxel_ccle_most_min_genes <- docetaxel_ccle_most_min_genes[-1]

docetaxel_ccle_most_1se_genes <- as.character(docetaxel_ccle_most_model_1se$name)
docetaxel_ccle_most_1se_genes <- docetaxel_ccle_most_1se_genes[-1]

docetaxel_ccle_least_min_genes <- as.character(docetaxel_ccle_least_model_min$name)
docetaxel_ccle_least_min_genes <- docetaxel_ccle_least_min_genes[-1]

docetaxel_ccle_least_1se_genes <- as.character(docetaxel_ccle_least_model_1se$name)
docetaxel_ccle_least_1se_genes <- docetaxel_ccle_least_1se_genes[-1]

fluorouracil_ccle_most_min_genes <- as.character(fluorouracil_ccle_most_model_min$name)
fluorouracil_ccle_most_min_genes <- fluorouracil_ccle_most_min_genes[-1]

fluorouracil_ccle_most_1se_genes <- as.character(fluorouracil_ccle_most_model_1se$name)
fluorouracil_ccle_most_1se_genes <- fluorouracil_ccle_most_1se_genes[-1]

fluorouracil_ccle_least_min_genes <- as.character(fluorouracil_ccle_least_model_min$name)
fluorouracil_ccle_least_min_genes <- fluorouracil_ccle_least_min_genes[-1]

fluorouracil_ccle_least_1se_genes <- as.character(fluorouracil_ccle_least_model_1se$name)
fluorouracil_ccle_least_1se_genes <- fluorouracil_ccle_least_1se_genes[-1]

gemcitabine_ccle_most_min_genes <- as.character(gemcitabine_ccle_most_model_min$name)
gemcitabine_ccle_most_min_genes <- gemcitabine_ccle_most_min_genes[-1]

gemcitabine_ccle_most_1se_genes <- as.character(gemcitabine_ccle_most_model_1se$name)
gemcitabine_ccle_most_1se_genes <- gemcitabine_ccle_most_1se_genes[-1]

gemcitabine_ccle_least_min_genes <- as.character(gemcitabine_ccle_least_model_min$name)
gemcitabine_ccle_least_min_genes <- gemcitabine_ccle_least_min_genes[-1]

gemcitabine_ccle_least_1se_genes <- as.character(gemcitabine_ccle_least_model_1se$name)
gemcitabine_ccle_least_1se_genes <- gemcitabine_ccle_least_1se_genes[-1]


# using: carboplatin most 1se, carboplatin least 1se, cyclophosphamide most 1se, cyclophosphamide least 1se, 
# docetaxel most 1se, docetaxel least min, etoposide most 1se, etoposide least 1se, 
# fluorouracil most 1se, fluorouracil least 1se, gemcitabine most 1se, gemcitabine least 1se, 
# methotrexate most 1se, methotrexate least 1se, oxalaplatin most min, oxalaplatin least 1se

length(intersect(carboplatin_ccle_most_1se_genes, carboplatin_ccle_least_1se_genes))
#0, gives total of 4 for carboplatin
length(intersect(cyclophosphamide_ccle_most_1se_genes, cyclophosphamide_least_1se_genes)) # no cyclo least models
#0, gives total of 2 for cyclophosphamide
length(intersect(docetaxel_ccle_most_1se_genes, docetaxel_ccle_least_min_genes))
#0, gives total of 45 for docetaxel
length(intersect(fluorouracil_ccle_most_1se_genes, fluorouracil_ccle_least_1se_genes))
#0, gives a total of 9 genes for fluorouracil
length(intersect(gemcitabine_ccle_most_1se_genes, gemcitabine_ccle_least_1se_genes))
#0, gives a total of 16 genes for gemcitabine

carboplatin_ccle_genes <- c(carboplatin_ccle_most_1se_genes, carboplatin_ccle_least_1se_genes)
carboplatin_ccle_genes <- unique(carboplatin_ccle_genes) #4

cyclophosphamide_ccle_genes <- c(cyclophosphamide_ccle_most_1se_genes, cyclophosphamide_ccle_least_1se_genes)
cyclophosphamide_ccle_genes <- unique(cyclophosphamide_ccle_genes) #2

docetaxel_ccle_genes <- c(docetaxel_ccle_most_1se_genes, docetaxel_ccle_least_min_genes)
docetaxel_ccle_genes <- unique(docetaxel_ccle_genes) #45

fluorouracil_ccle_genes <- c(fluorouracil_ccle_most_1se_genes, fluorouracil_ccle_least_1se_genes)
fluorouracil_ccle_genes <- unique(fluorouracil_ccle_genes) #9

gemcitabine_ccle_genes <- c(gemcitabine_ccle_most_1se_genes, gemcitabine_ccle_least_1se_genes)
gemcitabine_ccle_genes <- unique(gemcitabine_ccle_genes) #16


length(intersect(carboplatin_ccle_genes, cyclophosphamide_ccle_genes))
#0
length(intersect(carboplatin_ccle_genes, docetaxel_ccle_genes))
#0
length(intersect(carboplatin_ccle_genes, fluorouracil_ccle_genes))
#0
length(intersect(carboplatin_ccle_genes, gemcitabine_ccle_genes))
#0
length(intersect(cyclophosphamide_ccle_genes, docetaxel_ccle_genes))
#0
length(intersect(cyclophosphamide_ccle_genes, fluorouracil_ccle_genes))
#0
length(intersect(cyclophosphamide_ccle_genes, gemcitabine_ccle_genes))
#0
length(intersect(docetaxel_ccle_genes, fluorouracil_ccle_genes))
#0
length(intersect(docetaxel_ccle_genes, gemcitabine_ccle_genes))
#0
length(intersect(fluorouracil_ccle_genes, gemcitabine_ccle_genes))
#0

gene_list <- list(carboplatin_ccle_genes, cyclophosphamide_ccle_genes, docetaxel_ccle_genes, 
                  fluorouracil_ccle_genes, gemcitabine_ccle_genes)
gene_table <- sapply(gene_list, 'length<-', max(lengths(gene_list)))
colnames(gene_table) <- c('Carboplatin_Genes', 'Cyclophosphamide_Genes', 'Docetaxel_Genes', 
                          'Fluorouracil_Genes', 'Gemcitabine_Genes')
write.csv(gene_table, file = 'CCLE_model_genes.csv', row.names = FALSE)

# gene names converted in g:profiler
carboplatin_ccle_gene_names <- read.csv('carboplatin_ccle_genes.csv', stringsAsFactors = FALSE, header = TRUE)
carboplatin_ccle_gene_names <- carboplatin_ccle_gene_names$name
cyclophosphamide_ccle_gene_names <- read.csv('cyclophosphamide_ccle_genes.csv', stringsAsFactors = FALSE, header = TRUE)
cyclophosphamide_ccle_gene_names <- cyclophosphamide_ccle_gene_names$name
docetaxel_ccle_gene_names <- read.csv('docetaxel_ccle_genes.csv', stringsAsFactors = FALSE, header = TRUE)
docetaxel_ccle_gene_names <- docetaxel_ccle_gene_names$name
fluorouracil_ccle_gene_names <- read.csv('fluorouracil_ccle_genes.csv', stringsAsFactors = FALSE, header = TRUE)
fluorouracil_ccle_gene_names <- fluorouracil_ccle_gene_names$name
gemcitabine_ccle_gene_names <- read.csv('gemcitabine_ccle_genes.csv', stringsAsFactors = FALSE, header = TRUE)
gemcitabine_ccle_gene_names <- gemcitabine_ccle_gene_names$name


length(intersect(carboplatin_ccle_gene_names, cyclophosphamide_ccle_gene_names))
#0
length(intersect(carboplatin_ccle_gene_names, docetaxel_ccle_gene_names))
#0
length(intersect(carboplatin_ccle_gene_names, fluorouracil_ccle_gene_names))
#0
length(intersect(carboplatin_ccle_gene_names, gemcitabine_ccle_gene_names))
#0
length(intersect(cyclophosphamide_ccle_gene_names, docetaxel_ccle_gene_names))
#0
length(intersect(cyclophosphamide_ccle_gene_names, fluorouracil_ccle_gene_names))
#0
length(intersect(cyclophosphamide_ccle_gene_names, gemcitabine_ccle_gene_names))
#0
length(intersect(docetaxel_ccle_gene_names, fluorouracil_ccle_gene_names))
#0
length(intersect(docetaxel_ccle_gene_names, gemcitabine_ccle_gene_names))
#0
length(intersect(fluorouracil_ccle_gene_names, gemcitabine_ccle_gene_names))
#0

gene_list <- list(carboplatin_ccle_gene_names, cyclophosphamide_ccle_gene_names, docetaxel_ccle_gene_names, fluorouracil_ccle_gene_names, gemcitabine_ccle_gene_names)
gene_table <- sapply(gene_list, 'length<-', max(lengths(gene_list)))
colnames(gene_table) <- c('Carboplatin_Genes', 'Cyclophosphamide_Genes', 'Docetaxel_Genes', 'Fluorouracil_Genes', 'Gemcitabine_Genes')
flextable(as.data.frame(gene_table))

### treemap of all gene set numbers ----
drug_labels <- c('cisplatin (GDSC)\nn = 103', 'etoposide (GDSC)\nn = 65', 'methotrexate (GDSC)\nn = 6', 'carboplatin (CCLE)\nn = 4', 
                 'cyclophosphamide (CCLE)\nn = 2', 'docetaxel (CCLE)\nn = 45', 'fluorouracil (CCLE)\nn = 9', 'gemcitabine (CCLE)\nn = 16')
drug_numbers <- c(103, 65, 6, 4, 2, 45, 9, 16)
tree_df <- data.frame(drug_labels, drug_numbers)

png(filename = 'Images/geneset_treemap.png', width = 1000)
treemap(dtf = tree_df, vSize = 'drug_numbers', index = 'drug_labels', fontsize.labels = 10, palette = 'Dark2', title = 'Number of Genes per Drug Model (N = 248)', fontcolor.labels = 'black')
dev.off()

drug_labels <- c(rep('cisplatin (GDSC)\nn = 103', 2), rep('etoposide (GDSC)\nn = 65', 2), rep('methotrexate (GDSC)\nn = 6', 2), rep('carboplatin (CCLE)\nn = 4', 2), 
                 rep('cyclophosphamide (CCLE)\nn = 2', 2), rep('docetaxel (CCLE)\nn = 45', 2), rep('fluorouracil (CCLE)\nn = 9', 2), rep('gemcitabine (CCLE)\nn = 16', 2))
model_types <- rep(c('sensitive', 'resistant'), 8)
drug_numbers <- c(19, 88, 8, 62, 4, 3, 1, 3, 2, 0, 18, 27, 4, 5, 1, 15)
tree_df <- data.frame(drug_labels, model_types, drug_numbers)
treemap(dtf = tree_df, vSize = 'drug_numbers', index = c('drug_labels', 'model_types'), fontface.labels = c(2,1), fontsize.labels = c(15,12), bg.labels = 'transparent', overlap.labels = 0.5, align.labels = list(c('center', 'center'), c('right', 'bottom')), palette = 'Dark2', title = 'Number of Genes per Drug Model (N = 247)', fontcolor.labels = 'black')
