## THIS SCRIPT HANDLES AND PROCESSES ALL CLINICAL DATA NEEDED

### load necessary packages ----
if (!require ('openxlsx')) install.packages('openxlsx')
library(openxlsx) # to read in excel files
if (!require ('dplyr')) install.packages('dplyr')
library(dplyr) # to form some of the files

### functions needed ----
# function to read in GDSC clinical data
read_gdsc_clinical <- function(sheet_number) {
  drug_data <- read.xlsx(gdsc_files, sheet_number)
  return(drug_data[, c(3:6, 12)])
}

# function to read in CCLE clinical data
read_ccle_clinical <- function(sheet_number) {
  drug_data <- read.xlsx(ccle_files, sheet_number)
  drug_data <- drug_data[, c(4:6)]
  drug_data_new <- drug_data %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
  drug_data <- merge(drug_data, drug_data_new, by = 'Cell.line.name')
  drug_data <- drug_data[, -3]
  return(unique(drug_data))
}

# function to read in TCGA clinical data
read_tcga_clinical <- function(tcga_file_number) {
  tcga_phenos <- read.delim(tcga_files[tcga_file_number], sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  tcga_phenos_short <- tcga_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 
                                       'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
  tcga_phenos_not_lost <- tcga_phenos_short[tcga_phenos_short$lost_follow_up != 'YES', ]
  tcga_phenos_w_drug <- tcga_phenos_not_lost[tcga_phenos_not_lost$drug_name != '', ]
  tcga_phenos_w_drug$OS <- ifelse(tcga_phenos_w_drug$vital_status.diagnoses == 'dead', tcga_phenos_w_drug[, 5], tcga_phenos_w_drug[, 6])
  tcga_phenos_w_drug$PFS <- ifelse(is.na(tcga_phenos_w_drug[, 2]), tcga_phenos_w_drug$OS, tcga_phenos_w_drug[, 2])
  return(tcga_phenos_w_drug)
}

# functions for dummy variables in GDSC data
# most senistive
most_sens_bin_gdsc <- function(drug_data) {
  ifelse(drug_data$LN_IC50 < quantile(drug_data$LN_IC50, probs = 0.20, na.rm = TRUE), 1, 0)
}
# least sensitive
least_sens_bin_gdsc <- function(drug_data) {
  ifelse(drug_data$LN_IC50 > quantile(drug_data$LN_IC50, probs = 0.80, na.rm = TRUE), 1, 0)
}

# functions for dummy variables in CCLE data
# most sensitive
most_sens_bin_ccle <- function(drug_data) {
  ifelse(drug_data$AUC < quantile(drug_data$AUC, probs = 0.20, na.rm = TRUE), 1, 0)
}
# least sensitive
least_sens_bin_ccle <- function(drug_data) {
  ifelse(drug_data$AUC > quantile(drug_data$AUC, probs = 0.80, na.rm = TRUE), 1, 0)
}

# functions for dummy variables in TCGA data
# most sensitive
most_sens_bin_tcga <- function(drug_data) {
  ifelse(drug_data$PFS < quantile(drug_data$PFS, probs = 0.20, na.rm = TRUE), 1, 0)
}
# least sensitive
least_sens_bin_tcga <- function(drug_data) {
  ifelse(drug_data$PFS > quantile(drug_data$PFS, probs = 0.80, na.rm = TRUE), 1, 0)
}

### import data -----

# GDSC (ln(IC50))
gdsc_files <- 'Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx'

# bring in clinical data for drugs we're interested in
cisplatin     <- read_gdsc_clinical(4)

etoposide     <- read_gdsc_clinical(7)

methotrexate  <- read_gdsc_clinical(9)


#CCLE (AUC)
ccle_files <- 'Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx'

#bring in clinical data for drugs we're interested in
carboplatin_ccle      <- read_ccle_clinical(3)

cyclophosphamide_ccle <- read_ccle_clinical(6)

docetaxel_ccle        <- read_ccle_clinical(11)

fluorouracil_ccle     <- read_ccle_clinical(14)

gemcitabine_ccle      <- read_ccle_clinical(15)


#TCGA (RFS)
tcga_files <- list.files(path = "TCGA_files/", pattern = "*phenotype.tsv", full.names = T)

#process TCGA classes we're interested in
blca_clinical <- read_tcga_clinical(2)

brca_clinical <- read_tcga_clinical(3)

cesc_clinical <- read_tcga_clinical(4)

chol_clinical <- read_tcga_clinical(5)

coad_clinical <- read_tcga_clinical(6)

dlbc_clinical <- read_tcga_clinical(7)

esca_clinical <- read_tcga_clinical(8)

hnsc_clinical <- read_tcga_clinical(10)

kirc_clinical <- read_tcga_clinical(12)

lihc_clinical <- read_tcga_clinical(15)

luad_clinical <- read_tcga_clinical(16)

lusc_clinical <- read_tcga_clinical(17)

meso_clinical <- read_tcga_clinical(18)

ov_clinical   <- read_tcga_clinical(19)

paad_clinical <- read_tcga_clinical(20)

read_clinical <- read_tcga_clinical(23)

sarc_clinical <- read_tcga_clinical(24)

skcm_clinical <- read_tcga_clinical(25)

stad_clinical <- read_tcga_clinical(26)

tgct_clinical <- read_tcga_clinical(27)

ucec_clinical <- read_tcga_clinical(30)

ucs_clinical  <- read_tcga_clinical(31)


### create dummy binary variables for most/least sensitive -------

# GDSC
cisplatin$most_sensitive        <- most_sens_bin_gdsc(cisplatin)
cisplatin$least_sensitive       <- least_sens_bin_gdsc(cisplatin)

etoposide$most_sensitive        <- most_sens_bin_gdsc(etoposide)
etoposide$least_sensitive       <- least_sens_bin_gdsc(etoposide)

methotrexate$most_sensitive     <- most_sens_bin_gdsc(methotrexate)
methotrexate$least_sensitive    <- least_sens_bin_gdsc(methotrexate)

## CCLE
carboplatin_ccle$most_sensitive           <- most_sens_bin_ccle(carboplatin_ccle)
carboplatin_ccle$least_sensitive          <- least_sens_bin_ccle(carboplatin_ccle)

cyclophosphamide_ccle$most_sensitive      <- most_sens_bin_ccle(cyclophosphamide_ccle)
cyclophosphamide_ccle$least_sensitive     <- least_sens_bin_ccle(cyclophosphamide_ccle)

docetaxel_ccle$most_sensitive             <- most_sens_bin_ccle(docetaxel_ccle)
docetaxel_ccle$least_sensitive            <- least_sens_bin_ccle(docetaxel_ccle)

fluorouracil_ccle$most_sensitive          <- most_sens_bin_ccle(fluorouracil_ccle)
fluorouracil_ccle$least_sensitive         <- least_sens_bin_ccle(fluorouracil_ccle)

gemcitabine_ccle$most_sensitive           <- most_sens_bin_ccle(gemcitabine_ccle)
gemcitabine_ccle$least_sensitive          <- least_sens_bin_ccle(gemcitabine_ccle)

## TCGA
blca_clinical$most_sensitive        <- most_sens_bin_tcga(blca_clinical)
blca_clinical$least_sensitive       <- least_sens_bin_tcga(blca_clinical)

brca_clinical$most_sensitive        <- most_sens_bin_tcga(brca_clinical)
brca_clinical$least_sensitive       <- least_sens_bin_tcga(brca_clinical)

cesc_clinical$most_sensitive        <- most_sens_bin_tcga(cesc_clinical)
cesc_clinical$least_sensitive       <- least_sens_bin_tcga(cesc_clinical)

chol_clinical$most_sensitive        <- most_sens_bin_tcga(chol_clinical)
chol_clinical$least_sensitive       <- least_sens_bin_tcga(chol_clinical)

coad_clinical$most_sensitive        <- most_sens_bin_tcga(coad_clinical)
coad_clinical$least_sensitive       <- least_sens_bin_tcga(coad_clinical)

dlbc_clinical$most_sensitive        <- most_sens_bin_tcga(dlbc_clinical)
dlbc_clinical$least_sensitive       <- least_sens_bin_tcga(dlbc_clinical)

esca_clinical$most_sensitive        <- most_sens_bin_tcga(esca_clinical)
esca_clinical$least_sensitive       <- least_sens_bin_tcga(esca_clinical)

hnsc_clinical$most_sensitive        <- most_sens_bin_tcga(hnsc_clinical)
hnsc_clinical$least_sensitive       <- least_sens_bin_tcga(hnsc_clinical)

kirc_clinical$most_sensitive        <- most_sens_bin_tcga(kirc_clinical)
kirc_clinical$least_sensitive       <- least_sens_bin_tcga(kirc_clinical)

lihc_clinical$most_sensitive        <- most_sens_bin_tcga(lihc_clinical)
lihc_clinical$least_sensitive       <- least_sens_bin_tcga(lihc_clinical)

luad_clinical$most_sensitive        <- most_sens_bin_tcga(luad_clinical)
luad_clinical$least_sensitive       <- least_sens_bin_tcga(luad_clinical)

lusc_clinical$most_sensitive        <- most_sens_bin_tcga(lusc_clinical)
lusc_clinical$least_sensitive       <- least_sens_bin_tcga(lusc_clinical)

meso_clinical$most_sensitive        <- most_sens_bin_tcga(meso_clinical)
meso_clinical$least_sensitive       <- least_sens_bin_tcga(meso_clinical)

ov_clinical$most_sensitive          <- most_sens_bin_tcga(ov_clinical)
ov_clinical$least_sensitive         <- most_sens_bin_tcga(ov_clinical)

paad_clinical$most_sensitive        <- most_sens_bin_tcga(paad_clinical)
paad_clinical$least_sensitive       <- least_sens_bin_tcga(paad_clinical)

read_clinical$most_sensitive        <- most_sens_bin_tcga(read_clinical)
read_clinical$least_sensitive       <- least_sens_bin_tcga(read_clinical)

sarc_clinical$most_sensitive        <- most_sens_bin_tcga(sarc_clinical)
sarc_clinical$least_sensitive       <- least_sens_bin_tcga(sarc_clinical)

skcm_clinical$most_sensitive        <- most_sens_bin_tcga(skcm_clinical)
skcm_clinical$least_sensitive       <- least_sens_bin_tcga(skcm_clinical)

stad_clinical$most_sensitive        <- most_sens_bin_tcga(stad_clinical)
stad_clinical$least_sensitive       <- least_sens_bin_tcga(stad_clinical)

tgct_clinical$most_sensitive        <- most_sens_bin_tcga(tgct_clinical)
tgct_clinical$least_sensitive       <- least_sens_bin_tcga(tgct_clinical)

ucec_clinical$most_sensitive        <- most_sens_bin_tcga(ucec_clinical)
ucec_clinical$least_sensitive       <- least_sens_bin_tcga(ucec_clinical)

ucs_clinical$most_sensitive         <- most_sens_bin_tcga(ucs_clinical)
ucs_clinical$least_sensitive        <- least_sens_bin_tcga(ucs_clinical)


## write to files for later ----
# GDSC
write.csv(cisplatin, file = 'Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv')
write.csv(etoposide, file = 'Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv')
write.csv(methotrexate, file = 'Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv')

# CCLE
write.csv(carboplatin_ccle, file = 'Processed_Clinical_Data/carboplatin_ccle_clinical_processed.csv')
write.csv(cyclophosphamide_ccle, file = 'Processed_Clinical_Data/cyclophosphamide_ccle_clinical_processed.csv')
write.csv(docetaxel_ccle, file = 'Processed_Clinical_Data/docetaxel_ccle_clinical_processed.csv')
write.csv(fluorouracil_ccle, file = 'Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv')
write.csv(gemcitabine_ccle, file = 'Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv')

# TCGA
write.csv(blca_phenos_done, file = 'Processed_Clinical_Data/blca_tcga_clinical_processed.csv')
write.csv(brca_phenos_done, file = 'Processed_Clinical_Data/brca_tcga_clinical_processed.csv')
write.csv(cesc_phenos_done, file = 'Processed_Clinical_Data/cesc_tcga_clinical_processed.csv')
write.csv(chol_phenos_done, file = 'Processed_Clinical_Data/chol_tcga_clinical_processed.csv')
write.csv(coad_phenos_done, file = 'Processed_Clinical_Data/coad_tcga_clinical_processed.csv')
write.csv(dlbc_phenos_done, file = 'Processed_Clinical_Data/dlbc_tcga_clinical_processed.csv')
write.csv(esca_phenos_done, file = 'Processed_Clinical_Data/esca_tcga_clinical_processed.csv')
write.csv(hnsc_phenos_done, file = 'Processed_Clinical_Data/hnsc_tcga_clinical_processed.csv')
write.csv(kirc_phenos_done, file = 'Processed_Clinical_Data/kirc_tcga_clinical_processed.csv')
write.csv(lihc_phenos_done, file = 'Processed_Clinical_Data/lihc_tcga_clinical_processed.csv')
write.csv(luad_phenos_done, file = 'Processed_Clinical_Data/luad_tcga_clinical_processed.csv')
write.csv(lusc_phenos_done, file = 'Processed_Clinical_Data/lusc_tcga_clinical_processed.csv')
write.csv(meso_phenos_done, file = 'Processed_Clinical_Data/meso_tcga_clinical_processed.csv')
write.csv(ov_phenos_done, file = 'Processed_Clinical_Data/ov_tcga_clinical_processed.csv')
write.csv(paad_phenos_done, file = 'Processed_Clinical_Data/paad_tcga_clinical_processed.csv')
write.csv(read_phenos_done, file = 'Processed_Clinical_Data/read_tcga_clinical_processed.csv')
write.csv(sarc_phenos_done, file = 'Processed_Clinical_Data/sarc_tcga_clinical_processed.csv')
write.csv(skcm_phenos_done, file = 'Processed_Clinical_Data/skcm_tcga_clinical_processed.csv')
write.csv(stad_phenos_done, file = 'Processed_Clinical_Data/stad_tcga_clinical_processed.csv')
write.csv(tgct_phenos_done, file = 'Processed_Clinical_Data/tgct_tcga_clinical_processed.csv')
write.csv(ucec_phenos_done, file = 'Processed_Clinical_Data/ucec_tcga_clinical_processed.csv')
write.csv(ucs_phenos_done, file = 'Processed_Clinical_Data/ucs_tcga_clinical_processed.csv')
