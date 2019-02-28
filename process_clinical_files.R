## THIS SCRIPT HANDLES AND PROCESSES ALL CLINICAL DATA NEEDED

### load necessary packages ----
if (!require ('openxlsx')) install.packages('openxlsx')
library(openxlsx) # to read in excel files

### import data -----
# GDSC (ln(IC50))

cisplatin     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 4)
# 680 x 14
cisplatin     <- cisplatin[, c(3:6, 12)]

etoposide     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 7)
# 718 x 14
etoposide     <- etoposide[, c(3:6, 12)]

gemcitabine   <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 8)
# 707 x 14
gemcitabine   <- gemcitabine[, c(3:6, 12)]

methotrexate  <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 9)
# 679 x 14
methotrexate  <- methotrexate[, c(3:6, 12)]

# cell line info for breaking up into TCGA categories
gdsc_cell_line_info <- read.csv('Clinical_Files/Cell_listTue Dec 18 19_10_36 2018 (version 1).csv', header = TRUE, stringsAsFactors = FALSE)
gdsc_cell_line_info <- gdsc_cell_line_info[, c(1,2,5,6)]
colnames(gdsc_cell_line_info)[c(3,4)] <- c('TCGA_class', 'Tissue_subtype')

#CCLE (AUC)
carboplatin_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 3)
carboplatin_ccle          <- carboplatin_ccle[, c(4:6)]
carboplatin_ccle_new      <- carboplatin_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
carboplatin_ccle          <- merge(carboplatin_ccle, carboplatin_ccle_new, by = 'Cell.line.name')
carboplatin_ccle          <- carboplatin_ccle[, -3]
carboplatin_ccle          <- unique(carboplatin_ccle)

cyclophosphamide_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 6)
cyclophosphamide_ccle          <- cyclophosphamide_ccle[, c(4:6)]
cyclophosphamide_ccle_new      <- cyclophosphamide_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
cyclophosphamide_ccle          <- merge(cyclophosphamide_ccle, cyclophosphamide_ccle_new, by = 'Cell.line.name')
cyclophosphamide_ccle          <- cyclophosphamide_ccle[, -3]
cyclophosphamide_ccle          <- unique(cyclophosphamide_ccle)

dacarbazine_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 9)
dacarbazine_ccle          <- dacarbazine_ccle[, c(4:6)]
dacarbazine_ccle_new      <- dacarbazine_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
dacarbazine_ccle          <- merge(dacarbazine_ccle, dacarbazine_ccle_new, by = 'Cell.line.name')
dacarbazine_ccle          <- dacarbazine_ccle[, -3]
dacarbazine_ccle          <- unique(dacarbazine_ccle)

docetaxel_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 11)
docetaxel_ccle          <- docetaxel_ccle[, c(4:6)]
docetaxel_ccle_new      <- docetaxel_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
docetaxel_ccle          <- merge(docetaxel_ccle, docetaxel_ccle_new, by = 'Cell.line.name')
docetaxel_ccle          <- docetaxel_ccle[, -3]
docetaxel_ccle          <- unique(docetaxel_ccle)

doxorubicin_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 12)
doxorubicin_ccle          <- doxorubicin_ccle[, c(4:6)]
doxorubicin_ccle_new      <- doxorubicin_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
doxorubicin_ccle          <- merge(doxorubicin_ccle, doxorubicin_ccle_new, by = 'Cell.line.name')
doxorubicin_ccle          <- doxorubicin_ccle[, -3]
doxorubicin_ccle          <- unique(doxorubicin_ccle)

etoposide_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 13)
etoposide_ccle          <- etoposide_ccle[, c(4:6)]
etoposide_ccle_new      <- etoposide_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
etoposide_ccle          <- merge(etoposide_ccle, etoposide_ccle_new, by = 'Cell.line.name')
etoposide_ccle          <- etoposide_ccle[, -3]
etoposide_ccle          <- unique(etoposide_ccle)

fluorouracil_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 14)
fluorouracil_ccle          <- fluorouracil_ccle[, c(4:6)]
fluorouracil_ccle_new      <- fluorouracil_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
fluorouracil_ccle          <- merge(fluorouracil_ccle, fluorouracil_ccle_new, by = 'Cell.line.name')
fluorouracil_ccle          <- fluorouracil_ccle[, -3]
fluorouracil_ccle          <- unique(fluorouracil_ccle)

gemcitabine_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 15)
gemcitabine_ccle          <- gemcitabine_ccle[, c(4:6)]
gemcitabine_ccle_new      <- gemcitabine_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
gemcitabine_ccle          <- merge(gemcitabine_ccle, gemcitabine_ccle_new, by = 'Cell.line.name')
gemcitabine_ccle          <- gemcitabine_ccle[, -3]
gemcitabine_ccle          <- unique(gemcitabine_ccle)

methotrexate_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 17)
methotrexate_ccle          <- methotrexate_ccle[, c(4:6)]
methotrexate_ccle_new      <- methotrexate_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
methotrexate_ccle          <- merge(methotrexate_ccle, methotrexate_ccle_new, by = 'Cell.line.name')
methotrexate_ccle          <- methotrexate_ccle[, -3]
methotrexate_ccle          <- unique(methotrexate_ccle)

oxalaplatin_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 19)
oxalaplatin_ccle          <- oxalaplatin_ccle[, c(4:6)]
oxalaplatin_ccle_new      <- oxalaplatin_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
oxalaplatin_ccle          <- merge(oxalaplatin_ccle, oxalaplatin_ccle_new, by = 'Cell.line.name')
oxalaplatin_ccle          <- oxalaplatin_ccle[, -3]
oxalaplatin_ccle          <- unique(oxalaplatin_ccle)

paclitaxel_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 20)
paclitaxel_ccle          <- paclitaxel_ccle[, c(4:6)]
paclitaxel_ccle_new      <- paclitaxel_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
paclitaxel_ccle          <- merge(paclitaxel_ccle, paclitaxel_ccle_new, by = 'Cell.line.name')
paclitaxel_ccle          <- paclitaxel_ccle[, -3]
paclitaxel_ccle          <- unique(paclitaxel_ccle)

procarbazine_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 22)
procarbazine_ccle          <- procarbazine_ccle[, c(4:6)]
procarbazine_ccle_new      <- procarbazine_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
procarbazine_ccle          <- merge(procarbazine_ccle, procarbazine_ccle_new, by = 'Cell.line.name')
procarbazine_ccle          <- procarbazine_ccle[, -3]
procarbazine_ccle          <- unique(procarbazine_ccle)

temozolomide_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 24)
temozolomide_ccle          <- temozolomide_ccle[, c(4:6)]
temozolomide_ccle_new      <- temozolomide_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
temozolomide_ccle          <- merge(temozolomide_ccle, temozolomide_ccle_new, by = 'Cell.line.name')
temozolomide_ccle          <- temozolomide_ccle[, -3]
temozolomide_ccle          <- unique(temozolomide_ccle)

topotecan_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 26)
topotecan_ccle          <- topotecan_ccle[, c(4:6)]
topotecan_ccle_new      <- topotecan_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
topotecan_ccle          <- merge(topotecan_ccle, topotecan_ccle_new, by = 'Cell.line.name')
topotecan_ccle          <- topotecan_ccle[, -3]
topotecan_ccle          <- unique(topotecan_ccle)

vincristine_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 27)
vincristine_ccle          <- vincristine_ccle[, c(4:6)]
vincristine_ccle_new      <- vincristine_ccle %>% group_by(Cell.line.name) %>% summarise(AUC = mean(area_under_curve))
vincristine_ccle          <- merge(vincristine_ccle, vincristine_ccle_new, by = 'Cell.line.name')
vincristine_ccle          <- vincristine_ccle[, -3]
vincristine_ccle          <- unique(vincristine_ccle)


#TCGA (RFS)

#BLCA
blca_phenos <- read.delim('TCGA_files/TCGA-BLCA.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 453 TCGA samples
blca_phenos_short <- blca_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
blca_phenos_not_lost <- blca_phenos_short[blca_phenos_short$lost_follow_up != 'YES', ]
blca_phenos_w_drug <- blca_phenos_not_lost[blca_phenos_not_lost$drug_name != '', ]
blca_phenos_done <- blca_phenos_w_drug
blca_phenos_done$OS <- ifelse(blca_phenos_done$vital_status.diagnoses == 'dead', blca_phenos_done[, 5], blca_phenos_done[, 6])
blca_phenos_done$PFS <- ifelse(is.na(blca_phenos_done[, 2]), blca_phenos_done$OS, blca_phenos_done[, 2])
table(blca_phenos_done$drug_name)

#BRCA
brca_phenos <- read.delim('TCGA_files/TCGA-BRCA.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 1283 TCGA samples
brca_phenos_short <- brca_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
brca_phenos_not_lost <- brca_phenos_short[brca_phenos_short$lost_follow_up != 'YES', ]
brca_phenos_w_drug <- brca_phenos_not_lost[brca_phenos_not_lost$drug_name != '', ]
brca_phenos_done <- brca_phenos_w_drug
brca_phenos_done$OS <- ifelse(brca_phenos_done$vital_status.diagnoses == 'dead', brca_phenos_done[, 5], brca_phenos_done[, 6])
brca_phenos_done$PFS <- ifelse(is.na(brca_phenos_done[, 2]), brca_phenos_done$OS, brca_phenos_done[, 2])
table(brca_phenos_done$drug_name)

#CESC
cesc_phenos <- read.delim('TCGA_files/TCGA-CESC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 317 TCGA samples
cesc_phenos_short <- cesc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
cesc_phenos_not_lost <- cesc_phenos_short[cesc_phenos_short$lost_follow_up != 'YES', ]
cesc_phenos_w_drug <- cesc_phenos_not_lost[cesc_phenos_not_lost$drug_name != '', ]
cesc_phenos_done <- cesc_phenos_w_drug
cesc_phenos_done$OS <- ifelse(cesc_phenos_done$vital_status.diagnoses == 'dead', cesc_phenos_done[, 5], cesc_phenos_done[, 6])
cesc_phenos_done$PFS <- ifelse(is.na(cesc_phenos_done[, 2]), cesc_phenos_done$OS, cesc_phenos_done[, 2])
table(cesc_phenos_done$drug_name)

#CHOL
chol_phenos <- read.delim('TCGA_files/TCGA-CHOL.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 71 TCGA samples
chol_phenos_short <- chol_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
chol_phenos_not_lost <- chol_phenos_short[chol_phenos_short$lost_follow_up != 'YES', ]
chol_phenos_w_drug <- chol_phenos_not_lost[chol_phenos_not_lost$drug_name != '', ]
chol_phenos_done <- chol_phenos_w_drug
chol_phenos_done$OS <- ifelse(chol_phenos_done$vital_status.diagnoses == 'dead', chol_phenos_done[, 5], chol_phenos_done[, 6])
chol_phenos_done$PFS <- ifelse(is.na(chol_phenos_done[, 2]), chol_phenos_done$OS, chol_phenos_done[, 2])
table(chol_phenos_done$drug_name)

#COAD
coad_phenos <- read.delim('TCGA_files/TCGA-COAD.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 570 TCGA samples
coad_phenos_short <- coad_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
coad_phenos_not_lost <- coad_phenos_short[coad_phenos_short$lost_follow_up != 'YES', ]
coad_phenos_w_drug <- coad_phenos_not_lost[coad_phenos_not_lost$drug_name != '', ]
coad_phenos_done <- coad_phenos_w_drug
coad_phenos_done$OS <- ifelse(coad_phenos_done$vital_status.diagnoses == 'dead', coad_phenos_done[, 5], coad_phenos_done[, 6])
coad_phenos_done$PFS <- ifelse(is.na(coad_phenos_done[, 2]), coad_phenos_done$OS, coad_phenos_done[, 2])
table(coad_phenos_done$drug_name)

#DLBC
dlbc_phenos <- read.delim('TCGA_files/TCGA-DLBC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 52 TCGA samples
dlbc_phenos_short <- dlbc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
dlbc_phenos_not_lost <- dlbc_phenos_short[dlbc_phenos_short$lost_follow_up != 'YES', ]
dlbc_phenos_w_drug <- dlbc_phenos_not_lost[dlbc_phenos_not_lost$drug_name != '', ]
dlbc_phenos_done <- dlbc_phenos_w_drug
dlbc_phenos_done$OS <- ifelse(dlbc_phenos_done$vital_status.diagnoses == 'dead', dlbc_phenos_done[, 5], dlbc_phenos_done[, 6])
dlbc_phenos_done$PFS <- ifelse(is.na(dlbc_phenos_done[, 2]), dlbc_phenos_done$OS, dlbc_phenos_done[, 2])
table(dlbc_phenos_done$drug_name)

#ESCA
esca_phenos <- read.delim('TCGA_files/TCGA-ESCA.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 251 TCGA samples
esca_phenos_short <- esca_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
esca_phenos_not_lost <- esca_phenos_short[esca_phenos_short$lost_follow_up != 'YES', ]
esca_phenos_w_drug <- esca_phenos_not_lost[esca_phenos_not_lost$drug_name != '', ]
esca_phenos_done <- esca_phenos_w_drug
esca_phenos_done$OS <- ifelse(esca_phenos_done$vital_status.diagnoses == 'dead', esca_phenos_done[, 5], esca_phenos_done[, 6])
esca_phenos_done$PFS <- ifelse(is.na(esca_phenos_done[, 2]), esca_phenos_done$OS, esca_phenos_done[, 2])
table(esca_phenos_done$drug_name)

#HNSC
hnsc_phenos <- read.delim('TCGA_files/TCGA-HNSC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 612 TCGA samples
hnsc_phenos_short <- hnsc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
hnsc_phenos_not_lost <- hnsc_phenos_short[hnsc_phenos_short$lost_follow_up != 'YES', ]
hnsc_phenos_w_drug <- hnsc_phenos_not_lost[hnsc_phenos_not_lost$drug_name != '', ]
hnsc_phenos_done <- hnsc_phenos_w_drug
hnsc_phenos_done$OS <- ifelse(hnsc_phenos_done$vital_status.diagnoses == 'dead', hnsc_phenos_done[, 5], hnsc_phenos_done[, 6])
hnsc_phenos_done$PFS <- ifelse(is.na(hnsc_phenos_done[, 2]), hnsc_phenos_done$OS, hnsc_phenos_done[, 2])
table(hnsc_phenos_done$drug_name)

#KIRC
kirc_phenos <- read.delim('TCGA_files/TCGA-KIRC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 985 TCGA samples
kirc_phenos_short <- kirc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
kirc_phenos_not_lost <- kirc_phenos_short[kirc_phenos_short$lost_follow_up != 'YES', ]
kirc_phenos_w_drug <- kirc_phenos_not_lost[kirc_phenos_not_lost$drug_name != '', ]
kirc_phenos_done <- kirc_phenos_w_drug
kirc_phenos_done$OS <- ifelse(kirc_phenos_done$vital_status.diagnoses == 'dead', kirc_phenos_done[, 5], kirc_phenos_done[, 6])
kirc_phenos_done$PFS <- ifelse(is.na(kirc_phenos_done[, 2]), kirc_phenos_done$OS, kirc_phenos_done[, 2])
table(kirc_phenos_done$drug_name)

#LIHC
lihc_phenos <- read.delim('TCGA_files/TCGA-LIHC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 469 TCGA samples
lihc_phenos_short <- lihc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
lihc_phenos_not_lost <- lihc_phenos_short[lihc_phenos_short$lost_follow_up != 'YES', ]
lihc_phenos_w_drug <- lihc_phenos_not_lost[lihc_phenos_not_lost$drug_name != '', ]
lihc_phenos_done <- lihc_phenos_w_drug
lihc_phenos_done$OS <- ifelse(lihc_phenos_done$vital_status.diagnoses == 'dead', lihc_phenos_done[, 5], lihc_phenos_done[, 6])
lihc_phenos_done$PFS <- ifelse(is.na(lihc_phenos_done[, 2]), lihc_phenos_done$OS, lihc_phenos_done[, 2])
table(lihc_phenos_done$drug_name)

#LUAD
luad_phenos <- read.delim('TCGA_files/TCGA-LUAD.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 877 TCGA samples
luad_phenos_short <- luad_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
luad_phenos_not_lost <- luad_phenos_short[luad_phenos_short$lost_follow_up != 'YES', ]
luad_phenos_w_drug <- luad_phenos_not_lost[luad_phenos_not_lost$drug_name != '', ]
luad_phenos_done <- luad_phenos_w_drug
luad_phenos_done$OS <- ifelse(luad_phenos_done$vital_status.diagnoses == 'dead', luad_phenos_done[, 5], luad_phenos_done[, 6])
luad_phenos_done$PFS <- ifelse(is.na(luad_phenos_done[, 2]), luad_phenos_done$OS, luad_phenos_done[, 2])
table(luad_phenos_done$drug_name)

#LUSC
lusc_phenos <- read.delim('TCGA_files/TCGA-LUSC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 765 TCGA samples
lusc_phenos_short <- lusc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
lusc_phenos_not_lost <- lusc_phenos_short[lusc_phenos_short$lost_follow_up != 'YES', ]
lusc_phenos_w_drug <- lusc_phenos_not_lost[lusc_phenos_not_lost$drug_name != '', ]
lusc_phenos_done <- lusc_phenos_w_drug
lusc_phenos_done$OS <- ifelse(lusc_phenos_done$vital_status.diagnoses == 'dead', lusc_phenos_done[, 5], lusc_phenos_done[, 6])
lusc_phenos_done$PFS <- ifelse(is.na(lusc_phenos_done[, 2]), lusc_phenos_done$OS, lusc_phenos_done[, 2])
table(lusc_phenos_done$drug_name)

#MESO
meso_phenos <- read.delim('TCGA_files/TCGA-MESO.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 88 TCGA samples
meso_phenos_short <- meso_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
meso_phenos_not_lost <- meso_phenos_short[meso_phenos_short$lost_follow_up != 'YES', ]
meso_phenos_w_drug <- meso_phenos_not_lost[meso_phenos_not_lost$drug_name != '', ]
meso_phenos_done <- meso_phenos_w_drug
meso_phenos_done$OS <- ifelse(meso_phenos_done$vital_status.diagnoses == 'dead', meso_phenos_done[, 5], meso_phenos_done[, 6])
meso_phenos_done$PFS <- ifelse(is.na(meso_phenos_done[, 2]), meso_phenos_done$OS, meso_phenos_done[, 2])
table(meso_phenos_done$drug_name)

#OV
ov_phenos <- read.delim('TCGA_files/TCGA-OV.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 758 TCGA samples
ov_phenos_short <- ov_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
ov_phenos_not_lost <- ov_phenos_short[ov_phenos_short$lost_follow_up != 'YES', ]
ov_phenos_w_drug <- ov_phenos_not_lost[ov_phenos_not_lost$drug_name != '', ]
ov_phenos_done <- ov_phenos_w_drug
ov_phenos_done$OS <- ifelse(ov_phenos_done$vital_status.diagnoses == 'dead', ov_phenos_done[, 5], ov_phenos_done[, 6])
ov_phenos_done$PFS <- ifelse(is.na(ov_phenos_done[, 2]), ov_phenos_done$OS, ov_phenos_done[, 2])
table(ov_phenos_done$drug_name)

#PAAD
paad_phenos <- read.delim('TCGA_files/TCGA-PAAD.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 223 TCGA samples
paad_phenos_short <- paad_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
paad_phenos_not_lost <- paad_phenos_short[paad_phenos_short$lost_follow_up != 'YES', ]
paad_phenos_w_drug <- paad_phenos_not_lost[paad_phenos_not_lost$drug_name != '', ]
paad_phenos_done <- paad_phenos_w_drug
paad_phenos_done$OS <- ifelse(paad_phenos_done$vital_status.diagnoses == 'dead', paad_phenos_done[, 5], paad_phenos_done[, 6])
paad_phenos_done$PFS <- ifelse(is.na(paad_phenos_done[, 2]), paad_phenos_done$OS, paad_phenos_done[, 2])
table(paad_phenos_done$drug_name)

#READ
read_phenos <- read.delim('TCGA_files/TCGA-READ.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 192 TCGA samples
read_phenos_short <- read_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
read_phenos_not_lost <- read_phenos_short[read_phenos_short$lost_follow_up != 'YES', ]
read_phenos_w_drug <- read_phenos_not_lost[read_phenos_not_lost$drug_name != '', ]
read_phenos_done <- read_phenos_w_drug
read_phenos_done$OS <- ifelse(read_phenos_done$vital_status.diagnoses == 'dead', read_phenos_done[, 5], read_phenos_done[, 6])
read_phenos_done$PFS <- ifelse(is.na(read_phenos_done[, 2]), read_phenos_done$OS, read_phenos_done[, 2])
table(read_phenos_done$drug_name)

#SARC
sarc_phenos <- read.delim('TCGA_files/TCGA-SARC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 290 TCGA samples
sarc_phenos_short <- sarc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
sarc_phenos_not_lost <- sarc_phenos_short[sarc_phenos_short$lost_follow_up != 'YES', ]
sarc_phenos_w_drug <- sarc_phenos_not_lost[sarc_phenos_not_lost$drug_name != '', ]
sarc_phenos_done <- sarc_phenos_w_drug
sarc_phenos_done$OS <- ifelse(sarc_phenos_done$vital_status.diagnoses == 'dead', sarc_phenos_done[, 5], sarc_phenos_done[, 6])
sarc_phenos_done$PFS <- ifelse(is.na(sarc_phenos_done[, 2]), sarc_phenos_done$OS, sarc_phenos_done[, 2])
table(sarc_phenos_done$drug_name)

#SKCM
skcm_phenos <- read.delim('TCGA_files/TCGA-SKCM.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 477 TCGA samples
skcm_phenos_short <- skcm_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
skcm_phenos_not_lost <- skcm_phenos_short[skcm_phenos_short$lost_follow_up != 'YES', ]
skcm_phenos_w_drug <- skcm_phenos_not_lost[skcm_phenos_not_lost$drug_name != '', ]
skcm_phenos_done <- skcm_phenos_w_drug
skcm_phenos_done$OS <- ifelse(skcm_phenos_done$vital_status.diagnoses == 'dead', skcm_phenos_done[, 5], skcm_phenos_done[, 6])
skcm_phenos_done$PFS <- ifelse(is.na(skcm_phenos_done[, 2]), skcm_phenos_done$OS, skcm_phenos_done[, 2])
table(skcm_phenos_done$drug_name)

#STAD
stad_phenos <- read.delim('TCGA_files/TCGA-STAD.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 544 TCGA samples
stad_phenos_short <- stad_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
stad_phenos_not_lost <- stad_phenos_short[stad_phenos_short$lost_follow_up != 'YES', ]
stad_phenos_w_drug <- stad_phenos_not_lost[stad_phenos_not_lost$drug_name != '', ]
stad_phenos_done <- stad_phenos_w_drug
stad_phenos_done$OS <- ifelse(stad_phenos_done$vital_status.diagnoses == 'dead', stad_phenos_done[, 5], stad_phenos_done[, 6])
stad_phenos_done$PFS <- ifelse(is.na(stad_phenos_done[, 2]), stad_phenos_done$OS, stad_phenos_done[, 2])
table(stad_phenos_done$drug_name)

#TGCT
tgct_phenos <- read.delim('TCGA_files/TCGA-TGCT.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 156 TCGA samples
tgct_phenos_short <- tgct_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
tgct_phenos_not_lost <- tgct_phenos_short[tgct_phenos_short$lost_follow_up != 'YES', ]
tgct_phenos_w_drug <- tgct_phenos_not_lost[tgct_phenos_not_lost$drug_name != '', ]
tgct_phenos_done <- tgct_phenos_w_drug
tgct_phenos_done$OS <- ifelse(tgct_phenos_done$vital_status.diagnoses == 'dead', tgct_phenos_done[, 5], tgct_phenos_done[, 6])
tgct_phenos_done$PFS <- ifelse(is.na(tgct_phenos_done[, 2]), tgct_phenos_done$OS, tgct_phenos_done[, 2])
table(tgct_phenos_done$drug_name)

#UCEC
ucec_phenos <- read.delim('TCGA_files/TCGA-UCEC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 605 TCGA samples
ucec_phenos_short <- ucec_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
ucec_phenos_not_lost <- ucec_phenos_short[ucec_phenos_short$lost_follow_up != 'YES', ]
ucec_phenos_w_drug <- ucec_phenos_not_lost[ucec_phenos_not_lost$drug_name != '', ]
ucec_phenos_done <- ucec_phenos_w_drug
ucec_phenos_done$OS <- ifelse(ucec_phenos_done$vital_status.diagnoses == 'dead', ucec_phenos_done[, 5], ucec_phenos_done[, 6])
ucec_phenos_done$PFS <- ifelse(is.na(ucec_phenos_done[, 2]), ucec_phenos_done$OS, ucec_phenos_done[, 2])
table(ucec_phenos_done$drug_name)

#UCS
ucs_phenos <- read.delim('TCGA_files/TCGA-UCS.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 63 TCGA samples
ucs_phenos_short <- ucs_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
ucs_phenos_not_lost <- ucs_phenos_short[ucs_phenos_short$lost_follow_up != 'YES', ]
ucs_phenos_w_drug <- ucs_phenos_not_lost[ucs_phenos_not_lost$drug_name != '', ]
ucs_phenos_done <- ucs_phenos_w_drug
ucs_phenos_done$OS <- ifelse(ucs_phenos_done$vital_status.diagnoses == 'dead', ucs_phenos_done[, 5], ucs_phenos_done[, 6])
ucs_phenos_done$PFS <- ifelse(is.na(ucs_phenos_done[, 2]), ucs_phenos_done$OS, ucs_phenos_done[, 2])
table(ucs_phenos_done$drug_name)

### create dummy binary variables for most/least sensitive -------

# GDSC

cisplatin$most_sensitive        <- ifelse(cisplatin$LN_IC50 < quantile(cisplatin$LN_IC50, probs = .20), 1, 0)
cisplatin$least_sensitive       <- ifelse(cisplatin$LN_IC50 > quantile(cisplatin$LN_IC50, probs = .80), 1, 0)

etoposide$most_sensitive        <- ifelse(etoposide$LN_IC50 < quantile(etoposide$LN_IC50, probs = .20), 1, 0)
etoposide$least_sensitive       <- ifelse(etoposide$LN_IC50 > quantile(etoposide$LN_IC50, probs = .80), 1, 0)

gemcitabine$most_sensitive      <- ifelse(gemcitabine$LN_IC50 < quantile(gemcitabine$LN_IC50, probs = .20), 1, 0)
gemcitabine$least_sensitive     <- ifelse(gemcitabine$LN_IC50 > quantile(gemcitabine$LN_IC50, probs = .80), 1, 0)

methotrexate$most_sensitive     <- ifelse(methotrexate$LN_IC50 < quantile(methotrexate$LN_IC50, probs = .20), 1, 0)
methotrexate$least_sensitive    <- ifelse(methotrexate$LN_IC50 > quantile(methotrexate$LN_IC50, probs = .80), 1, 0)

## CCLE
carboplatin_ccle$most_sensitive           <- ifelse(carboplatin_ccle$AUC < quantile(carboplatin_ccle$AUC, probs = .20), 1, 0)
carboplatin_ccle$least_sensitive          <- ifelse(carboplatin_ccle$AUC > quantile(carboplatin_ccle$AUC, probs = .80), 1, 0)

cyclophosphamide_ccle$most_sensitive          <- ifelse(cyclophosphamide_ccle$AUC < quantile(cyclophosphamide_ccle$AUC, probs = .20), 1, 0)
cyclophosphamide_ccle$least_sensitive         <- ifelse(cyclophosphamide_ccle$AUC > quantile(cyclophosphamide_ccle$AUC, probs = .80), 1, 0)

dacarbazine_ccle$most_sensitive           <- ifelse(dacarbazine_ccle$AUC < quantile(dacarbazine_ccle$AUC, probs = .20), 1, 0)
dacarbazine_ccle$least_sensitive          <- ifelse(dacarbazine_ccle$AUC > quantile(dacarbazine_ccle$AUC, probs = .80), 1, 0)

docetaxel_ccle$most_sensitive           <- ifelse(docetaxel_ccle$AUC < quantile(docetaxel_ccle$AUC, probs = .20), 1, 0)
docetaxel_ccle$least_sensitive          <- ifelse(docetaxel_ccle$AUC > quantile(docetaxel_ccle$AUC, probs = .80), 1, 0)

doxorubicin_ccle$most_sensitive                 <- ifelse(doxorubicin_ccle$AUC < quantile(doxorubicin_ccle$AUC, probs = .20), 1, 0)
doxorubicin_ccle$least_sensitive                <- ifelse(doxorubicin_ccle$AUC > quantile(doxorubicin_ccle$AUC, probs = .80), 1, 0)

etoposide_ccle$most_sensitive           <- ifelse(etoposide_ccle$AUC < quantile(etoposide_ccle$AUC, probs = .20), 1, 0)
etoposide_ccle$least_sensitive          <- ifelse(etoposide_ccle$AUC > quantile(etoposide_ccle$AUC, probs = .80), 1, 0)

fluorouracil_ccle$most_sensitive           <- ifelse(fluorouracil_ccle$AUC < quantile(fluorouracil_ccle$AUC, probs = .20), 1, 0)
fluorouracil_ccle$least_sensitive          <- ifelse(fluorouracil_ccle$AUC > quantile(fluorouracil_ccle$AUC, probs = .80), 1, 0)

gemcitabine_ccle$most_sensitive         <- ifelse(gemcitabine_ccle$AUC < quantile(gemcitabine_ccle$AUC, probs = .20), 1, 0)
gemcitabine_ccle$least_sensitive        <- ifelse(gemcitabine_ccle$AUC > quantile(gemcitabine_ccle$AUC, probs = .80), 1, 0)

methotrexate_ccle$most_sensitive        <- ifelse(methotrexate_ccle$AUC < quantile(methotrexate_ccle$AUC, probs = .20), 1, 0)
methotrexate_ccle$least_sensitive       <- ifelse(methotrexate_ccle$AUC > quantile(methotrexate_ccle$AUC, probs = .80), 1, 0)

oxalaplatin_ccle$most_sensitive           <- ifelse(oxalaplatin_ccle$AUC < quantile(oxalaplatin_ccle$AUC, probs = .20), 1, 0)
oxalaplatin_ccle$least_sensitive          <- ifelse(oxalaplatin_ccle$AUC > quantile(oxalaplatin_ccle$AUC, probs = .80), 1, 0)

paclitaxel_ccle$most_sensitive                <- ifelse(paclitaxel_ccle$AUC < quantile(paclitaxel_ccle$AUC, probs = .20), 1, 0)
paclitaxel_ccle$least_sensitive               <- ifelse(paclitaxel_ccle$AUC > quantile(paclitaxel_ccle$AUC, probs = .80), 1, 0)

procarbazine_ccle$most_sensitive           <- ifelse(procarbazine_ccle$AUC < quantile(procarbazine_ccle$AUC, probs = .20), 1, 0)
procarbazine_ccle$least_sensitive          <- ifelse(procarbazine_ccle$AUC > quantile(procarbazine_ccle$AUC, probs = .80), 1, 0)

temozolomide_ccle$most_sensitive        <- ifelse(temozolomide_ccle$AUC < quantile(temozolomide_ccle$AUC, probs = .20), 1, 0)
temozolomide_ccle$least_sensitive       <- ifelse(temozolomide_ccle$AUC > quantile(temozolomide_ccle$AUC, probs = .80), 1, 0)

topotecan_ccle$most_sensitive           <- ifelse(topotecan_ccle$AUC < quantile(topotecan_ccle$AUC, probs = .20), 1, 0)
topotecan_ccle$least_sensitive          <- ifelse(topotecan_ccle$AUC > quantile(topotecan_ccle$AUC, probs = .80), 1, 0)

vincristine_ccle$most_sensitive           <- ifelse(vincristine_ccle$AUC < quantile(vincristine_ccle$AUC, probs = .20), 1, 0)
vincristine_ccle$least_sensitive          <- ifelse(vincristine_ccle$AUC > quantile(vincristine_ccle$AUC, probs = .80), 1, 0)

## TCGA
blca_phenos_done$most_sensitive        <- ifelse(blca_phenos_done$PFS < quantile(blca_phenos_done$PFS, probs = .20), 1, 0)
blca_phenos_done$least_sensitive       <- ifelse(blca_phenos_done$PFS > quantile(blca_phenos_done$PFS, probs = .80), 1, 0)

brca_phenos_done$most_sensitive        <- ifelse(brca_phenos_done$PFS < quantile(brca_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
brca_phenos_done$least_sensitive       <- ifelse(brca_phenos_done$PFS > quantile(brca_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

cesc_phenos_done$most_sensitive        <- ifelse(cesc_phenos_done$PFS < quantile(cesc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
cesc_phenos_done$least_sensitive       <- ifelse(cesc_phenos_done$PFS > quantile(cesc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

chol_phenos_done$most_sensitive        <- ifelse(chol_phenos_done$PFS < quantile(chol_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
chol_phenos_done$least_sensitive       <- ifelse(chol_phenos_done$PFS > quantile(chol_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

coad_phenos_done$most_sensitive        <- ifelse(coad_phenos_done$PFS < quantile(coad_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
coad_phenos_done$least_sensitive       <- ifelse(coad_phenos_done$PFS > quantile(coad_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

dlbc_phenos_done$most_sensitive        <- ifelse(dlbc_phenos_done$PFS < quantile(dlbc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
dlbc_phenos_done$least_sensitive       <- ifelse(dlbc_phenos_done$PFS > quantile(dlbc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

esca_phenos_done$most_sensitive        <- ifelse(esca_phenos_done$PFS < quantile(esca_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
esca_phenos_done$least_sensitive       <- ifelse(esca_phenos_done$PFS > quantile(esca_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

hnsc_phenos_done$most_sensitive        <- ifelse(hnsc_phenos_done$PFS < quantile(hnsc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
hnsc_phenos_done$least_sensitive       <- ifelse(hnsc_phenos_done$PFS > quantile(hnsc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

kirc_phenos_done$most_sensitive        <- ifelse(kirc_phenos_done$PFS < quantile(kirc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
kirc_phenos_done$least_sensitive       <- ifelse(kirc_phenos_done$PFS > quantile(kirc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

lihc_phenos_done$most_sensitive        <- ifelse(lihc_phenos_done$PFS < quantile(lihc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
lihc_phenos_done$least_sensitive       <- ifelse(lihc_phenos_done$PFS > quantile(lihc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

luad_phenos_done$most_sensitive        <- ifelse(luad_phenos_done$PFS < quantile(luad_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
luad_phenos_done$least_sensitive       <- ifelse(luad_phenos_done$PFS > quantile(luad_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

lusc_phenos_done$most_sensitive        <- ifelse(lusc_phenos_done$PFS < quantile(lusc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
lusc_phenos_done$least_sensitive       <- ifelse(lusc_phenos_done$PFS > quantile(lusc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

meso_phenos_done$most_sensitive        <- ifelse(meso_phenos_done$PFS < quantile(meso_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
meso_phenos_done$least_sensitive       <- ifelse(meso_phenos_done$PFS > quantile(meso_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

ov_phenos_done$most_sensitive        <- ifelse(ov_phenos_done$PFS < quantile(ov_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
ov_phenos_done$least_sensitive       <- ifelse(ov_phenos_done$PFS > quantile(ov_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

paad_phenos_done$most_sensitive        <- ifelse(paad_phenos_done$PFS < quantile(paad_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
paad_phenos_done$least_sensitive       <- ifelse(paad_phenos_done$PFS > quantile(paad_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

read_phenos_done$most_sensitive        <- ifelse(read_phenos_done$PFS < quantile(read_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
read_phenos_done$least_sensitive       <- ifelse(read_phenos_done$PFS > quantile(read_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

sarc_phenos_done$most_sensitive        <- ifelse(sarc_phenos_done$PFS < quantile(sarc_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
sarc_phenos_done$least_sensitive       <- ifelse(sarc_phenos_done$PFS > quantile(sarc_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

skcm_phenos_done$most_sensitive        <- ifelse(skcm_phenos_done$PFS < quantile(skcm_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
skcm_phenos_done$least_sensitive       <- ifelse(skcm_phenos_done$PFS > quantile(skcm_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

stad_phenos_done$most_sensitive        <- ifelse(stad_phenos_done$PFS < quantile(stad_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
stad_phenos_done$least_sensitive       <- ifelse(stad_phenos_done$PFS > quantile(stad_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

tgct_phenos_done$most_sensitive        <- ifelse(tgct_phenos_done$PFS < quantile(tgct_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
tgct_phenos_done$least_sensitive       <- ifelse(tgct_phenos_done$PFS > quantile(tgct_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

ucec_phenos_done$most_sensitive        <- ifelse(ucec_phenos_done$PFS < quantile(ucec_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
ucec_phenos_done$least_sensitive       <- ifelse(ucec_phenos_done$PFS > quantile(ucec_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

ucs_phenos_done$most_sensitive        <- ifelse(ucs_phenos_done$PFS < quantile(ucs_phenos_done$PFS, probs = .20, na.rm = TRUE), 1, 0)
ucs_phenos_done$least_sensitive       <- ifelse(ucs_phenos_done$PFS > quantile(ucs_phenos_done$PFS, probs = .80, na.rm = TRUE), 1, 0)

## write to files for later ----
# GDSC
write.csv(cisplatin, file = 'Processed_Clinical_Data/cisplatin_gdsc_clinical_processed.csv')
write.csv(etoposide, file = 'Processed_Clinical_Data/etoposide_gdsc_clinical_processed.csv')
write.csv(gemcitabine, file = 'Processed_Clinical_Data/gemcitabine_gdsc_clinical_processed.csv')
write.csv(methotrexate, file = 'Processed_Clinical_Data/methotrexate_gdsc_clinical_processed.csv')

# CCLE
write.csv(carboplatin_ccle, file = 'Processed_Clinical_Data/carboplatin_ccle_clinical_processed.csv')
write.csv(cyclophosphamide_ccle, file = 'Processed_Clinical_Data/cyclophosphamide_ccle_clinical_processed.csv')
write.csv(dacarbazine_ccle, file = 'Processed_Clinical_Data/dacarbazine_ccle_clinical_processed.csv')
write.csv(docetaxel_ccle, file = 'Processed_Clinical_Data/docetaxel_ccle_clinical_processed.csv')
write.csv(doxorubicin_ccle, file = 'Processed_Clinical_Data/doxorubicin_ccle_clinical_processed.csv')
write.csv(etoposide_ccle, file = 'Processed_Clinical_Data/etoposide_ccle_clinical_processed.csv')
write.csv(fluorouracil_ccle, file = 'Processed_Clinical_Data/fluorouracil_ccle_clinical_processed.csv')
write.csv(gemcitabine_ccle, file = 'Processed_Clinical_Data/gemcitabine_ccle_clinical_processed.csv')
write.csv(methotrexate_ccle, file = 'Processed_Clinical_Data/methotrexate_ccle_clinical_processed.csv')
write.csv(oxalaplatin_ccle, file = 'Processed_Clinical_Data/oxalaplatin_ccle_clinical_processed.csv')
write.csv(paclitaxel_ccle, file = 'Processed_Clinical_Data/paclitaxel_ccle_clinical_processed.csv')
write.csv(procarbazine_ccle, file = 'Processed_Clinical_Data/procarbazine_ccle_clinical_processed.csv')
write.csv(temozolomide_ccle, file = 'Processed_Clinical_Data/temozolomide_ccle_clinical_processed.csv')
write.csv(topotecan_ccle, file = 'Processed_Clinical_Data/topotecan_ccle_clinical_processed.csv')
write.csv(vincristine_ccle, file = 'Processed_Clinical_Data/vincristine_ccle_clinical_processed.csv')

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
