###this is the script for running the entire project
###this is completely self-contained along with necessary data files MAKE SURE  THIS IS TRUE!!!!

##load (and install if needed) required packages -------
if (!require ('openxlsx')) install.packages('openxlsx')
library(openxlsx) # to read in excel files

if (!require ('biomaRt')) BiocInstaller::biocLite('biomaRt')
library(biomaRt) # for conversion between gene names, Entrez IDs, and ENSG IDs

if (!require ('ggplot2')) install.packages('ggplot2')
library(ggplot2) # general plotting functions

if (!require ('glmnet')) install.packages('glmnet')
library(glmnet) # for model building

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

### load all gene expression data --------
### CELL LINES
# GDSC
gdsc_rna_seq <- read.csv('Gene_Expression_Files/gdsc_rnaseq.csv', header = TRUE, stringsAsFactors = FALSE)
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$X, unique = TRUE) # ENSG IDs as row names
gdsc_rna_seq <- gdsc_rna_seq[, -1] # remove gene names
gdsc_rna_seq <- gdsc_rna_seq[-c(17738,17739), ] # removes NAs that were there for some reason
# dimensions now: 17737 genes (by ENSG ID) x 962 cell lines (by COSMIC ID)

# CCLE
ccle_microarray <- read.xlsx('Gene_Expression_Files/Gene_expression_data_no blood_good.xlsx', 1)
rownames(ccle_microarray) <- make.names(ccle_microarray$Description, unique = TRUE)
ccle_microarray <- ccle_microarray[, -1]
# dimensions now: 18988 (by gene name) x 855 cell lines (by name and type)

# NCI 60
nci60_microarray <- read.xlsx('Gene_Expression_Files/GE and DS together no blood.xlsx', 1)
rownames(nci60_microarray) <- nci60_microarray$Entrez.gene.id.e
nci60_microarray <- nci60_microarray[, -1]
# dimensions now: 25675 genes (by entrez ID) x 52 cell lines (name and type)

### TUMOR SETS (TCGA)
# ACC
acc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-ACC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(acc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(acc_tcga_rnaseq))
# dimesnions now: 60483 genes (by ENSG ID) x 79 cases (TCGA ID)

# BLCA
blca_tcga_rnaseq <- read.delim('TCGA_files/TCGA-BLCA.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(blca_tcga_rnaseq) <- gsub('\\..*$', '', rownames(blca_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 430 cases (TCGA ID)

#BRCA
brca_tcga_rnaseq <- read.delim('TCGA_files/TCGA-BRCA.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(brca_tcga_rnaseq) <- gsub('\\..*$', '', rownames(brca_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 1217 cases (TCGA ID)

#CESC
cesc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-CESC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(cesc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(cesc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 309 cases (TCGA ID)

#CHOL
chol_tcga_rnaseq <- read.delim('TCGA_files/TCGA-CHOL.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(chol_tcga_rnaseq) <- gsub('\\..*$', '', rownames(chol_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 45 cases (TCGA ID)

#COAD
coad_tcga_rnaseq <- read.delim('TCGA_files/TCGA-COAD.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(coad_tcga_rnaseq) <- gsub('\\..*$', '', rownames(coad_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 512 cases (TCGA ID)

#DLBC
dlbc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-DLBC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(dlbc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(dlbc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 48 cases (TCGA ID)

#ESCA
esca_tcga_rnaseq <- read.delim('TCGA_files/TCGA-ESCA.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(esca_tcga_rnaseq) <- gsub('\\..*$', '', rownames(esca_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 173 cases (TCGA ID)

#GBM
gbm_tcga_rnaseq <- read.delim('TCGA_files/TCGA-GBM.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(gbm_tcga_rnaseq) <- gsub('\\..*$', '', rownames(gbm_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 173 cases (TCGA ID)

#HNSC
hnsc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-HNSC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(hnsc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(hnsc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 546 cases (TCGA ID)

#KICH
kich_tcga_rnaseq <- read.delim('TCGA_files/TCGA-KICH.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(kich_tcga_rnaseq) <- gsub('\\..*$', '', rownames(kich_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 89 cases (TCGA ID)

#KIRC
kirc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-KIRC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(kirc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(kirc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 607 cases (TCGA ID)

#KIRP
kirp_tcga_rnaseq <- read.delim('TCGA_files/TCGA-KIRP.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(kirp_tcga_rnaseq) <- gsub('\\..*$', '', rownames(kirp_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 321 cases (TCGA ID)

#LAML
laml_tcga_rnaseq <- read.delim('TCGA_files/TCGA-LAML.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(laml_tcga_rnaseq) <- gsub('\\..*$', '', rownames(laml_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 151 cases (TCGA ID)

#LIHC
lihc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-LIHC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(lihc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(lihc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 424 cases (TCGA ID)

#LUAD
luad_tcga_rnaseq <- read.delim('TCGA_files/TCGA-LUAD.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(luad_tcga_rnaseq) <- gsub('\\..*$', '', rownames(luad_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 585 cases (TCGA ID)

#LUSC
lusc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-LUSC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(lusc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(lusc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 550 cases (TCGA ID)

#MESO
meso_tcga_rnaseq <- read.delim('TCGA_files/TCGA-MESO.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(meso_tcga_rnaseq) <- gsub('\\..*$', '', rownames(meso_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 86 cases (TCGA ID)

#OV
ov_tcga_rnaseq <- read.delim('TCGA_files/TCGA-OV.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(ov_tcga_rnaseq) <- gsub('\\..*$', '', rownames(ov_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 379 cases (TCGA ID)

#PAAD
paad_tcga_rnaseq <- read.delim('TCGA_files/TCGA-PAAD.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(paad_tcga_rnaseq) <- gsub('\\..*$', '', rownames(paad_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 182 cases (TCGA ID)

#PCPG
pcpg_tcga_rnaseq <- read.delim('TCGA_files/TCGA-PCPG.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(pcpg_tcga_rnaseq) <- gsub('\\..*$', '', rownames(pcpg_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 186 cases (TCGA ID)

#PRAD
prad_tcga_rnaseq <- read.delim('TCGA_files/TCGA-PRAD.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(prad_tcga_rnaseq) <- gsub('\\..*$', '', rownames(prad_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 551 cases (TCGA ID)

#READ
read_tcga_rnaseq <- read.delim('TCGA_files/TCGA-READ.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(read_tcga_rnaseq) <- gsub('\\..*$', '', rownames(read_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 178 cases (TCGA ID)

#SARC
sarc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-SARC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(sarc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(sarc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 265 cases (TCGA ID)

#SKCM
skcm_tcga_rnaseq <- read.delim('TCGA_files/TCGA-SKCM.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(skcm_tcga_rnaseq) <- gsub('\\..*$', '', rownames(skcm_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 472 cases (TCGA ID)

#STAD
stad_tcga_rnaseq <- read.delim('TCGA_files/TCGA-STAD.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(stad_tcga_rnaseq) <- gsub('\\..*$', '', rownames(stad_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 407 cases (TCGA ID)

#TGCT
tgct_tcga_rnaseq <- read.delim('TCGA_files/TCGA-TGCT.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(tgct_tcga_rnaseq) <- gsub('\\..*$', '', rownames(tgct_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 156 cases (TCGA ID)

#THCA
thca_tcga_rnaseq <- read.delim('TCGA_files/TCGA-THCA.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(thca_tcga_rnaseq) <- gsub('\\..*$', '', rownames(thca_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 568 cases (TCGA ID)

#THYM
thym_tcga_rnaseq <- read.delim('TCGA_files/TCGA-THYM.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(thym_tcga_rnaseq) <- gsub('\\..*$', '', rownames(thym_tcga_rnaseq))
# dimensions now: 6043 genes (by ENSG ID) x 121 cases (TCGA ID)

#UCEC
ucec_tcga_rnaseq <- read.delim('TCGA_files/TCGA-UCEC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(ucec_tcga_rnaseq) <- gsub('\\..*$', '', rownames(ucec_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 583 cases (TCGA ID)

#UCS
ucs_tcga_rnaseq <- read.delim('TCGA_files/TCGA-UCS.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(ucs_tcga_rnaseq) <- gsub('\\..*$', '', rownames(ucs_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 56 cases (TCGA ID)

#UVM
uvm_tcga_rnaseq <- read.delim('TCGA_files/TCGA-UVM.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(uvm_tcga_rnaseq) <- gsub('\\..*$', '', rownames(uvm_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 80 cases (TCGA ID)

## make sure all gene expression sets have the same features ------
# makes sure model doesn't include something not in another set
# gdsc is by ENSG ID
# ccle is by gene name
# nci60 is by entrez ID
# TCGA is by ENSG ID
# convert ccle and nci60 to ENSG

# get gene names and entrez IDs
ccle_genes <- rownames(ccle_microarray)
nci60_genes <- rownames(nci60_microarray)

# convert them to ENSG IDs
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
ccle_annots <- getBM(mart = mart, attributes = c('hgnc_symbol', 'ensembl_gene_id'), filters = 'hgnc_symbol', values = ccle_genes, uniqueRows = TRUE)
# this can be imported as ccle_annots.csv
ccle_annots <- ccle_annots[!duplicated(ccle_annots[, 2]), ]
ccle_microarray <- ccle_microarray[which(rownames(ccle_microarray) %in% ccle_annots[, 1]), ]
ccle_annots <- ccle_annots[which(ccle_annots[, 1] %in% rownames(ccle_microarray)), ]
ccle_microarray <- ccle_microarray[match(ccle_annots[, 1], rownames(ccle_microarray)), ]
rownames(ccle_microarray) <- make.names(ccle_annots[, 2], unique = TRUE)
# dimensions now: 17846 genes (by ENSG ID) x 855 cell lines (name and type)

nci60_annots <- getBM(mart = mart, attributes = c('entrezgene', 'ensembl_gene_id'), filters = 'entrezgene', values = nci60_genes, uniqueRows = TRUE)
# this ca be imported as nci60_annots.csv
nci60_annots <- nci60_annots[!duplicated(nci60_annots[, 2]), ]
nci60_microarray <- nci60_microarray[which(rownames(nci60_microarray) %in% nci60_annots[, 1]), ]
nci60_annots <- nci60_annots[which(nci60_annots[, 1] %in% rownames(nci60_microarray)), ]
nci60_microarray <- nci60_microarray[match(nci60_annots[, 1], rownames(nci60_microarray)), ]
rownames(nci60_microarray) <- make.names(nci60_annots[, 2], unique = TRUE)
# dimensions now: 22788 genes (by ENSG ID) x 52 cell lines (name and type)

# limit each gene expression set to common ENSG IDs
ensg_i_need <- Reduce(intersect, list(rownames(gdsc_rna_seq), rownames(ccle_microarray), rownames(nci60_microarray)))
#14242 ENSG ID in common

gdsc_rna_seq <- gdsc_rna_seq[which(rownames(gdsc_rna_seq) %in% ensg_i_need), ]
#14242 x 962

ccle_microarray <- ccle_microarray[which(rownames(ccle_microarray) %in% ensg_i_need), ]
#14242 x 855

nci60_microarray <- nci60_microarray[which(rownames(nci60_microarray) %in% ensg_i_need), ]
#14242 x 52

acc_tcga_rnaseq <- acc_tcga_rnaseq[which(rownames(acc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 79

blca_tcga_rnaseq <- blca_tcga_rnaseq[which(rownames(blca_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 430

brca_tcga_rnaseq <- brca_tcga_rnaseq[which(rownames(brca_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 1217

cesc_tcga_rnaseq <- cesc_tcga_rnaseq[which(rownames(cesc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 309

chol_tcga_rnaseq <- chol_tcga_rnaseq[which(rownames(chol_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 45

coad_tcga_rnaseq <- coad_tcga_rnaseq[which(rownames(coad_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 512

dlbc_tcga_rnaseq <- dlbc_tcga_rnaseq[which(rownames(dlbc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 48

esca_tcga_rnaseq <- esca_tcga_rnaseq[which(rownames(esca_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 173

gbm_tcga_rnaseq <- gbm_tcga_rnaseq[which(rownames(gbm_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 173

hnsc_tcga_rnaseq <- hnsc_tcga_rnaseq[which(rownames(hnsc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 546

kich_tcga_rnaseq <- kich_tcga_rnaseq[which(rownames(kich_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 89

kirc_tcga_rnaseq <- kirc_tcga_rnaseq[which(rownames(kirc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 607

kirp_tcga_rnaseq <- kirp_tcga_rnaseq[which(rownames(kirp_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 321

laml_tcga_rnaseq <- laml_tcga_rnaseq[which(rownames(laml_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 151

lihc_tcga_rnaseq <- lihc_tcga_rnaseq[which(rownames(lihc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 424

luad_tcga_rnaseq <- luad_tcga_rnaseq[which(rownames(luad_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 585

lusc_tcga_rnaseq <- lusc_tcga_rnaseq[which(rownames(lusc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 550

meso_tcga_rnaseq <- meso_tcga_rnaseq[which(rownames(meso_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 86

ov_tcga_rnaseq <- ov_tcga_rnaseq[which(rownames(ov_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 379

paad_tcga_rnaseq <- paad_tcga_rnaseq[which(rownames(paad_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 182

pcpg_tcga_rnaseq <- pcpg_tcga_rnaseq[which(rownames(pcpg_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 186

prad_tcga_rnaseq <- prad_tcga_rnaseq[which(rownames(prad_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 551

read_tcga_rnaseq <- read_tcga_rnaseq[which(rownames(read_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 177

sarc_tcga_rnaseq <- sarc_tcga_rnaseq[which(rownames(sarc_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 265

skcm_tcga_rnaseq <- skcm_tcga_rnaseq[which(rownames(skcm_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 472

stad_tcga_rnaseq <- stad_tcga_rnaseq[which(rownames(stad_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 407

tgct_tcga_rnaseq <- tgct_tcga_rnaseq[which(rownames(tgct_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 156

thca_tcga_rnaseq <- thca_tcga_rnaseq[which(rownames(thca_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 568

thym_tcga_rnaseq <- thym_tcga_rnaseq[which(rownames(thym_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 121

ucec_tcga_rnaseq <- ucec_tcga_rnaseq[which(rownames(ucec_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 583

ucs_tcga_rnaseq <- ucs_tcga_rnaseq[which(rownames(ucs_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 56

uvm_tcga_rnaseq <- uvm_tcga_rnaseq[which(rownames(uvm_tcga_rnaseq) %in% ensg_i_need), ]
#14242 x 80

### import clinical data ---------
# GDSC (ln(IC50))
bleomycin     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 2)
# 766 x 14
bleomycin     <- bleomycin[, c(3:6, 12)]

camptothecan  <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 3)
# 678 x 14
camptothecan  <- camptothecan[, c(3:6, 12)]

cisplatin     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 4)
# 680 x 14
cisplatin     <- cisplatin[, c(3:6, 12)]

cytarabine    <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 5)
# 676 x 14
cytarabine    <- cytarabine[, c(3:6, 12)]

dox           <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 6)
# 711 x 14
dox           <- dox[, c(3:6, 12)]

etoposide     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 7)
# 718 x 14
etoposide     <- etoposide[, c(3:6, 12)]

gemcitabine   <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 8)
# 707 x 14
gemcitabine   <- gemcitabine[, c(3:6, 12)]

methotrexate  <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 9)
# 679 x 14
methotrexate  <- methotrexate[, c(3:6, 12)]

mitomycin     <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 10)
# 712 x 14
mitomycin     <- mitomycin[, c(3:6, 12)]

sn38          <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 11)
# 761 x 14
sn38          <- sn38[, c(3:6, 12)]

temozolomide  <- read.xlsx('Clinical_Files/v17_fitted_dose_response_noblood_breakdown_DNAdamageagents.xlsx', 12)
# 739 x 14
temozolomide  <- temozolomide[, c(3:6, 12)]

# cell line info for breaking up into TCGA categories
gdsc_cell_line_info <- read.csv('Clinical_Files/Cell_listTue Dec 18 19_10_36 2018 (version 1).csv', header = TRUE, stringsAsFactors = FALSE)
gdsc_cell_line_info <- gdsc_cell_line_info[, c(1,2,5,6)]
colnames(gdsc_cell_line_info)[c(3,4)] <- c('TCGA_class', 'Tissue_subtype')

#CCLE (AUC)
bleomycin_ccle    <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 2)
# 615 x 16
bleomycin_ccle    <- bleomycin_ccle[, c(5:7)]

cytarabine_ccle   <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 7)
# 653 x 15
cytarabine_ccle   <- cytarabine_ccle[, c(4:6)]

dox_ccle          <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 12)
# 660 x 15
dox_ccle          <- dox_ccle[, c(4:6)]

etoposide_ccle    <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 13)
# 649 x 15
etoposide_ccle    <- etoposide_ccle[, c(4:6)]

gemcitabine_ccle  <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 15)
# 589 x 15
gemcitabine_ccle  <- gemcitabine_ccle[, c(4:6)]

methotrexate_ccle <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 17)
# 622 x 15
methotrexate_ccle <- methotrexate_ccle[, c(4:6)]

mitomycin_ccle    <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 18)
# 647 x 15
mitomycin_ccle    <- mitomycin_ccle[, c(4:6)]

sn38_ccle         <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 23)
# 586 x 15
sn38_ccle         <- sn38_ccle[, c(4:6)]

temozolomide_ccle <- read.xlsx('Clinical_Files/Drug sensitivity data - cytotoxics only - no heme.xlsx', 24)
# 646 x 15
temozolomide_ccle <- temozolomide_ccle[, c(4:6)]


#NCI60 (z-score)
nci60_sensitivity <- read.csv('Clinical_Files/DTP_NCI60_ZSCORE.csv', header = TRUE, stringsAsFactors = FALSE)
nci60_sensitivity <- nci60_sensitivity[, -c(1,63:71)]
# 21739 x 61

bleomycin_nci60     <- nci60_sensitivity[2064, ]
bleomycin_nci60     <- bleomycin_nci60[, -1]

cisplatin_nci60     <- nci60_sensitivity[1977, ]
cisplatin_nci60     <- cisplatin_nci60[, -1]

cytarabine_nci60    <- nci60_sensitivity[3245, ]
cytarabine_nci60    <- cytarabine_nci60[, -1]

dox_nci60           <- nci60_sensitivity[20683, ]
dox_nci60           <- dox_nci60[, -1]

etoposide_nci60     <- nci60_sensitivity[2301, ]
etoposide_nci60     <- etoposide_nci60[, -1]

gemcitabine_nci60   <- nci60_sensitivity[4533, ]
gemcitabine_nci60   <- gemcitabine_nci60[, -1]

methotrexate_nci60  <- nci60_sensitivity[13, ]
methotrexate_nci60  <- methotrexate_nci60[, -1]

mitomycin_nci60     <- nci60_sensitivity[581, ]
mitomycin_nci60     <- mitomycin_nci60[, -1]


#TCGA (RFS)
#ACC
acc_phenos <- read.delim('TCGA_files/TCGA-ACC.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 97 TCGA samples
acc_phenos_short <- acc_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
acc_phenos_not_lost <- acc_phenos_short[acc_phenos_short$lost_follow_up != 'YES', ]
acc_phenos_w_drug <- acc_phenos_not_lost[acc_phenos_not_lost$drug_name != '', ]
acc_phenos_done <- acc_phenos_w_drug
acc_phenos_done$OS <- ifelse(acc_phenos_done$vital_status.diagnoses == 'dead', acc_phenos_done[, 5], acc_phenos_done[, 6])
acc_phenos_done$PFS <- ifelse(is.na(acc_phenos_done[, 2]), acc_phenos_done$OS, acc_phenos_done[, 2])
table(acc_phenos_done$drug_name)

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

#GBM
gbm_phenos <- read.delim('TCGA_files/TCGA-GBM.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 670 TCGA samples
gbm_phenos_short <- gbm_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
gbm_phenos_not_lost <- gbm_phenos_short[gbm_phenos_short$lost_follow_up != 'YES', ]
gbm_phenos_w_drug <- gbm_phenos_not_lost[gbm_phenos_not_lost$drug_name != '', ]
gbm_phenos_done <- gbm_phenos_w_drug
gbm_phenos_done$OS <- ifelse(gbm_phenos_done$vital_status.diagnoses == 'dead', gbm_phenos_done[, 5], gbm_phenos_done[, 6])
gbm_phenos_done$PFS <- ifelse(is.na(gbm_phenos_done[, 2]), gbm_phenos_done$OS, gbm_phenos_done[, 2])
table(gbm_phenos_done$drug_name)

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

#KICH
kich_phenos <- read.delim('TCGA_files/TCGA-KICH.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 184 TCGA samples
kich_phenos_short <- kich_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
kich_phenos_not_lost <- kich_phenos_short[kich_phenos_short$lost_follow_up != 'YES', ]
kich_phenos_w_drug <- kich_phenos_not_lost[kich_phenos_not_lost$drug_name != '', ]
kich_phenos_done <- kich_phenos_w_drug
kich_phenos_done$OS <- ifelse(kich_phenos_done$vital_status.diagnoses == 'dead', kich_phenos_done[, 5], kich_phenos_done[, 6])
kich_phenos_done$PFS <- ifelse(is.na(kich_phenos_done[, 2]), kich_phenos_done$OS, kich_phenos_done[, 2])
table(kich_phenos_done$drug_name)

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

#KIRP
kirp_phenos <- read.delim('TCGA_files/TCGA-KIRP.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 380 TCGA samples
kirp_phenos_short <- kirp_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
kirp_phenos_not_lost <- kirp_phenos_short[kirp_phenos_short$lost_follow_up != 'YES', ]
kirp_phenos_w_drug <- kirp_phenos_not_lost[kirp_phenos_not_lost$drug_name != '', ]
kirp_phenos_done <- kirp_phenos_w_drug
kirp_phenos_done$OS <- ifelse(kirp_phenos_done$vital_status.diagnoses == 'dead', kirp_phenos_done[, 5], kirp_phenos_done[, 6])
kirp_phenos_done$PFS <- ifelse(is.na(kirp_phenos_done[, 2]), kirp_phenos_done$OS, kirp_phenos_done[, 2])
table(kirp_phenos_done$drug_name)

#LAML
laml_phenos <- read.delim('TCGA_files/TCGA-LAML.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 697 TCGA samples
laml_phenos_short <- laml_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
laml_phenos_not_lost <- laml_phenos_short[laml_phenos_short$lost_follow_up != 'YES', ]
laml_phenos_w_drug <- laml_phenos_not_lost[laml_phenos_not_lost$drug_name != '', ]
laml_phenos_done <- laml_phenos_w_drug
laml_phenos_done$OS <- ifelse(laml_phenos_done$vital_status.diagnoses == 'dead', laml_phenos_done[, 5], laml_phenos_done[, 6])
laml_phenos_done$PFS <- ifelse(is.na(laml_phenos_done[, 2]), laml_phenos_done$OS, laml_phenos_done[, 2])
table(laml_phenos_done$drug_name)

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

#PCPG
pcpg_phenos <- read.delim('TCGA_files/TCGA-PCPG.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 189 TCGA samples
pcpg_phenos_short <- pcpg_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
pcpg_phenos_not_lost <- pcpg_phenos_short[pcpg_phenos_short$lost_follow_up != 'YES', ]
pcpg_phenos_w_drug <- pcpg_phenos_not_lost[pcpg_phenos_not_lost$drug_name != '', ]
pcpg_phenos_done <- pcpg_phenos_w_drug
pcpg_phenos_done$OS <- ifelse(pcpg_phenos_done$vital_status.diagnoses == 'dead', pcpg_phenos_done[, 5], pcpg_phenos_done[, 6])
pcpg_phenos_done$PFS <- ifelse(is.na(pcpg_phenos_done[, 2]), pcpg_phenos_done$OS, pcpg_phenos_done[, 2])
table(pcpg_phenos_done$drug_name)

#PRAD
prad_phenos <- read.delim('TCGA_files/TCGA-PRAD.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 623 TCGA samples
prad_phenos_short <- prad_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
prad_phenos_not_lost <- prad_phenos_short[prad_phenos_short$lost_follow_up != 'YES', ]
prad_phenos_w_drug <- prad_phenos_not_lost[prad_phenos_not_lost$drug_name != '', ]
prad_phenos_done <- prad_phenos_w_drug
prad_phenos_done$OS <- ifelse(prad_phenos_done$vital_status.diagnoses == 'dead', prad_phenos_done[, 5], prad_phenos_done[, 6])
prad_phenos_done$PFS <- ifelse(is.na(prad_phenos_done[, 2]), prad_phenos_done$OS, prad_phenos_done[, 2])
table(prad_phenos_done$drug_name)

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

#THCA
thca_phenos <- read.delim('TCGA_files/TCGA-THCA.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 615 TCGA samples
thca_phenos_short <- thca_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
thca_phenos_not_lost <- thca_phenos_short[thca_phenos_short$lost_follow_up != 'YES', ]
thca_phenos_w_drug <- thca_phenos_not_lost[thca_phenos_not_lost$drug_name != '', ]
thca_phenos_done <- thca_phenos_w_drug
thca_phenos_done$OS <- ifelse(thca_phenos_done$vital_status.diagnoses == 'dead', thca_phenos_done[, 5], thca_phenos_done[, 6])
thca_phenos_done$PFS <- ifelse(is.na(thca_phenos_done[, 2]), thca_phenos_done$OS, thca_phenos_done[, 2])
table(thca_phenos_done$drug_name)

#THYM
thym_phenos <- read.delim('TCGA_files/TCGA-THYM.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 139 TCGA samples
thym_phenos_short <- thym_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
thym_phenos_not_lost <- thym_phenos_short[thym_phenos_short$lost_follow_up != 'YES', ]
thym_phenos_w_drug <- thym_phenos_not_lost[thym_phenos_not_lost$drug_name != '', ]
thym_phenos_done <- thym_phenos_w_drug
thym_phenos_done$OS <- ifelse(thym_phenos_done$vital_status.diagnoses == 'dead', thym_phenos_done[, 5], thym_phenos_done[, 6])
thym_phenos_done$PFS <- ifelse(is.na(thym_phenos_done[, 2]), thym_phenos_done$OS, thym_phenos_done[, 2])
table(thym_phenos_done$drug_name)

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

#UVM
uvm_phenos <- read.delim('TCGA_files/TCGA-UVM.GDC_phenotype.tsv', sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# 80 TCGA samples
uvm_phenos_short <- uvm_phenos[, c('submitter_id.samples', 'days_to_new_tumor_event_after_initial_treatment', 'drug_name', 'lost_follow_up', 'days_to_death.diagnoses', 'days_to_last_follow_up.diagnoses', 'vital_status.diagnoses')]
uvm_phenos_not_lost <- uvm_phenos_short[uvm_phenos_short$lost_follow_up != 'YES', ]
uvm_phenos_w_drug <- uvm_phenos_not_lost[uvm_phenos_not_lost$drug_name != '', ]
uvm_phenos_done <- uvm_phenos_w_drug
uvm_phenos_done$OS <- ifelse(uvm_phenos_done$vital_status.diagnoses == 'dead', uvm_phenos_done[, 5], uvm_phenos_done[, 6])
uvm_phenos_done$PFS <- ifelse(is.na(uvm_phenos_done[, 2]), uvm_phenos_done$OS, uvm_phenos_done[, 2])
table(uvm_phenos_done$drug_name)
## INSPECT THESE DRUG TABLES AND REMOVE WHAT CAN'T BE USED AND ALTER CODE AS NEEDED

### create dummy binary variables for most/least sensitive -------

# GDSC
bleomycin$most_sensitive        <- ifelse(bleomycin$LN_IC50 < quantile(bleomycin$LN_IC50, probs = .20), 1, 0)
bleomycin$least_sensitive       <- ifelse(bleomycin$LN_IC50 > quantile(bleomycin$LN_IC50, probs = .80), 1, 0)

camptothecan$most_sensitive     <- ifelse(camptothecan$LN_IC50 < quantile(camptothecan$LN_IC50, probs = .20), 1, 0)
camptothecan$least_sensitive    <- ifelse(camptothecan$LN_IC50 > quantile(camptothecan$LN_IC50, probs = .80), 1, 0)

cisplatin$most_sensitive        <- ifelse(cisplatin$LN_IC50 < quantile(cisplatin$LN_IC50, probs = .20), 1, 0)
cisplatin$least_sensitive       <- ifelse(cisplatin$LN_IC50 > quantile(cisplatin$LN_IC50, probs = .80), 1, 0)

cytarabine$most_sensitive       <- ifelse(cytarabine$LN_IC50 < quantile(cytarabine$LN_IC50, probs = .20), 1, 0)
cytarabine$least_sensitive      <- ifelse(cytarabine$LN_IC50 > quantile(cytarabine$LN_IC50, probs = .80), 1, 0)

dox$most_sensitive              <- ifelse(dox$LN_IC50 < quantile(dox$LN_IC50, probs = .20), 1, 0)
dox$least_sensitive             <- ifelse(dox$LN_IC50 > quantile(dox$LN_IC50, probs = .80), 1, 0)

etoposide$most_sensitive        <- ifelse(etoposide$LN_IC50 < quantile(etoposide$LN_IC50, probs = .20), 1, 0)
etoposide$least_sensitive       <- ifelse(etoposide$LN_IC50 > quantile(etoposide$LN_IC50, probs = .80), 1, 0)

gemcitabine$most_sensitive      <- ifelse(gemcitabine$LN_IC50 < quantile(gemcitabine$LN_IC50, probs = .20), 1, 0)
gemcitabine$least_sensitive     <- ifelse(gemcitabine$LN_IC50 > quantile(gemcitabine$LN_IC50, probs = .80), 1, 0)

methotrexate$most_sensitive     <- ifelse(methotrexate$LN_IC50 < quantile(methotrexate$LN_IC50, probs = .20), 1, 0)
methotrexate$least_sensitive    <- ifelse(methotrexate$LN_IC50 > quantile(methotrexate$LN_IC50, probs = .80), 1, 0)

mitomycin$most_sensitive        <- ifelse(mitomycin$LN_IC50 < quantile(mitomycin$LN_IC50, probs = .20), 1, 0)
mitomycin$least_sensitive       <- ifelse(mitomycin$LN_IC50 > quantile(mitomycin$LN_IC50, probs = .80), 1, 0)

sn38$most_sensitive             <- ifelse(sn38$LN_IC50 < quantile(sn38$LN_IC50, probs = .20), 1, 0)
sn38$least_sensitive            <- ifelse(sn38$LN_IC50 > quantile(sn38$LN_IC50, probs = .80), 1, 0)

temozolomide$most_sensitive     <- ifelse(temozolomide$LN_IC50 < quantile(temozolomide$LN_IC50, probs = .20), 1, 0)
temozolomide$least_sensitive    <- ifelse(temozolomide$LN_IC50 > quantile(temozolomide$LN_IC50, probs = .80), 1, 0)

### set up data for model building ----------
# get names of GDSC cell lines treated with each drug
bleomycin_lines           <- bleomycin$COSMIC_ID # 766
camptothecan_lines        <- camptothecan$COSMIC_ID #678
cisplatin_lines           <- cisplatin$COSMIC_ID #680
cytarabine_lines          <- cytarabine$COSMIC_ID #676
dox_lines                 <- dox$COSMIC_ID #711
etoposide_lines           <- etoposide$COSMIC_ID #718
gemcitabine_lines         <- gemcitabine$COSMIC_ID #707
methotrexate_lines        <- methotrexate$COSMIC_ID # 679
mitomycin_lines           <- mitomycin$COSMIC_ID #712
sn38_lines                <- sn38$COSMIC_ID #761
temozolomide_lines        <- temozolomide$COSMIC_ID #379

# set GDSC to usable format
gdsc <- data.frame(t(gdsc_rna_seq)) # puts predictors in columns
rownames(gdsc) <- gsub('X', '', rownames(gdsc))
# dim: 962 x 14242 

# split GDSC in half randomly
set.seed(5)
# get random numbers to use for split
random_sample <- sample(x = rownames(gdsc), size = nrow(gdsc)/2)
# create function opposite of %in%
'%ni%' <- Negate('%in%')

# get training and testing sets
gdsc_train         <- gdsc[random_sample, ] #481 x 14242

gdsc_test          <- gdsc[which(rownames(gdsc) %ni% random_sample), ] #481 x 14242

# make sure zero overlap
intersect(rownames(gdsc_train), rownames(gdsc_test))

# create training/testing sets for each drug
bleomycin_rna_seq_train    <- gdsc_train[intersect(bleomycin_lines, rownames(gdsc_train)), ]
# 381 x 14242
bleomycin_rna_seq_test     <- gdsc_test[intersect(bleomycin_lines, rownames(gdsc_test)), ]
# 385 x 14242
camptothecan_rna_seq_train  <- gdsc_train[intersect(camptothecan_lines, rownames(gdsc_train)), ]
# 336 x 14242
camptothecan_rna_seq_test   <- gdsc_test[intersect(camptothecan_lines, rownames(gdsc_test)), ]
# 342 x 14242
cisplatin_rna_seq_train     <- gdsc_train[intersect(cisplatin_lines, rownames(gdsc_train)), ]
# 338 x 14242
cisplatin_rna_seq_test      <- gdsc_test[intersect(cisplatin_lines, rownames(gdsc_test)), ]
# 342 x 14242
cytarabine_rna_seq_train    <- gdsc_train[intersect(cytarabine_lines, rownames(gdsc_train)), ]
# 337 x 14242
cytarabine_rna_seq_test     <- gdsc_test[intersect(cytarabine_lines, rownames(gdsc_test)), ]
# 339 x 14242
dox_rna_seq_train           <- gdsc_train[intersect(dox_lines, rownames(gdsc_train)), ]
# 348 x 14242
dox_rna_seq_test            <- gdsc_test[intersect(dox_lines, rownames(gdsc_test)), ]
# 363 x 14242
etoposide_rna_seq_train     <- gdsc_train[intersect(etoposide_lines, rownames(gdsc_train)), ]
# 352 x 14242
etoposide_rna_seq_test      <- gdsc_test[intersect(etoposide_lines, rownames(gdsc_test)), ]
# 366 x 14242
gemcitabine_rna_seq_train   <- gdsc_train[intersect(gemcitabine_lines, rownames(gdsc_train)), ]
# 347 x 14242
gemcitabine_rna_seq_test    <- gdsc_test[intersect(gemcitabine_lines, rownames(gdsc_test)), ]
# 360 x 14242
methotrexate_rna_seq_train  <- gdsc_train[intersect(methotrexate_lines, rownames(gdsc_train)), ]
# 337 x 14242
methotrexate_rna_seq_test   <- gdsc_test[intersect(methotrexate_lines, rownames(gdsc_test)), ]
# 342 x 14242
mitomycin_rna_seq_train     <- gdsc_train[intersect(mitomycin_lines, rownames(gdsc_train)), ]
# 349 x 14242
mitomycin_rna_seq_test      <- gdsc_test[intersect(mitomycin_lines, rownames(gdsc_test)), ]
# 363 x 14242
sn38_rna_seq_train          <- gdsc_train[intersect(sn38_lines, rownames(gdsc_train)), ]
# 379 x 14242
sn38_rna_seq_test           <- gdsc_test[intersect(sn38_lines, rownames(gdsc_test)), ]
# 382 x 14242
temozolomide_rna_seq_train  <- gdsc_train[intersect(temozolomide_lines, rownames(gdsc_train)), ]
# 374 x 14242
temozolomide_rna_seq_test   <- gdsc_test[intersect(temozolomide_lines, rownames(gdsc_test)), ]
# 365 x 14242

# split clinical data
bleomycin_train        <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_train)), ]
bleomycin_test         <- bleomycin[which(bleomycin$COSMIC_ID %in% rownames(bleomycin_rna_seq_test)), ]

camptothecan_train     <- camptothecan[which(camptothecan$COSMIC_ID %in% rownames(camptothecan_rna_seq_train)), ]
camptothecan_test      <- camptothecan[which(camptothecan$COSMIC_ID %in% rownames(camptothecan_rna_seq_test)), ]

cisplatin_train        <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_train)), ]
cisplatin_test         <- cisplatin[which(cisplatin$COSMIC_ID %in% rownames(cisplatin_rna_seq_test)), ]

cytarabine_train       <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_train)), ]
cytarabine_test        <- cytarabine[which(cytarabine$COSMIC_ID %in% rownames(cytarabine_rna_seq_test)), ]

dox_train              <- dox[which(dox$COSMIC_ID %in% rownames(dox_rna_seq_train)), ]
dox_test               <- dox[which(dox$COSMIC_ID %in% rownames(dox_rna_seq_test)), ]

etoposide_train        <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_train)), ]
etoposide_test         <- etoposide[which(etoposide$COSMIC_ID %in% rownames(etoposide_rna_seq_test)), ]

gemcitabine_train      <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_train)), ]
gemcitabine_test       <- gemcitabine[which(gemcitabine$COSMIC_ID %in% rownames(gemcitabine_rna_seq_test)), ]

methotrexate_train     <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_train)), ]
methotrexate_test      <- methotrexate[which(methotrexate$COSMIC_ID %in% rownames(methotrexate_rna_seq_test)), ]

mitomycin_train        <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_train)), ]
mitomycin_test         <- mitomycin[which(mitomycin$COSMIC_ID %in% rownames(mitomycin_rna_seq_test)), ]

sn38_train             <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_train)), ]
sn38_test              <- sn38[which(sn38$COSMIC_ID %in% rownames(sn38_rna_seq_test)), ]

temozolomide_train     <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_train)), ]
temozolomide_test      <- temozolomide[which(temozolomide$COSMIC_ID %in% rownames(temozolomide_rna_seq_test)), ]

# scale data
bleomycin_rna_seq_train_scaled          <- apply(bleomycin_rna_seq_train, 2, scale)
bleomycin_rna_seq_test_scaled           <- as.data.frame(apply(bleomycin_rna_seq_test, 2, scale))

camptothecan_rna_seq_train_scaled       <- apply(camptothecan_rna_seq_train, 2, scale)
camptothecan_rna_seq_test_scaled        <- as.data.frame(apply(camptothecan_rna_seq_test, 2, scale))

cisplatin_rna_seq_train_scaled          <- apply(cisplatin_rna_seq_train, 2, scale)
cisplatin_rna_seq_test_scaled           <- as.data.frame(apply(cisplatin_rna_seq_test, 2, scale))

cytarabine_rna_seq_train_scaled         <- apply(cytarabine_rna_seq_train, 2, scale)
cytarabine_rna_seq_test_scaled          <- as.data.frame(apply(cytarabine_rna_seq_test, 2, scale))

dox_rna_seq_train_scaled                <- apply(dox_rna_seq_train, 2, scale)
dox_rna_seq_test_scaled                 <- as.data.frame(apply(dox_rna_seq_test, 2, scale))

etoposide_rna_seq_train_scaled          <- apply(etoposide_rna_seq_train, 2, scale)
etoposide_rna_seq_test_scaled           <- as.data.frame(apply(etoposide_rna_seq_test, 2, scale))

gemcitabine_rna_seq_train_scaled        <- apply(gemcitabine_rna_seq_train, 2, scale)
gemcitabine_rna_seq_test_scaled         <- as.data.frame(apply(gemcitabine_rna_seq_test, 2, scale))

methotrexate_rna_seq_train_scaled       <- apply(methotrexate_rna_seq_train, 2, scale)
methotrexate_rna_seq_test_scaled        <- as.data.frame(apply(methotrexate_rna_seq_test, 2, scale))

mitomycin_rna_seq_train_scaled          <- apply(mitomycin_rna_seq_train, 2, scale)
mitomycin_rna_seq_test_scaled           <- as.data.frame(apply(mitomycin_rna_seq_test, 2, scale))

sn38_rna_seq_train_scaled               <- apply(sn38_rna_seq_train, 2, scale)
sn38_rna_seq_test_scaled                <- as.data.frame(apply(sn38_rna_seq_test, 2, scale))

temozolomide_rna_seq_train_scaled       <- apply(temozolomide_rna_seq_train, 2, scale)
temozolomide_rna_seq_test_scaled        <- as.data.frame(apply(temozolomide_rna_seq_test, 2, scale))

### fit models --------
## BLEOMYCIN
# create model for predicting most sensitive samples 
bleo_most_fit_elnet <- cv.glmnet(x = as.matrix(bleomycin_rna_seq_train_scaled), y = bleomycin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')
# plot image of auc values to file
png(filename = 'Images/bleomycin_most_auc.png')
plot(bleo_most_fit_elnet)
dev.off()
# how well does this predict in a pan-cancer fashion on testing set?
new_bleo_most_sensitive_min <- predict(bleo_most_fit_elnet, newx = as.matrix(bleomycin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', type = 'response')
# get AUC values
bleo_most_test_gdsc_auc_min <- auc(bleomycin_test$most_sensitive, new_bleo_most_sensitive_min)
bleo_most_test_gdsc_auc_min <- round(bleo_most_test_gdsc_auc_min, digits = 2)

new_bleo_most_sensitive_1se <- predict(bleo_most_fit_elnet, newx = as.matrix(bleomycin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', type = 'response')
# get AUC values
bleo_most_test_gdsc_auc_1se <- auc(bleomycin_test$most_sensitive, new_bleo_most_sensitive_1se)
bleo_most_test_gdsc_auc_1se <- round(bleo_most_test_gdsc_auc_1se, digits = 2)

#save models in GLM_Models folder when complete and then program them to load
saveRDS(file = 'GLM_Models/bleomycin_most_model.rds', bleo_most_fit_elnet)

# create model for predicting least sensitive samples
bleo_least_fit_elnet <- cv.glmnet(x = as.matrix(bleomycin_rna_seq_train_scaled), y = bleomycin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/bleomycin_least_auc.png')
plot(bleo_least_fit_elnet)
dev.off()

new_bleo_least_sensitive_min <- predict(bleo_least_fit_elnet, newx = as.matrix(bleomycin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

bleo_least_test_gdsc_auc_min <- auc(bleomycin_test$least_sensitive, new_bleo_least_sensitive_min)
bleo_least_test_gdsc_auc_min <- round(bleo_least_test_gdsc_auc_min, digits = 2)

new_bleo_least_sensitive_1se <- predict(bleo_least_fit_elnet, newx = as.matrix(bleomycin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

bleo_least_test_gdsc_auc_1se <- auc(bleomycin_test$least_sensitive, new_bleo_least_sensitive_1se)
bleo_least_test_gdsc_auc_1se <- round(bleo_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/bleomycin_least_model.rds', bleo_least_fit_elnet)

## CAMPTOTHECAN
campto_most_fit_elnet <- cv.glmnet(x = as.matrix(camptothecan_rna_seq_train_scaled), y = camptothecan_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/camptothecan_most_auc.png')
plot(campto_most_fit_elnet)
dev.off()

new_campto_most_sensitive_min <- predict(campto_most_fit_elnet, newx = as.matrix(camptothecan_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

campto_most_test_gdsc_auc_min <- auc(camptothecan_test$most_sensitive, new_campto_most_sensitive_min)
campto_most_test_gdsc_auc_min <- round(campto_most_test_gdsc_auc_min, digits = 2)

new_campto_most_sensitive_1se <- predict(campto_most_fit_elnet, newx = as.matrix(camptothecan_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

campto_most_test_gdsc_auc_1se <- auc(camptothecan_test$most_sensitive, new_campto_most_sensitive_1se)
campto_most_test_gdsc_auc_1se <- round(campto_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/camptothecan_most_model.rds', campto_most_fit_elnet)

campto_least_fit_elnet <- cv.glmnet(x = as.matrix(camptothecan_rna_seq_train_scaled), y = camptothecan_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/camptothecan_least_auc.png')
plot(campto_least_fit_elnet)
dev.off()

new_campto_least_sensitive_min <- predict(campto_least_fit_elnet, newx = as.matrix(camptothecan_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

campto_least_test_gdsc_auc_min <- auc(camptothecan_test$least_sensitive, new_campto_least_sensitive_min)
campto_least_test_gdsc_auc_min <- round(campto_least_test_gdsc_auc_min, digits = 2)

new_campto_least_sensitive_1se <- predict(campto_least_fit_elnet, newx = as.matrix(camptothecan_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

campto_least_test_gdsc_auc_1se <- auc(camptothecan_test$least_sensitive, new_campto_least_sensitive_1se)
campto_least_test_gdsc_auc_1se <- round(campto_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/camptothecan_least_model.rds', campto_least_fit_elnet)

## CISPLATIN
cisplatin_most_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rna_seq_train_scaled), y = cisplatin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cisplatin_most_auc.png')
plot(cisplatin_most_fit_elnet)
dev.off()

new_cisplatin_most_sensitive_min <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_most_test_gdsc_auc_min <- auc(cisplatin_test$most_sensitive, new_cisplatin_most_sensitive_min)
cisplatin_most_test_gdsc_auc_min <- round(cisplatin_most_test_gdsc_auc_min, digits = 2)

new_cisplatin_most_sensitive_1se <- predict(cisplatin_most_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_most_test_gdsc_auc_1se <- auc(cisplatin_test$most_sensitive, new_cisplatin_most_sensitive_1se)
cisplatin_most_test_gdsc_auc_1se <- round(cisplatin_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/cisplatin_most_model.rds', cisplatin_most_fit_elnet)

cisplatin_least_fit_elnet <- cv.glmnet(x = as.matrix(cisplatin_rna_seq_train_scaled), y = cisplatin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cisplatin_least_auc.png')
plot(cisplatin_least_fit_elnet)
dev.off()

new_cisplatin_least_sensitive_min <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_least_test_gdsc_auc_min <- auc(cisplatin_test$least_sensitive, new_cisplatin_least_sensitive_min)
cisplatin_least_test_gdsc_auc_min <- round(cisplatin_least_test_gdsc_auc_min, digits = 2)

new_cisplatin_least_sensitive_1se <- predict(cisplatin_least_fit_elnet, newx = as.matrix(cisplatin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

cisplatin_least_test_gdsc_auc_1se <- auc(cisplatin_test$least_sensitive, new_cisplatin_least_sensitive_1se)
cisplatin_least_test_gdsc_auc_1se <- round(cisplatin_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/cisplatin_least_model.rds', cisplatin_least_fit_elnet)

## CYTARABINE
cytarabine_most_fit_elnet <- cv.glmnet(x = as.matrix(cytarabine_rna_seq_train_scaled), y = cytarabine_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cytarabine_most_auc.png')
plot(cytarabine_most_fit_elnet)
dev.off()

new_cytarabine_most_sensitive_min <- predict(cytarabine_most_fit_elnet, newx = as.matrix(cytarabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = TRUE, type = 'response')

cytarabine_most_test_gdsc_auc_min <- auc(cytarabine_test$most_sensitive, new_cytarabine_most_sensitive_min)
cytarabine_most_test_gdsc_auc_min <- round(cytarabine_most_test_gdsc_auc_min, digits = 2)

new_cytarabine_most_sensitive_1se <- predict(cytarabine_most_fit_elnet, newx = as.matrix(cytarabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = TRUE, type = 'response')

cytarabine_most_test_gdsc_auc_1se <- auc(cytarabine_test$most_sensitive, new_cytarabine_most_sensitive_1se)
cytarabine_most_test_gdsc_auc_1se <- round(cytarabine_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/cytarabine_most_model.rds', cytarabine_most_fit_elnet)

cytarabine_least_fit_elnet <- cv.glmnet(x = as.matrix(cytarabine_rna_seq_train_scaled), y = cytarabine_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/cytarabine_least_auc.png')
plot(cytarabine_least_fit_elnet)
dev.off()

new_cytarabine_least_sensitive_min <- predict(cytarabine_least_fit_elnet, newx = as.matrix(cytarabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = TRUE, type = 'response')

cytarabine_least_test_gdsc_auc_min <- auc(cytarabine_test$least_sensitive, new_cytarabine_least_sensitive_min)
cytarabine_least_test_gdsc_auc_min <- round(cytarabine_least_test_gdsc_auc_min, digits = 2)

new_cytarabine_least_sensitive_1se <- predict(cytarabine_least_fit_elnet, newx = as.matrix(cytarabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = TRUE, type = 'response')

cytarabine_least_test_gdsc_auc_1se <- auc(cytarabine_test$least_sensitive, new_cytarabine_least_sensitive_1se)
cytarabine_least_test_gdsc_auc_1se <- round(cytarabine_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/cytarabine_least_model.rds', cytarabine_least_fit_elnet)

## DOXORUBICIN
dox_most_fit_elnet <- cv.glmnet(x = as.matrix(dox_rna_seq_train_scaled), y = dox_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/dox_most_auc.png')
plot(dox_most_fit_elnet)
dev.off()

new_dox_most_sensitive_min <- predict(dox_most_fit_elnet, newx = as.matrix(dox_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

dox_most_test_gdsc_auc_min <- auc(dox_test$most_sensitive, new_dox_most_sensitive_min)
dox_most_test_gdsc_auc_min <- round(dox_most_test_gdsc_auc_min, digits = 2)

new_dox_most_sensitive_1se <- predict(dox_most_fit_elnet, newx = as.matrix(dox_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

dox_most_test_gdsc_auc_1se <- auc(dox_test$most_sensitive, new_dox_most_sensitive_1se)
dox_most_test_gdsc_auc_1se <- round(dox_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/dox_most_model.rds', dox_most_fit_elnet)

dox_least_fit_elnet <- cv.glmnet(x = as.matrix(dox_rna_seq_train_scaled), y = dox_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/dox_least_auc.png')
plot(dox_least_fit_elnet)
dev.off()

new_dox_least_sensitive_min <- predict(dox_least_fit_elnet, newx = as.matrix(dox_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

dox_least_test_gdsc_auc_min <- auc(dox_test$least_sensitive, new_dox_least_sensitive_min)
dox_least_test_gdsc_auc_min <- round(dox_least_test_gdsc_auc_min, digits = 2)

new_dox_least_sensitive_1se <- predict(dox_least_fit_elnet, newx = as.matrix(dox_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

dox_least_test_gdsc_auc_1se <- auc(dox_test$least_sensitive, new_dox_least_sensitive_1se)
dox_least_test_gdsc_auc_1se <- round(dox_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/dox_least_model.rds', dox_least_fit_elnet)

## ETOPOSIDE
etoposide_most_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rna_seq_train_scaled), y = etoposide_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/etoposide_most_auc.png')
plot(etoposide_most_fit_elnet)
dev.off()

new_etoposide_most_sensitive_min <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_most_test_gdsc_auc_min <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_min)
etoposide_most_test_gdsc_auc_min <- round(etoposide_most_test_gdsc_auc_min, digits = 2)

new_etoposide_most_sensitive_1se <- predict(etoposide_most_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_most_test_gdsc_auc_1se <- auc(etoposide_test$most_sensitive, new_etoposide_most_sensitive_1se)
etoposide_most_test_gdsc_auc_1se <- round(etoposide_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/etoposide_most_model.rds', etoposide_most_fit_elnet)

etoposide_least_fit_elnet <- cv.glmnet(x = as.matrix(etoposide_rna_seq_train_scaled), y = etoposide_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/etoposide_least_auc.png')
plot(etoposide_least_fit_elnet)
dev.off()

new_etoposide_least_sensitive_min <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_least_test_gdsc_auc_min <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_min)
etoposide_least_test_gdsc_auc_min <- round(etoposide_least_test_gdsc_auc_min, digits = 2)

new_etoposide_least_sensitive_1se <- predict(etoposide_least_fit_elnet, newx = as.matrix(etoposide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

etoposide_least_test_gdsc_auc_1se <- auc(etoposide_test$least_sensitive, new_etoposide_least_sensitive_1se)
etoposide_least_test_gdsc_auc_1se <- round(etoposide_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/etoposide_least_model.rds', etoposide_least_fit_elnet)

## GEMCITABINE
gemcitabine_most_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/gemcitabine_most_auc.png')
plot(gemcitabine_most_fit_elnet)
dev.off()

new_gemcitabine_most_sensitive_min <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_most_test_gdsc_auc_min <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_min)
gemcitabine_most_test_gdsc_auc_min <- round(gemcitabine_most_test_gdsc_auc_min, digits = 2)

new_gemcitabine_most_sensitive_1se <- predict(gemcitabine_most_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_most_test_gdsc_auc_1se <- auc(gemcitabine_test$most_sensitive, new_gemcitabine_most_sensitive_1se)
gemcitabine_most_test_gdsc_auc_1se <- round(gemcitabine_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/gemcitabine_most_model.rds', gemcitabine_most_fit_elnet)

gemcitabine_least_fit_elnet <- cv.glmnet(x = as.matrix(gemcitabine_rna_seq_train_scaled), y = gemcitabine_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/gemcitabine_least_auc.png')
plot(gemcitabine_least_fit_elnet)
dev.off()

new_gemcitabine_least_sensitive_min <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_least_test_gdsc_auc_min <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_min)
gemcitabine_least_test_gdsc_auc_min <- round(gemcitabine_least_test_gdsc_auc_min, digits = 2)

new_gemcitabine_least_sensitive_1se <- predict(gemcitabine_least_fit_elnet, newx = as.matrix(gemcitabine_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

gemcitabine_least_test_gdsc_auc_1se <- auc(gemcitabine_test$least_sensitive, new_gemcitabine_least_sensitive_1se)
gemcitabine_least_test_gdsc_auc_1se <- round(gemcitabine_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/gemcitabine_least_model.rds', gemcitabine_least_fit_elnet)

## METHOTREXATE
methotrexate_most_fit_elnet <- cv.glmnet(x = as.matrix(methotrexate_rna_seq_train_scaled), y = methotrexate_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/methotrexate_most_auc.png')
plot(methotrexate_most_fit_elnet)
dev.off()

new_methotrexate_most_sensitive_min <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_most_test_gdsc_auc_min <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_min)
methotrexate_most_test_gdsc_auc_min <- round(methotrexate_most_test_gdsc_auc_min, digits = 2)

new_methotrexate_most_sensitive_1se <- predict(methotrexate_most_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_most_test_gdsc_auc_1se <- auc(methotrexate_test$most_sensitive, new_methotrexate_most_sensitive_1se)
methotrexate_most_test_gdsc_auc_1se <- round(methotrexate_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/methotrexate_most_model.rds', methotrexate_most_fit_elnet)

methotrexate_least_fit_elnet <- cv.glmnet(x = as.matrix(methotrexate_rna_seq_train_scaled), y = methotrexate_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/methotrexate_least_auc.png')
plot(methotrexate_least_fit_elnet)
dev.off()

new_methotrexate_least_sensitive_min <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_least_test_gdsc_auc_min <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_min)
methotrexate_least_test_gdsc_auc_min <- round(methotrexate_least_test_gdsc_auc_min, digits = 2)

new_methotrexate_least_sensitive_1se <- predict(methotrexate_least_fit_elnet, newx = as.matrix(methotrexate_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

methotrexate_least_test_gdsc_auc_1se <- auc(methotrexate_test$least_sensitive, new_methotrexate_least_sensitive_1se)
methotrexate_least_test_gdsc_auc_1se <- round(methotrexate_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/methotrexate_least_model.rds', methotrexate_least_fit_elnet)

## MITOMYCIN
mitomycin_most_fit_elnet <- cv.glmnet(x = as.matrix(mitomycin_rna_seq_train_scaled), y = mitomycin_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/mitomycin_most_auc.png')
plot(mitomycin_most_fit_elnet)
dev.off()

new_mitomycin_most_sensitive_min <- predict(mitomycin_most_fit_elnet, newx = as.matrix(mitomycin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

mitomycin_most_test_gdsc_auc_min <- auc(mitomycin_test$most_sensitive, new_mitomycin_most_sensitive_min)
mitomycin_most_test_gdsc_auc_min <- round(mitomycin_most_test_gdsc_auc_min, digits = 2)

new_mitomycin_most_sensitive_1se <- predict(mitomycin_most_fit_elnet, newx = as.matrix(mitomycin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

mitomycin_most_test_gdsc_auc_1se <- auc(mitomycin_test$most_sensitive, new_mitomycin_most_sensitive_1se)
mitomycin_most_test_gdsc_auc_1se <- round(mitomycin_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/mitomycin_most_model.rds', mitomycin_most_fit_elnet)

mitomycin_least_fit_elnet <- cv.glmnet(x = as.matrix(mitomycin_rna_seq_train_scaled), y = mitomycin_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/mitomycin_least_auc.png')
plot(mitomycin_least_fit_elnet)
dev.off()

new_mitomycin_least_sensitive_min <- predict(mitomycin_least_fit_elnet, newx = as.matrix(mitomycin_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

mitomycin_least_test_gdsc_auc_min <- auc(mitomycin_test$least_sensitive, new_mitomycin_least_sensitive_min)
mitomycin_least_test_gdsc_auc_min <- round(mitomycin_least_test_gdsc_auc_min, digits = 2)

new_mitomycin_least_sensitive_1se <- predict(mitomycin_least_fit_elnet, newx = as.matrix(mitomycin_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

mitomycin_least_test_gdsc_auc_1se <- auc(mitomycin_test$least_sensitive, new_mitomycin_least_sensitive_1se)
mitomycin_least_test_gdsc_auc_1se <- round(mitomycin_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/mitomycin_least_model.rds', mitomycin_least_fit_elnet)

## SN38
sn38_most_fit_elnet <- cv.glmnet(x = as.matrix(sn38_rna_seq_train_scaled), y = sn38_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/sn38_most_auc.png')
plot(sn38_most_fit_elnet)
dev.off()

new_sn38_most_sensitive_min <- predict(sn38_most_fit_elnet, newx = as.matrix(sn38_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

sn38_most_test_gdsc_auc_min <- auc(sn38_test$most_sensitive, new_sn38_most_sensitive_min)
sn38_most_test_gdsc_auc_min <- round(sn38_most_test_gdsc_auc_min, digits = 2)

new_sn38_most_sensitive_1se <- predict(sn38_most_fit_elnet, newx = as.matrix(sn38_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

sn38_most_test_gdsc_auc_1se <- auc(sn38_test$most_sensitive, new_sn38_most_sensitive_1se)
sn38_most_test_gdsc_auc_1se <- round(sn38_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/sn38_most_model.rds', sn38_most_fit_elnet)

sn38_least_fit_elnet <- cv.glmnet(x = as.matrix(sn38_rna_seq_train_scaled), y = sn38_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/sn38_least_auc.png')
plot(sn38_least_fit_elnet)
dev.off()

new_sn38_least_sensitive_min <- predict(sn38_least_fit_elnet, newx = as.matrix(sn38_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

sn38_least_test_gdsc_auc_min <- auc(sn38_test$least_sensitive, new_sn38_least_sensitive_min)
sn38_least_test_gdsc_auc_min <- round(sn38_least_test_gdsc_auc_min, digits = 2)

new_sn38_least_sensitive_1se <- predict(sn38_least_fit_elnet, newx = as.matrix(sn38_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

sn38_least_test_gdsc_auc_1se <- auc(sn38_test$least_sensitive, new_sn38_least_sensitive_1se)
sn38_least_test_gdsc_auc_1se <- round(sn38_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/sn38_least_model.rds', sn38_least_fit_elnet)

## TEMOZOLOMIDE
temozolomide_most_fit_elnet <- cv.glmnet(x = as.matrix(temozolomide_rna_seq_train_scaled), y = temozolomide_train$most_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/temozolomide_most_auc.png')
plot(temozolomide_most_fit_elnet)
dev.off()

new_temozolomide_most_sensitive_min <- predict(temozolomide_most_fit_elnet, newx = as.matrix(temozolomide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

temozolomide_most_test_gdsc_auc_min <- auc(temozolomide_test$most_sensitive, new_temozolomide_most_sensitive_min)
temozolomide_most_test_gdsc_auc_min <- round(temozolomide_most_test_gdsc_auc_min, digits = 2)

new_temozolomide_most_sensitive_1se <- predict(temozolomide_most_fit_elnet, newx = as.matrix(temozolomide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

temozolomide_most_test_gdsc_auc_1se <- auc(temozolomide_test$most_sensitive, new_temozolomide_most_sensitive_1se)
temozolomide_most_test_gdsc_auc_1se <- round(temozolomide_most_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/temozolomide_most_model.rds', temozolomide_most_fit_elnet)

temozolomide_least_fit_elnet <- cv.glmnet(x = as.matrix(temozolomide_rna_seq_train_scaled), y = temozolomide_train$least_sensitive, family = 'binomial', alpha = 0.5, type.measure = 'auc')

png(filename = 'Images/temozolomide_least_auc.png')
plot(temozolomide_least_fit_elnet)
dev.off()

new_temozolomide_least_sensitive_min <- predict(temozolomide_least_fit_elnet, newx = as.matrix(temozolomide_rna_seq_test_scaled), s = 'lambda.min', interval = 'confidence', probability = FALSE, type = 'response')

temozolomide_least_test_gdsc_auc_min <- auc(temozolomide_test$least_sensitive, new_temozolomide_least_sensitive_min)
temozolomide_least_test_gdsc_auc_min <- round(temozolomide_least_test_gdsc_auc_min, digits = 2)

new_temozolomide_least_sensitive_1se <- predict(temozolomide_least_fit_elnet, newx = as.matrix(temozolomide_rna_seq_test_scaled), s = 'lambda.1se', interval = 'confidence', probability = FALSE, type = 'response')

temozolomide_least_test_gdsc_auc_1se <- auc(temozolomide_test$least_sensitive, new_temozolomide_least_sensitive_1se)
temozolomide_least_test_gdsc_auc_1se <- round(temozolomide_least_test_gdsc_auc_1se, digits = 2)

saveRDS(file = 'GLM_Models/temozolomide_least_model.rds', temozolomide_least_fit_elnet)

#### capture genes used in models --------
### BLEOMYCIN MOST SENSITIVE
bleo_most_min_tmp_coeffs  <- coef(bleo_most_fit_elnet, s = "lambda.min")
bleo_most_model_min       <- data.frame(name = bleo_most_min_tmp_coeffs@Dimnames[[1]][bleo_most_min_tmp_coeffs@i + 1], coefficient = bleo_most_min_tmp_coeffs@x)
write.csv(bleo_most_model_min, file = 'GLM_Models/bleo_most_model_min.csv', row.names = FALSE)

bleo_most_1se_tmp_coeffs  <- coef(bleo_most_fit_elnet, s = "lambda.1se")
bleo_most_model_1se       <- data.frame(name = bleo_most_1se_tmp_coeffs@Dimnames[[1]][bleo_most_1se_tmp_coeffs@i + 1], coefficient = bleo_most_1se_tmp_coeffs@x)
write.csv(bleo_most_model_1se, file = 'GLM_Models/bleo_most_model_1se.csv', row.names = FALSE)

## BLEOMYCIN LEAST SENSITIVE
bleo_least_min_tmp_coeffs <- coef(bleo_least_fit_elnet, s = "lambda.min")
bleo_least_model_min      <- data.frame(name = bleo_least_min_tmp_coeffs@Dimnames[[1]][bleo_least_min_tmp_coeffs@i + 1], coefficient = bleo_least_min_tmp_coeffs@x)
write.csv(bleo_least_model_min, file = 'GLM_Models/bleo_least_model_min.csv', row.names = FALSE)

bleo_least_1se_tmp_coeffs <- coef(bleo_least_fit_elnet, s = "lambda.1se")
bleo_least_model_1se      <- data.frame(name = bleo_least_1se_tmp_coeffs@Dimnames[[1]][bleo_least_1se_tmp_coeffs@i + 1], coefficient = bleo_least_1se_tmp_coeffs@x)
write.csv(bleo_least_model_1se, file = 'GLM_Models/bleo_least_model_1se.csv', row.names = FALSE)

## CAMPTOTHECAN MOST SENISTIVE
campto_most_min_tmp_coeffs  <- coef(campto_most_fit_elnet, s = "lambda.min")
campto_most_model_min       <- data.frame(name = campto_most_min_tmp_coeffs@Dimnames[[1]][campto_most_min_tmp_coeffs@i + 1], coefficient = campto_most_min_tmp_coeffs@x)
write.csv(campto_most_model_min, file = 'GLM_Models/campto_most_model_min.csv', row.names = FALSE)

campto_most_1se_tmp_coeffs  <- coef(campto_most_fit_elnet, s = "lambda.1se")
campto_most_model_1se       <- data.frame(name = campto_most_1se_tmp_coeffs@Dimnames[[1]][campto_most_1se_tmp_coeffs@i + 1], coefficient = campto_most_1se_tmp_coeffs@x)
write.csv(campto_most_model_1se, file = 'GLM_Models/campto_most_model_1se.csv', row.names = FALSE)

## CAMPTOTHECAN LEAST SENSITIVE
campto_least_min_tmp_coeffs <- coef(campto_least_fit_elnet, s = "lambda.min")
campto_least_model_min      <- data.frame(name = campto_least_min_tmp_coeffs@Dimnames[[1]][campto_least_min_tmp_coeffs@i + 1], coefficient = campto_least_min_tmp_coeffs@x)
write.csv(campto_least_model_min, file = 'GLM_Models/campto_least_model_min.csv', row.names = FALSE)

campto_least_1se_tmp_coeffs <- coef(campto_least_fit_elnet, s = "lambda.1se")
campto_least_model_1se      <- data.frame(name = campto_least_1se_tmp_coeffs@Dimnames[[1]][campto_least_1se_tmp_coeffs@i + 1], coefficient = campto_least_1se_tmp_coeffs@x)
write.csv(campto_least_model_1se, file = 'GLM_Models/campto_least_model_1se.csv', row.names = FALSE)

## CISPLATIN MOST SENSITIVE
cisplatin_most_min_tmp_coeffs <- coef(cisplatin_most_fit_elnet, s = "lambda.min")
cisplatin_most_model_min      <- data.frame(name = cisplatin_most_min_tmp_coeffs@Dimnames[[1]][cisplatin_most_min_tmp_coeffs@i + 1], coefficient = cisplatin_most_min_tmp_coeffs@x)
write.csv(cisplatin_most_model_min, file = 'GLM_Models/cisplatin_most_model_min.csv', row.names = FALSE)

cisplatin_most_1se_tmp_coeffs <- coef(cisplatin_most_fit_elnet, s = "lambda.1se")
cisplatin_most_model_1se      <- data.frame(name = cisplatin_most_1se_tmp_coeffs@Dimnames[[1]][cisplatin_most_1se_tmp_coeffs@i + 1], coefficient = cisplatin_most_1se_tmp_coeffs@x)
write.csv(cisplatin_most_model_1se, file = 'GLM_Models/cisplatin_most_model_1se.csv', row.names = FALSE)

## CISPLATINE LEAST SENSITIVE
cisplatin_least_min_tmp_coeffs  <- coef(cisplatin_least_fit_elnet, s = "lambda.min")
cisplatin_least_model_min       <- data.frame(name = cisplatin_least_min_tmp_coeffs@Dimnames[[1]][cisplatin_least_min_tmp_coeffs@i + 1], coefficient = cisplatin_least_min_tmp_coeffs@x)
write.csv(cisplatin_least_model_min, file = 'GLM_Models/cisplatin_least_model_min.csv', row.names = FALSE)

cisplatin_least_1se_tmp_coeffs  <- coef(cisplatin_least_fit_elnet, s = "lambda.1se")
cisplatin_least_model_1se       <- data.frame(name = cisplatin_least_1se_tmp_coeffs@Dimnames[[1]][cisplatin_least_1se_tmp_coeffs@i + 1], coefficient = cisplatin_least_1se_tmp_coeffs@x)
write.csv(cisplatin_least_model_1se, file = 'GLM_Models/cisplatin_least_model_1se.csv', row.names = FALSE)

## CYTARABINE MOST SENSITIVE
cytarabine_most_min_tmp_coeffs  <- coef(cytarabine_most_fit_elnet, s = "lambda.min")
cytarabine_most_model_min       <- data.frame(name = cytarabine_most_min_tmp_coeffs@Dimnames[[1]][cytarabine_most_min_tmp_coeffs@i + 1], coefficient = cytarabine_most_min_tmp_coeffs@x)
write.csv(cytarabine_most_model_min, file = 'GLM_Models/cytarabine_most_model_min.csv', row.names = FALSE)

cytarabine_most_1se_tmp_coeffs  <- coef(cytarabine_most_fit_elnet, s = "lambda.1se")
cytarabine_most_model_1se       <- data.frame(name = cytarabine_most_1se_tmp_coeffs@Dimnames[[1]][cytarabine_most_1se_tmp_coeffs@i + 1], coefficient = cytarabine_most_1se_tmp_coeffs@x)
write.csv(cytarabine_most_model_1se, file = 'GLM_Models/cytarabine_most_model_1se.csv', row.names = FALSE)

## CYTARABINE LEAST SENSITIVE
cytarabine_least_min_tmp_coeffs <- coef(cytarabine_least_fit_elnet, s = "lambda.min")
cytarabine_least_model_min      <- data.frame(name = cytarabine_least_min_tmp_coeffs@Dimnames[[1]][cytarabine_least_min_tmp_coeffs@i + 1], coefficient = cytarabine_least_min_tmp_coeffs@x)
write.csv(cytarabine_least_model_min, file = 'GLM_Models/cytarabine_least_model_min.csv', row.names = FALSE)

cytarabine_least_1se_tmp_coeffs <- coef(cytarabine_least_fit_elnet, s = "lambda.1se")
cytarabine_least_model_1se      <- data.frame(name = cytarabine_least_1se_tmp_coeffs@Dimnames[[1]][cytarabine_least_1se_tmp_coeffs@i + 1], coefficient = cytarabine_least_1se_tmp_coeffs@x)
write.csv(cytarabine_least_model_1se, file = 'GLM_Models/cytarabine_least_model_1se.csv', row.names = FALSE)

## DOX MOST SENSITIVE
dox_most_min_tmp_coeffs <- coef(dox_most_fit_elnet, s = "lambda.min")
dox_most_model_min      <- data.frame(name = dox_most_min_tmp_coeffs@Dimnames[[1]][dox_most_min_tmp_coeffs@i + 1], coefficient = dox_most_min_tmp_coeffs@x)
write.csv(dox_most_model_min, file = 'GLM_Models/dox_most_model_min.csv', row.names = FALSE)

dox_most_1se_tmp_coeffs <- coef(dox_most_fit_elnet, s = "lambda.1se")
dox_most_model_1se      <- data.frame(name = dox_most_1se_tmp_coeffs@Dimnames[[1]][dox_most_1se_tmp_coeffs@i + 1], coefficient = dox_most_1se_tmp_coeffs@x)
write.csv(dox_most_model_1se, file = 'GLM_Models/dox_most_model_1se.csv', row.names = FALSE)

## DOX LEAST SENSITIVE
dox_least_min_tmp_coeffs  <- coef(dox_least_fit_elnet, s = "lambda.min")
dox_least_model_min       <- data.frame(name = dox_least_min_tmp_coeffs@Dimnames[[1]][dox_least_min_tmp_coeffs@i + 1], coefficient = dox_least_min_tmp_coeffs@x)
write.csv(dox_least_model_min, file = 'GLM_Models/dox_least_model_min.csv', row.names = FALSE)

dox_least_1se_tmp_coeffs  <- coef(dox_least_fit_elnet, s = "lambda.1se")
dox_least_model_1se       <- data.frame(name = dox_least_1se_tmp_coeffs@Dimnames[[1]][dox_least_1se_tmp_coeffs@i + 1], coefficient = dox_least_1se_tmp_coeffs@x)
write.csv(dox_least_model_1se, file = 'GLM_Models/dox_least_model_1se.csv', row.names = FALSE)

## ETOPOSIDE MOST SENSITIVE
etoposide_most_min_tmp_coeffs <- coef(etoposide_most_fit_elnet, s = "lambda.min")
etoposide_most_model_min      <- data.frame(name = etoposide_most_min_tmp_coeffs@Dimnames[[1]][etoposide_most_min_tmp_coeffs@i + 1], coefficient = etoposide_most_min_tmp_coeffs@x)
write.csv(etoposide_most_model_min, file = 'GLM_Models/etoposide_most_model_min.csv', row.names = FALSE)

etoposide_most_1se_tmp_coeffs <- coef(etoposide_most_fit_elnet, s = "lambda.1se")
etoposide_most_model_1se      <- data.frame(name = etoposide_most_1se_tmp_coeffs@Dimnames[[1]][etoposide_most_1se_tmp_coeffs@i + 1], coefficient = etoposide_most_1se_tmp_coeffs@x)
write.csv(etoposide_most_model_1se, file = 'GLM_Models/etoposide_most_model_1se.csv', row.names = FALSE)

## ETOPOSIDE LEAST SENSITIVE
etoposide_least_min_tmp_coeffs  <- coef(etoposide_least_fit_elnet, s = "lambda.min")
etoposide_least_model_min       <- data.frame(name = etoposide_least_min_tmp_coeffs@Dimnames[[1]][etoposide_least_min_tmp_coeffs@i + 1], coefficient = etoposide_least_min_tmp_coeffs@x)
write.csv(etoposide_least_model_min, file = 'GLM_Models/etoposide_least_model_min.csv', row.names = FALSE)

etoposide_least_1se_tmp_coeffs  <- coef(etoposide_least_fit_elnet, s = "lambda.1se")
etoposide_least_model_1se       <- data.frame(name = etoposide_least_1se_tmp_coeffs@Dimnames[[1]][etoposide_least_1se_tmp_coeffs@i + 1], coefficient = etoposide_least_1se_tmp_coeffs@x)
write.csv(etoposide_least_model_1se, file = 'GLM_Models/etoposide_least_model_1se.csv', row.names = FALSE)

## GEMCITABINE MOST SENSITIVE
gemcitabine_most_min_tmp_coeffs <- coef(gemcitabine_most_fit_elnet, s = "lambda.min")
gemcitabine_most_model_min      <- data.frame(name = gemcitabine_most_min_tmp_coeffs@Dimnames[[1]][gemcitabine_most_min_tmp_coeffs@i + 1], coefficient = gemcitabine_most_min_tmp_coeffs@x)
write.csv(gemcitabine_most_model_min, file = 'GLM_Models/gemcitabine_most_model_min.csv', row.names = FALSE)

gemcitabine_most_1se_tmp_coeffs <- coef(gemcitabine_most_fit_elnet, s = "lambda.1se")
gemcitabine_most_model_1se      <- data.frame(name = gemcitabine_most_1se_tmp_coeffs@Dimnames[[1]][gemcitabine_most_1se_tmp_coeffs@i + 1], coefficient = gemcitabine_most_1se_tmp_coeffs@x)
write.csv(gemcitabine_most_model_1se, file = 'GLM_Models/gemcitabine_most_model_1se.csv', row.names = FALSE)

## GEMCITABINE LEAST SENSITIVE
gemcitabine_least_min_tmp_coeffs  <- coef(gemcitabine_least_fit_elnet, s = "lambda.min")
gemcitabine_least_model_min       <- data.frame(name = gemcitabine_least_min_tmp_coeffs@Dimnames[[1]][gemcitabine_least_min_tmp_coeffs@i + 1], coefficient = gemcitabine_least_min_tmp_coeffs@x)
write.csv(gemcitabine_least_model_min, file = 'GLM_Models/gemcitabine_least_model_min.csv', row.names = FALSE)

gemcitabine_least_1se_tmp_coeffs  <- coef(gemcitabine_least_fit_elnet, s = "lambda.1se")
gemcitabine_least_model_1se       <- data.frame(name = gemcitabine_least_1se_tmp_coeffs@Dimnames[[1]][gemcitabine_least_1se_tmp_coeffs@i + 1], coefficient = gemcitabine_least_1se_tmp_coeffs@x)
write.csv(gemcitabine_least_model_1se, file = 'GLM_Models/gemcitabine_least_model_1se.csv', row.names = FALSE)

## METHOTREXATE MOST SENSITIVE
methotrexate_most_min_tmp_coeffs  <- coef(methotrexate_most_fit_elnet, s = "lambda.min")
methotrexate_most_model_min       <- data.frame(name = methotrexate_most_min_tmp_coeffs@Dimnames[[1]][methotrexate_most_min_tmp_coeffs@i + 1], coefficient = methotrexate_most_min_tmp_coeffs@x)
write.csv(methotrexate_most_model_min, file = 'GLM_Models/methotrexate_most_model_min.csv', row.names = FALSE)

methotrexate_most_1se_tmp_coeffs  <- coef(methotrexate_most_fit_elnet, s = "lambda.1se")
methotrexate_most_model_1se       <- data.frame(name = methotrexate_most_1se_tmp_coeffs@Dimnames[[1]][methotrexate_most_1se_tmp_coeffs@i + 1], coefficient = methotrexate_most_1se_tmp_coeffs@x)
write.csv(methotrexate_most_model_1se, file = 'GLM_Models/methotrexate_most_model_1se.csv', row.names = FALSE)

## METHOTREXATE LEAST SENSITIVE
methotrexate_least_min_tmp_coeffs <- coef(methotrexate_least_fit_elnet, s = "lambda.min")
methotrexate_least_model_min      <- data.frame(name = methotrexate_least_min_tmp_coeffs@Dimnames[[1]][methotrexate_least_min_tmp_coeffs@i + 1], coefficient = methotrexate_least_min_tmp_coeffs@x)
write.csv(methotrexate_least_model_min, file = 'GLM_Models/methotrexate_least_model_min.csv', row.names = FALSE)

methotrexate_least_1se_tmp_coeffs <- coef(methotrexate_least_fit_elnet, s = "lambda.1se")
methotrexate_least_model_1se      <- data.frame(name = methotrexate_least_1se_tmp_coeffs@Dimnames[[1]][methotrexate_least_1se_tmp_coeffs@i + 1], coefficient = methotrexate_least_1se_tmp_coeffs@x)
write.csv(methotrexate_least_model_1se, file = 'GLM_Models/methotrexate_least_model_1se.csv', row.names = FALSE)

## MITOMYCIN MOST SENSITIVE
mitomycin_most_min_tmp_coeffs <- coef(mitomycin_most_fit_elnet, s = "lambda.min")
mitomycin_most_model_min      <- data.frame(name = mitomycin_most_min_tmp_coeffs@Dimnames[[1]][mitomycin_most_min_tmp_coeffs@i + 1], coefficient = mitomycin_most_min_tmp_coeffs@x)
write.csv(mitomycin_most_model_min, file = 'GLM_Models/mitomycin_most_model_min.csv', row.names = FALSE)

mitomycin_most_1se_tmp_coeffs <- coef(mitomycin_most_fit_elnet, s = "lambda.1se")
mitomycin_most_model_1se      <- data.frame(name = mitomycin_most_1se_tmp_coeffs@Dimnames[[1]][mitomycin_most_1se_tmp_coeffs@i + 1], coefficient = mitomycin_most_1se_tmp_coeffs@x)
write.csv(mitomycin_most_model_1se, file = 'GLM_Models/mitomycin_most_model_1se.csv', row.names = FALSE)

## MITOMYCIN LEAST SENSITIVE
mitomycin_least_min_tmp_coeffs  <- coef(mitomycin_least_fit_elnet, s = "lambda.min")
mitomycin_least_model_min       <- data.frame(name = mitomycin_least_min_tmp_coeffs@Dimnames[[1]][mitomycin_least_min_tmp_coeffs@i + 1], coefficient = mitomycin_least_min_tmp_coeffs@x)
write.csv(mitomycin_least_model_min, file = 'GLM_Models/mitomycin_least_model_min.csv', row.names = FALSE)

mitomycin_least_1se_tmp_coeffs  <- coef(mitomycin_least_fit_elnet, s = "lambda.1se")
mitomycin_least_model_1se       <- data.frame(name = mitomycin_least_1se_tmp_coeffs@Dimnames[[1]][mitomycin_least_1se_tmp_coeffs@i + 1], coefficient = mitomycin_least_1se_tmp_coeffs@x)
write.csv(mitomycin_least_model_1se, file = 'GLM_Models/mitomycin_least_model_1se.csv', row.names = FALSE)

## SN38 MOST SENSITIVE
sn38_most_min_tmp_coeffs  <- coef(sn38_most_fit_elnet, s = "lambda.min")
sn38_most_model_min       <- data.frame(name = sn38_most_min_tmp_coeffs@Dimnames[[1]][sn38_most_min_tmp_coeffs@i + 1], coefficient = sn38_most_min_tmp_coeffs@x)
write.csv(sn38_most_model_min, file = 'GLM_Models/sn38_most_model_min.csv', row.names = FALSE)

sn38_most_1se_tmp_coeffs  <- coef(sn38_most_fit_elnet, s = "lambda.1se")
sn38_most_model_1se       <- data.frame(name = sn38_most_1se_tmp_coeffs@Dimnames[[1]][sn38_most_1se_tmp_coeffs@i + 1], coefficient = sn38_most_1se_tmp_coeffs@x)
write.csv(sn38_most_model_1se, file = 'GLM_Models/sn38_most_model_1se.csv', row.names = FALSE)

## SN38 LEAST SENSITIVE
sn38_least_min_tmp_coeffs <- coef(sn38_least_fit_elnet, s = "lambda.min")
sn38_least_model_min      <- data.frame(name = sn38_least_min_tmp_coeffs@Dimnames[[1]][sn38_least_min_tmp_coeffs@i + 1], coefficient = sn38_least_min_tmp_coeffs@x)
write.csv(sn38_least_model_min, file = 'GLM_Models/sn38_least_model_min.csv', row.names = FALSE)

sn38_least_1se_tmp_coeffs <- coef(sn38_least_fit_elnet, s = "lambda.1se")
sn38_least_model_1se      <- data.frame(name = sn38_least_1se_tmp_coeffs@Dimnames[[1]][sn38_least_1se_tmp_coeffs@i + 1], coefficient = sn38_least_1se_tmp_coeffs@x)
write.csv(sn38_least_model_1se, file = 'GLM_Models/sn38_least_model_1se.csv', row.names = FALSE)

## TEMOZOLOMIDE MOST SENSITIVE
temozolomide_most_min_tmp_coeffs  <- coef(temozolomide_most_fit_elnet, s = "lambda.min")
temozolomide_most_model_min       <- data.frame(name = temozolomide_most_min_tmp_coeffs@Dimnames[[1]][temozolomide_most_min_tmp_coeffs@i + 1], coefficient = temozolomide_most_min_tmp_coeffs@x)
write.csv(temozolomide_most_model_min, file = 'GLM_Models/temozolomide_most_model_min.csv', row.names = FALSE)

temozolomide_most_1se_tmp_coeffs  <- coef(temozolomide_most_fit_elnet, s = "lambda.1se")
temozolomide_most_model_1se       <- data.frame(name = temozolomide_most_1se_tmp_coeffs@Dimnames[[1]][temozolomide_most_1se_tmp_coeffs@i + 1], coefficient = temozolomide_most_1se_tmp_coeffs@x)
write.csv(temozolomide_most_model_1se, file = 'GLM_Models/temozolomide_most_model_1se.csv', row.names = FALSE)

## TEMOZOLOMIDE LEAST SENSITIVE
temozolomide_least_min_tmp_coeffs <- coef(temozolomide_least_fit_elnet, s = "lambda.min")
temozolomide_least_model_min      <- data.frame(name = temozolomide_least_min_tmp_coeffs@Dimnames[[1]][temozolomide_least_min_tmp_coeffs@i + 1], coefficient = temozolomide_least_min_tmp_coeffs@x)
write.csv(temozolomide_least_model_min, file = 'GLM_Models/temozolomide_least_model_min.csv', row.names = FALSE)

temozolomide_least_1se_tmp_coeffs <- coef(temozolomide_least_fit_elnet, s = "lambda.1se")
temozolomide_least_model_1se      <- data.frame(name = temozolomide_least_1se_tmp_coeffs@Dimnames[[1]][temozolomide_least_1se_tmp_coeffs@i + 1], coefficient = temozolomide_least_1se_tmp_coeffs@x)
write.csv(temozolomide_least_model_1se, file = 'GLM_Models/temozolomide_least_model_1se.csv', row.names = FALSE)

## create table containing number of genes in each model and overall AUC ----
# the number of genes used in min models
min_model_genes <- c(nrow(bleo_most_model_min), nrow(bleo_least_model_min), nrow(campto_most_model_min), nrow(campto_least_model_min), 
                     nrow(cisplatin_most_model_min), nrow(cisplatin_least_model_min), nrow(cytarabine_most_model_min), nrow(cytarabine_least_model_min), 
                     nrow(dox_most_model_min), nrow(dox_least_model_min), nrow(etoposide_most_model_min), nrow(etoposide_least_model_min), 
                     nrow(gemcitabine_most_model_min), nrow(gemcitabine_least_model_min), nrow(methotrexate_most_model_min), nrow(methotrexate_least_model_min), 
                     nrow(mitomycin_most_model_min), nrow(mitomycin_least_model_min), nrow(sn38_most_model_min), nrow(sn38_least_model_min), 
                     nrow(temozolomide_most_model_min), nrow(temozolomide_least_model_min))

# the number of genes used in 1se models
one_se_model_genes <- c(nrow(bleo_most_model_1se), nrow(bleo_least_model_1se), nrow(campto_most_model_1se), nrow(campto_least_model_1se), 
                     nrow(cisplatin_most_model_1se), nrow(cisplatin_least_model_1se), nrow(cytarabine_most_model_1se), nrow(cytarabine_least_model_1se), 
                     nrow(dox_most_model_1se), nrow(dox_least_model_1se), nrow(etoposide_most_model_1se), nrow(etoposide_least_model_1se), 
                     nrow(gemcitabine_most_model_1se), nrow(gemcitabine_least_model_1se), nrow(methotrexate_most_model_1se), nrow(methotrexate_least_model_1se), 
                     nrow(mitomycin_most_model_1se), nrow(mitomycin_least_model_1se), nrow(sn38_most_model_1se), nrow(sn38_least_model_1se), 
                     nrow(temozolomide_most_model_1se), nrow(temozolomide_least_model_1se))

# overall AUC values for min models
auc_min_models <- c(bleo_most_test_gdsc_auc_min, bleo_least_test_gdsc_auc_min, campto_most_test_gdsc_auc_min, campto_least_test_gdsc_auc_min, 
                    cisplatin_most_test_gdsc_auc_min, cisplatin_least_test_gdsc_auc_min, cytarabine_most_test_gdsc_auc_min, cytarabine_least_test_gdsc_auc_min, 
                    dox_most_test_gdsc_auc_min, dox_least_test_gdsc_auc_min, etoposide_most_test_gdsc_auc_min, etoposide_least_test_gdsc_auc_min, 
                    gemcitabine_most_test_gdsc_auc_min, gemcitabine_least_test_gdsc_auc_min, methotrexate_most_test_gdsc_auc_min, methotrexate_least_test_gdsc_auc_min, 
                    mitomycin_most_test_gdsc_auc_min, mitomycin_least_test_gdsc_auc_min, sn38_most_test_gdsc_auc_min, sn38_least_test_gdsc_auc_min, 
                    temozolomide_most_test_gdsc_auc_min, temozolomide_least_test_gdsc_auc_min)

# overall AUC values for 1se models
auc_1se_models <- c(bleo_most_test_gdsc_auc_1se, bleo_least_test_gdsc_auc_1se, campto_most_test_gdsc_auc_1se, campto_least_test_gdsc_auc_1se, 
                    cisplatin_most_test_gdsc_auc_1se, cisplatin_least_test_gdsc_auc_1se, cytarabine_most_test_gdsc_auc_1se, cytarabine_least_test_gdsc_auc_1se, 
                    dox_most_test_gdsc_auc_1se, dox_least_test_gdsc_auc_1se, etoposide_most_test_gdsc_auc_1se, etoposide_least_test_gdsc_auc_1se, 
                    gemcitabine_most_test_gdsc_auc_1se, gemcitabine_least_test_gdsc_auc_1se, methotrexate_most_test_gdsc_auc_1se, methotrexate_least_test_gdsc_auc_1se, 
                    mitomycin_most_test_gdsc_auc_1se, mitomycin_least_test_gdsc_auc_1se, sn38_most_test_gdsc_auc_1se, sn38_least_test_gdsc_auc_1se, 
                    temozolomide_most_test_gdsc_auc_1se, temozolomide_least_test_gdsc_auc_1se)

# put together
genes_and_auc_df <- data.frame(min_model_genes, auc_min_models, one_se_model_genes, auc_1se_models)
rownames(genes_and_auc_df) <- c('Bleomycin_most_sensitive', 'Bleomycin_least_senistive', 'Camptothecan_most_sensitive', 'Camptothecan_least_sensitive', 
                                'Cisplatin_most_sensitive', 'Cisplatin_least_sensitive', 'Cytarabine_most_sensitive', 'Cytarabine_least_sensitive', 
                                'Doxorubicin_most_sensitive', 'Doxorubicin_least_sensitive', 'Etoposide_most_sensitive', 'Etoposide_least_sensitive', 
                                'Gemcitabine_most_sensitive', 'Gemcitabine_least_sensitive', 'Methotrexate_most_sensitive', 'Methotrexate_least_sensitive', 
                                'Mitomycin_most_sensitive', 'Mitomycin_least_sensitive', 'Sn38_most_sensitive', 'Sn38_least_sensitive', 
                                'Temozolomide_most_sensitive', 'Temozolomide_least_sensitive')
colnames(genes_and_auc_df) <- c('Number_genes_min_model', 'AUC_min_model', 'Number_genes_1se_model', 'AUC_1se_model')
genes_and_auc_df$Number_genes_lost <- genes_and_auc_df$Number_genes_min_model - genes_and_auc_df$Number_genes_1se_model
# creating table
genes_auc_formattable <- formattable(genes_and_auc_df, list(
  AUC_min_model = color_tile('dodgerblue', 'firebrick1'),
  AUC_1se_model = color_tile('dodgerblue', 'firebrick1'),
  area(col= c(Number_genes_min_model, Number_genes_1se_model, Number_genes_lost)) ~ normalize_bar('springgreen')
))
 # and exporting it
webshot::install_phantomjs()
# from: https://github.com/renkun-ken/formattable/issues/26
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

export_formattable(genes_auc_formattable, file = 'Images/Genes_AUC_table.png')

### get auc for individual cancer types -------

# add TCGA classes to all testing data
# BLEOMYCIN
bleo_idx        <- gdsc_cell_line_info$COSMIC_ID %in% bleomycin_test$COSMIC_ID
gdsc_info_bleo  <- gdsc_cell_line_info[bleo_idx, ]
bleomycin_test  <- merge(bleomycin_test, gdsc_info_bleo, by = 'COSMIC_ID')

# CAMPTOTHECAN
campto_idx        <- gdsc_cell_line_info$COSMIC_ID %in% camptothecan_test$COSMIC_ID
gdsc_info_campto  <- gdsc_cell_line_info[campto_idx, ]
camptothecan_test <- merge(camptothecan_test, gdsc_info_campto, by = 'COSMIC_ID')

# CISPLATIN
cisplatin_idx       <- gdsc_cell_line_info$COSMIC_ID %in% cisplatin_test$COSMIC_ID
gdsc_info_cisplatin <- gdsc_cell_line_info[cisplatin_idx, ]
cisplatin_test      <- merge(cisplatin_test, gdsc_info_cisplatin, by = 'COSMIC_ID')

# CYTARABINE
cytarabine_idx        <- gdsc_cell_line_info$COSMIC_ID %in% cisplatin_test$COSMIC_ID
gdsc_info_cytarabine  <- gdsc_cell_line_info[cytarabine_idx, ]
cytarabine_test       <- merge(cytarabine_test, gdsc_info_cytarabine, by = 'COSMIC_ID')

# DOX
dox_idx       <- gdsc_cell_line_info$COSMIC_ID %in% dox_test$COSMIC_ID
gdsc_info_dox <- gdsc_cell_line_info[dox_idx, ]
dox_test      <- merge(dox_test, gdsc_info_dox, by = 'COSMIC_ID')

# ETOPOSIDE
etoposide_idx       <- gdsc_cell_line_info$COSMIC_ID %in% etoposide_test$COSMIC_ID
gdsc_info_etoposide <- gdsc_cell_line_info[etoposide_idx, ]
etoposide_test      <- merge(etoposide_test, gdsc_info_etoposide, by = 'COSMIC_ID')

# GEMCITABINE
gemcitabine_idx       <- gdsc_cell_line_info$COSMIC_ID %in% gemcitabine_test$COSMIC_ID
gdsc_info_gemcitabine <- gdsc_cell_line_info[gemcitabine_idx, ]
gemcitabine_test      <- merge(gemcitabine_test, gdsc_info_gemcitabine, by = 'COSMIC_ID')

# METHOTREXATE
methotrexate_idx        <- gdsc_cell_line_info$COSMIC_ID %in% methotrexate_test$COSMIC_ID
gdsc_info_methotrexate  <- gdsc_cell_line_info[methotrexate_idx, ]
methotrexate_test       <- merge(methotrexate_test, gdsc_info_methotrexate, by = 'COSMIC_ID')

# MITOMYCIN
mitomycin_idx       <- gdsc_cell_line_info$COSMIC_ID %in% mitomycin_test$COSMIC_ID
gdsc_info_mitomycin <- gdsc_cell_line_info[mitomycin_idx, ]
mitomycin_test      <- merge(mitomycin_test, gdsc_info_mitomycin, by = 'COSMIC_ID')

# SN38
sn38_idx        <- gdsc_cell_line_info$COSMIC_ID %in% sn38_test$COSMIC_ID
gdsc_info_sn38  <- gdsc_cell_line_info[sn38_idx, ]
sn38_test       <- merge(sn38_test, gdsc_info_sn38, by = 'COSMIC_ID')

# TEMOZOLOMIDE
temozolomide_idx        <- gdsc_cell_line_info$COSMIC_ID %in% temozolomide_test$COSMIC_ID
gdsc_info_temozolomide  <- gdsc_cell_line_info[temozolomide_idx, ]
temozolomide_test       <- merge(temozolomide_test, gdsc_info_temozolomide, by = 'COSMIC_ID')

# AUC by TCGA class
## BLEOMYCIN
table(bleomycin_test$TCGA_class)

# BLCA
bleo_blca_idx <- bleomycin_test$TCGA_class == 'BLCA' #7

bleo_blca_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_blca_idx]
actual_bleo_blca_most         <- bleomycin_test$most_sensitive[bleo_blca_idx]
bleo_blca_most_min_auc        <- auc(actual_bleo_blca_most, bleo_blca_most_sensitive_min)

bleo_blca_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_blca_idx]
actual_bleo_blca_least <- bleomycin_test$least_sensitive[bleo_blca_idx]
bleo_blca_least_min_auc <- auc(actual_bleo_blca_least, bleo_blca_least_sensitive_min)

bleo_blca_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_blca_idx]
actual_bleo_blca_most <- bleomycin_test$most_sensitive[bleo_blca_idx]
bleo_blca_most_1se_auc <- auc(actual_bleo_blca_most, bleo_blca_most_sensitive_1se)

bleo_blca_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_blca_idx]
actual_bleo_blca_least <- bleomycin_test$least_sensitive[bleo_blca_idx]
bleo_blca_least_1se_auc <- auc(actual_bleo_blca_least, bleo_blca_least_sensitive_1se)

# BRCA
bleo_brca_idx <- bleomycin_test$TCGA_class == 'BRCA' #29

bleo_brca_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_brca_idx]
actual_bleo_brca_most         <- bleomycin_test$most_sensitive[bleo_brca_idx]
bleo_brca_most_min_auc        <- auc(actual_bleo_brca_most, bleo_brca_most_sensitive_min)

bleo_brca_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_brca_idx]
actual_bleo_brca_least <- bleomycin_test$least_sensitive[bleo_brca_idx]
bleo_brca_least_min_auc <- auc(actual_bleo_brca_least, bleo_brca_least_sensitive_min)

bleo_brca_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_brca_idx]
actual_bleo_brca_most <- bleomycin_test$most_sensitive[bleo_brca_idx]
bleo_brca_most_1se_auc <- auc(actual_bleo_brca_most, bleo_brca_most_sensitive_1se)

bleo_brca_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_brca_idx]
actual_bleo_brca_least <- bleomycin_test$least_sensitive[bleo_brca_idx]
bleo_brca_least_1se_auc <- auc(actual_bleo_brca_least, bleo_brca_least_sensitive_1se)

# CESC
bleo_cesc_idx <- bleomycin_test$TCGA_class == 'CESC' #4

bleo_cesc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_cesc_idx]
actual_bleo_cesc_most         <- bleomycin_test$most_sensitive[bleo_cesc_idx]
bleo_cesc_most_min_auc        <- auc(actual_bleo_cesc_most, bleo_cesc_most_sensitive_min)

bleo_cesc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_cesc_idx]
actual_bleo_cesc_least <- bleomycin_test$least_sensitive[bleo_cesc_idx]
bleo_cesc_least_min_auc <- auc(actual_bleo_cesc_least, bleo_cesc_least_sensitive_min)

bleo_cesc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_cesc_idx]
actual_bleo_cesc_most <- bleomycin_test$most_sensitive[bleo_cesc_idx]
bleo_cesc_most_1se_auc <- auc(actual_bleo_cesc_most, bleo_cesc_most_sensitive_1se)

bleo_cesc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_cesc_idx]
actual_bleo_cesc_least <- bleomycin_test$least_sensitive[bleo_cesc_idx]
bleo_cesc_least_1se_auc <- auc(actual_bleo_cesc_least, bleo_cesc_least_sensitive_1se)

# COREAD
bleo_coread_idx <- bleomycin_test$TCGA_class == 'COREAD' #21

bleo_coread_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_coread_idx]
actual_bleo_coread_most         <- bleomycin_test$most_sensitive[bleo_coread_idx]
bleo_coread_most_min_auc        <- auc(actual_bleo_coread_most, bleo_coread_most_sensitive_min)

bleo_coread_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_coread_idx]
actual_bleo_coread_least <- bleomycin_test$least_sensitive[bleo_coread_idx]
bleo_coread_least_min_auc <- auc(actual_bleo_coread_least, bleo_coread_least_sensitive_min)

bleo_coread_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_coread_idx]
actual_bleo_coread_most <- bleomycin_test$most_sensitive[bleo_coread_idx]
bleo_coread_most_1se_auc <- auc(actual_bleo_coread_most, bleo_coread_most_sensitive_1se)

bleo_coread_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_coread_idx]
actual_bleo_coread_least <- bleomycin_test$least_sensitive[bleo_coread_idx]
bleo_coread_least_1se_auc <- auc(actual_bleo_coread_least, bleo_coread_least_sensitive_1se)

# ESCA
bleo_esca_idx <- bleomycin_test$TCGA_class == 'ESCA' #15

bleo_esca_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_esca_idx]
actual_bleo_esca_most         <- bleomycin_test$most_sensitive[bleo_esca_idx]
bleo_esca_most_min_auc        <- auc(actual_bleo_esca_most, bleo_esca_most_sensitive_min)

bleo_esca_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_esca_idx]
actual_bleo_esca_least <- bleomycin_test$least_sensitive[bleo_esca_idx]
bleo_esca_least_min_auc <- auc(actual_bleo_esca_least, bleo_esca_least_sensitive_min)

bleo_esca_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_esca_idx]
actual_bleo_esca_most <- bleomycin_test$most_sensitive[bleo_esca_idx]
bleo_esca_most_1se_auc <- auc(actual_bleo_esca_most, bleo_esca_most_sensitive_1se)

bleo_esca_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_esca_idx]
actual_bleo_esca_least <- bleomycin_test$least_sensitive[bleo_esca_idx]
bleo_esca_least_1se_auc <- auc(actual_bleo_esca_least, bleo_esca_least_sensitive_1se)

# GBM
bleo_gbm_idx <- bleomycin_test$TCGA_class == 'GBM' #19

bleo_gbm_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_gbm_idx]
actual_bleo_gbm_most         <- bleomycin_test$most_sensitive[bleo_gbm_idx]
bleo_gbm_most_min_auc        <- auc(actual_bleo_gbm_most, bleo_gbm_most_sensitive_min)

bleo_gbm_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_gbm_idx]
actual_bleo_gbm_least <- bleomycin_test$least_sensitive[bleo_gbm_idx]
bleo_gbm_least_min_auc <- auc(actual_bleo_gbm_least, bleo_gbm_least_sensitive_min)

bleo_gbm_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_gbm_idx]
actual_bleo_gbm_most <- bleomycin_test$most_sensitive[bleo_gbm_idx]
bleo_gbm_most_1se_auc <- auc(actual_bleo_gbm_most, bleo_gbm_most_sensitive_1se)

bleo_gbm_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_gbm_idx]
actual_bleo_gbm_least <- bleomycin_test$least_sensitive[bleo_gbm_idx]
bleo_gbm_least_1se_auc <- auc(actual_bleo_gbm_least, bleo_gbm_least_sensitive_1se)


# HNSC
bleo_hnsc_idx <- bleomycin_test$TCGA_class == 'HNSC' #22

bleo_hnsc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_hnsc_idx]
actual_bleo_hnsc_most         <- bleomycin_test$most_sensitive[bleo_hnsc_idx]
bleo_hnsc_most_min_auc        <- auc(actual_bleo_hnsc_most, bleo_hnsc_most_sensitive_min)

bleo_hnsc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_hnsc_idx]
actual_bleo_hnsc_least <- bleomycin_test$least_sensitive[bleo_hnsc_idx]
bleo_hnsc_least_min_auc <- auc(actual_bleo_hnsc_least, bleo_hnsc_least_sensitive_min)

bleo_hnsc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_hnsc_idx]
actual_bleo_hnsc_most <- bleomycin_test$most_sensitive[bleo_hnsc_idx]
bleo_hnsc_most_1se_auc <- auc(actual_bleo_hnsc_most, bleo_hnsc_most_sensitive_1se)

bleo_hnsc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_hnsc_idx]
actual_bleo_hnsc_least <- bleomycin_test$least_sensitive[bleo_hnsc_idx]
bleo_hnsc_least_1se_auc <- auc(actual_bleo_hnsc_least, bleo_hnsc_least_sensitive_1se)


# KIRC
bleo_kirc_idx <- bleomycin_test$TCGA_class == 'KIRC' #18

bleo_kirc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_kirc_idx]
actual_bleo_kirc_most         <- bleomycin_test$most_sensitive[bleo_kirc_idx]
bleo_kirc_most_min_auc        <- auc(actual_bleo_kirc_most, bleo_kirc_most_sensitive_min)

bleo_kirc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_kirc_idx]
actual_bleo_kirc_least <- bleomycin_test$least_sensitive[bleo_kirc_idx]
bleo_kirc_least_min_auc <- auc(actual_bleo_kirc_least, bleo_kirc_least_sensitive_min)

bleo_kirc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_kirc_idx]
actual_bleo_kirc_most <- bleomycin_test$most_sensitive[bleo_kirc_idx]
bleo_kirc_most_1se_auc <- auc(actual_bleo_kirc_most, bleo_kirc_most_sensitive_1se)

bleo_kirc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_kirc_idx]
actual_bleo_kirc_least <- bleomycin_test$least_sensitive[bleo_kirc_idx]
bleo_kirc_least_1se_auc <- auc(actual_bleo_kirc_least, bleo_kirc_least_sensitive_1se)

#LGG
bleo_lgg_idx <- bleomycin_test$TCGA_class == 'LGG' #9

bleo_lgg_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_lgg_idx]
actual_bleo_lgg_most         <- bleomycin_test$most_sensitive[bleo_lgg_idx]
bleo_lgg_most_min_auc        <- auc(actual_bleo_lgg_most, bleo_lgg_most_sensitive_min)

bleo_lgg_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_lgg_idx]
actual_bleo_lgg_least <- bleomycin_test$least_sensitive[bleo_lgg_idx]
bleo_lgg_least_min_auc <- auc(actual_bleo_lgg_least, bleo_lgg_least_sensitive_min)

bleo_lgg_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_lgg_idx]
actual_bleo_lgg_most <- bleomycin_test$most_sensitive[bleo_lgg_idx]
bleo_lgg_most_1se_auc <- auc(actual_bleo_lgg_most, bleo_lgg_most_sensitive_1se)

bleo_lgg_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_lgg_idx]
actual_bleo_lgg_least <- bleomycin_test$least_sensitive[bleo_lgg_idx]
bleo_lgg_least_1se_auc <- auc(actual_bleo_lgg_least, bleo_lgg_least_sensitive_1se)

# LIHC
bleo_lihc_idx <- bleomycin_test$TCGA_class == 'LIHC' #8

bleo_lihc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_lihc_idx]
actual_bleo_lihc_most         <- bleomycin_test$most_sensitive[bleo_lihc_idx]
bleo_lihc_most_min_auc        <- auc(actual_bleo_lihc_most, bleo_lihc_most_sensitive_min)

bleo_lihc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_lihc_idx]
actual_bleo_lihc_least <- bleomycin_test$least_sensitive[bleo_lihc_idx]
bleo_lihc_least_min_auc <- auc(actual_bleo_lihc_least, bleo_lihc_least_sensitive_min)

bleo_lihc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_lihc_idx]
actual_bleo_lihc_most <- bleomycin_test$most_sensitive[bleo_lihc_idx]
bleo_lihc_most_1se_auc <- auc(actual_bleo_lihc_most, bleo_lihc_most_sensitive_1se)

bleo_lihc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_lihc_idx]
actual_bleo_lihc_least <- bleomycin_test$least_sensitive[bleo_lihc_idx]
bleo_lihc_least_1se_auc <- auc(actual_bleo_lihc_least, bleo_lihc_least_sensitive_1se)

# LUAD
bleo_luad_idx <- bleomycin_test$TCGA_class == 'LUAD' #32

bleo_luad_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_luad_idx]
actual_bleo_luad_most         <- bleomycin_test$most_sensitive[bleo_luad_idx]
bleo_luad_most_min_auc        <- auc(actual_bleo_luad_most, bleo_luad_most_sensitive_min)

bleo_luad_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_luad_idx]
actual_bleo_luad_least <- bleomycin_test$least_sensitive[bleo_luad_idx]
bleo_luad_least_min_auc <- auc(actual_bleo_luad_least, bleo_luad_least_sensitive_min)

bleo_luad_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_luad_idx]
actual_bleo_luad_most <- bleomycin_test$most_sensitive[bleo_luad_idx]
bleo_luad_most_1se_auc <- auc(actual_bleo_luad_most, bleo_luad_most_sensitive_1se)

bleo_luad_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_luad_idx]
actual_bleo_luad_least <- bleomycin_test$least_sensitive[bleo_luad_idx]
bleo_luad_least_1se_auc <- auc(actual_bleo_luad_least, bleo_luad_least_sensitive_1se)

# LUSC
bleo_lusc_idx <- bleomycin_test$TCGA_class == 'LUSC' #6

bleo_lusc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_lusc_idx]
actual_bleo_lusc_most         <- bleomycin_test$most_sensitive[bleo_lusc_idx]
bleo_lusc_most_min_auc        <- auc(actual_bleo_lusc_most, bleo_lusc_most_sensitive_min)

bleo_lusc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_lusc_idx]
actual_bleo_lusc_least <- bleomycin_test$least_sensitive[bleo_lusc_idx]
bleo_lusc_least_min_auc <- auc(actual_bleo_lusc_least, bleo_lusc_least_sensitive_min)

bleo_lusc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_lusc_idx]
actual_bleo_lusc_most <- bleomycin_test$most_sensitive[bleo_lusc_idx]
bleo_lusc_most_1se_auc <- auc(actual_bleo_lusc_most, bleo_lusc_most_sensitive_1se)

bleo_lusc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_lusc_idx]
actual_bleo_lusc_least <- bleomycin_test$least_sensitive[bleo_lusc_idx]
bleo_lusc_least_1se_auc <- auc(actual_bleo_lusc_least, bleo_lusc_least_sensitive_1se)

# MB
bleo_mb_idx <- bleomycin_test$TCGA_class == 'MB' #3

bleo_mb_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_mb_idx]
actual_bleo_mb_most         <- bleomycin_test$most_sensitive[bleo_mb_idx]
bleo_mb_most_min_auc        <- auc(actual_bleo_mb_most, bleo_mb_most_sensitive_min)

bleo_mb_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_mb_idx]
actual_bleo_mb_least <- bleomycin_test$least_sensitive[bleo_mb_idx]
bleo_mb_least_min_auc <- auc(actual_bleo_mb_least, bleo_mb_least_sensitive_min)

bleo_mb_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_mb_idx]
actual_bleo_mb_most <- bleomycin_test$most_sensitive[bleo_mb_idx]
bleo_mb_most_1se_auc <- auc(actual_bleo_mb_most, bleo_mb_most_sensitive_1se)

bleo_mb_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_mb_idx]
actual_bleo_mb_least <- bleomycin_test$least_sensitive[bleo_mb_idx]
bleo_mb_least_1se_auc <- auc(actual_bleo_mb_least, bleo_mb_least_sensitive_1se)

# MESO
bleo_meso_idx <- bleomycin_test$TCGA_class == 'MESO' #7

bleo_meso_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_meso_idx]
actual_bleo_meso_most         <- bleomycin_test$most_sensitive[bleo_meso_idx]
bleo_meso_most_min_auc        <- auc(actual_bleo_meso_most, bleo_meso_most_sensitive_min)

bleo_meso_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_meso_idx]
actual_bleo_meso_least <- bleomycin_test$least_sensitive[bleo_meso_idx]
bleo_meso_least_min_auc <- auc(actual_bleo_meso_least, bleo_meso_least_sensitive_min)

bleo_meso_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_meso_idx]
actual_bleo_meso_most <- bleomycin_test$most_sensitive[bleo_meso_idx]
bleo_meso_most_1se_auc <- auc(actual_bleo_meso_most, bleo_meso_most_sensitive_1se)

bleo_meso_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_meso_idx]
actual_bleo_meso_least <- bleomycin_test$least_sensitive[bleo_meso_idx]
bleo_meso_least_1se_auc <- auc(actual_bleo_meso_least, bleo_meso_least_sensitive_1se)

# NB
bleo_nb_idx <- bleomycin_test$TCGA_class == 'NB' #14

bleo_nb_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_nb_idx]
actual_bleo_nb_most         <- bleomycin_test$most_sensitive[bleo_nb_idx]
bleo_nb_most_min_auc        <- auc(actual_bleo_nb_most, bleo_nb_most_sensitive_min)

bleo_nb_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_nb_idx]
actual_bleo_nb_least <- bleomycin_test$least_sensitive[bleo_nb_idx]
bleo_nb_least_min_auc <- auc(actual_bleo_nb_least, bleo_nb_least_sensitive_min)

bleo_nb_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_nb_idx]
actual_bleo_nb_most <- bleomycin_test$most_sensitive[bleo_nb_idx]
bleo_nb_most_1se_auc <- auc(actual_bleo_nb_most, bleo_nb_most_sensitive_1se)

bleo_nb_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_nb_idx]
actual_bleo_nb_least <- bleomycin_test$least_sensitive[bleo_nb_idx]
bleo_nb_least_1se_auc <- auc(actual_bleo_nb_least, bleo_nb_least_sensitive_1se)

# OV
bleo_ov_idx <- bleomycin_test$TCGA_class == 'OV' #18

bleo_ov_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_ov_idx]
actual_bleo_ov_most         <- bleomycin_test$most_sensitive[bleo_ov_idx]
bleo_ov_most_min_auc        <- auc(actual_bleo_ov_most, bleo_ov_most_sensitive_min)

bleo_ov_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_ov_idx]
actual_bleo_ov_least <- bleomycin_test$least_sensitive[bleo_ov_idx]
bleo_ov_least_min_auc <- auc(actual_bleo_ov_least, bleo_ov_least_sensitive_min)

bleo_ov_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_ov_idx]
actual_bleo_ov_most <- bleomycin_test$most_sensitive[bleo_ov_idx]
bleo_ov_most_1se_auc <- auc(actual_bleo_ov_most, bleo_ov_most_sensitive_1se)

bleo_ov_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_ov_idx]
actual_bleo_ov_least <- bleomycin_test$least_sensitive[bleo_ov_idx]
bleo_ov_least_1se_auc <- auc(actual_bleo_ov_least, bleo_ov_least_sensitive_1se)

# PAAD
bleo_paad_idx <- bleomycin_test$TCGA_class == 'PAAD' #14

bleo_paad_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_paad_idx]
actual_bleo_paad_most         <- bleomycin_test$most_sensitive[bleo_paad_idx]
bleo_paad_most_min_auc        <- auc(actual_bleo_paad_most, bleo_paad_most_sensitive_min)

bleo_paad_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_paad_idx]
actual_bleo_paad_least <- bleomycin_test$least_sensitive[bleo_paad_idx]
bleo_paad_least_min_auc <- auc(actual_bleo_paad_least, bleo_paad_least_sensitive_min)

bleo_paad_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_paad_idx]
actual_bleo_paad_most <- bleomycin_test$most_sensitive[bleo_paad_idx]
bleo_paad_most_1se_auc <- auc(actual_bleo_paad_most, bleo_paad_most_sensitive_1se)

bleo_paad_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_paad_idx]
actual_bleo_paad_least <- bleomycin_test$least_sensitive[bleo_paad_idx]
bleo_paad_least_1se_auc <- auc(actual_bleo_paad_least, bleo_paad_least_sensitive_1se)

# PRAD
bleo_prad_idx <- bleomycin_test$TCGA_class == 'PRAD' #3

bleo_prad_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_prad_idx]
actual_bleo_prad_most         <- bleomycin_test$most_sensitive[bleo_prad_idx]
bleo_prad_most_min_auc        <- auc(actual_bleo_prad_most, bleo_prad_most_sensitive_min)

bleo_prad_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_prad_idx]
actual_bleo_prad_least <- bleomycin_test$least_sensitive[bleo_prad_idx]
bleo_prad_least_min_auc <- auc(actual_bleo_prad_least, bleo_prad_least_sensitive_min)

bleo_prad_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_prad_idx]
actual_bleo_prad_most <- bleomycin_test$most_sensitive[bleo_prad_idx]
bleo_prad_most_1se_auc <- auc(actual_bleo_prad_most, bleo_prad_most_sensitive_1se)

bleo_prad_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_prad_idx]
actual_bleo_prad_least <- bleomycin_test$least_sensitive[bleo_prad_idx]
bleo_prad_least_1se_auc <- auc(actual_bleo_prad_least, bleo_prad_least_sensitive_1se)

# SCLC
bleo_sclc_idx <- bleomycin_test$TCGA_class == 'SCLC' #30

bleo_sclc_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_sclc_idx]
actual_bleo_sclc_most         <- bleomycin_test$most_sensitive[bleo_sclc_idx]
bleo_sclc_most_min_auc        <- auc(actual_bleo_sclc_most, bleo_sclc_most_sensitive_min)

bleo_sclc_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_sclc_idx]
actual_bleo_sclc_least <- bleomycin_test$least_sensitive[bleo_sclc_idx]
bleo_sclc_least_min_auc <- auc(actual_bleo_sclc_least, bleo_sclc_least_sensitive_min)

bleo_sclc_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_sclc_idx]
actual_bleo_sclc_most <- bleomycin_test$most_sensitive[bleo_sclc_idx]
bleo_sclc_most_1se_auc <- auc(actual_bleo_sclc_most, bleo_sclc_most_sensitive_1se)

bleo_sclc_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_sclc_idx]
actual_bleo_sclc_least <- bleomycin_test$least_sensitive[bleo_sclc_idx]
bleo_sclc_least_1se_auc <- auc(actual_bleo_sclc_least, bleo_sclc_least_sensitive_1se)

# SKCM
bleo_skcm_idx <- bleomycin_test$TCGA_class == 'SKCM' #31

bleo_skcm_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_skcm_idx]
actual_bleo_skcm_most         <- bleomycin_test$most_sensitive[bleo_skcm_idx]
bleo_skcm_most_min_auc        <- auc(actual_bleo_skcm_most, bleo_skcm_most_sensitive_min)

bleo_skcm_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_skcm_idx]
actual_bleo_skcm_least <- bleomycin_test$least_sensitive[bleo_skcm_idx]
bleo_skcm_least_min_auc <- auc(actual_bleo_skcm_least, bleo_skcm_least_sensitive_min)

bleo_skcm_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_skcm_idx]
actual_bleo_skcm_most <- bleomycin_test$most_sensitive[bleo_skcm_idx]
bleo_skcm_most_1se_auc <- auc(actual_bleo_skcm_most, bleo_skcm_most_sensitive_1se)

bleo_skcm_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_skcm_idx]
actual_bleo_skcm_least <- bleomycin_test$least_sensitive[bleo_skcm_idx]
bleo_skcm_least_1se_auc <- auc(actual_bleo_skcm_least, bleo_skcm_least_sensitive_1se)

# STAD
bleo_stad_idx <- bleomycin_test$TCGA_class == 'STAD' #8

bleo_stad_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_stad_idx]
actual_bleo_stad_most         <- bleomycin_test$most_sensitive[bleo_stad_idx]
bleo_stad_most_min_auc        <- auc(actual_bleo_stad_most, bleo_stad_most_sensitive_min)

bleo_stad_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_stad_idx]
actual_bleo_stad_least <- bleomycin_test$least_sensitive[bleo_stad_idx]
bleo_stad_least_min_auc <- auc(actual_bleo_stad_least, bleo_stad_least_sensitive_min)

bleo_stad_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_stad_idx]
actual_bleo_stad_most <- bleomycin_test$most_sensitive[bleo_stad_idx]
bleo_stad_most_1se_auc <- auc(actual_bleo_stad_most, bleo_stad_most_sensitive_1se)

bleo_stad_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_stad_idx]
actual_bleo_stad_least <- bleomycin_test$least_sensitive[bleo_stad_idx]
bleo_stad_least_1se_auc <- auc(actual_bleo_stad_least, bleo_stad_least_sensitive_1se)

# THCA
bleo_thca_idx <- bleomycin_test$TCGA_class == 'THCA' #5

bleo_thca_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_thca_idx]
actual_bleo_thca_most         <- bleomycin_test$most_sensitive[bleo_thca_idx]
bleo_thca_most_min_auc        <- auc(actual_bleo_thca_most, bleo_thca_most_sensitive_min)

bleo_thca_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_thca_idx]
actual_bleo_thca_least <- bleomycin_test$least_sensitive[bleo_thca_idx]
bleo_thca_least_min_auc <- auc(actual_bleo_thca_least, bleo_thca_least_sensitive_min)

bleo_thca_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_thca_idx]
actual_bleo_thca_most <- bleomycin_test$most_sensitive[bleo_thca_idx]
bleo_thca_most_1se_auc <- auc(actual_bleo_thca_most, bleo_thca_most_sensitive_1se)

bleo_thca_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_thca_idx]
actual_bleo_thca_least <- bleomycin_test$least_sensitive[bleo_thca_idx]
bleo_thca_least_1se_auc <- auc(actual_bleo_thca_least, bleo_thca_least_sensitive_1se)

# UCEC
bleo_ucec_idx <- bleomycin_test$TCGA_class == 'UCEC' #4

bleo_ucec_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_ucec_idx]
actual_bleo_ucec_most         <- bleomycin_test$most_sensitive[bleo_ucec_idx]
bleo_ucec_most_min_auc        <- auc(actual_bleo_ucec_most, bleo_ucec_most_sensitive_min)

bleo_ucec_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_ucec_idx]
actual_bleo_ucec_least <- bleomycin_test$least_sensitive[bleo_ucec_idx]
bleo_ucec_least_min_auc <- auc(actual_bleo_ucec_least, bleo_ucec_least_sensitive_min)

bleo_ucec_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_ucec_idx]
actual_bleo_ucec_most <- bleomycin_test$most_sensitive[bleo_ucec_idx]
bleo_ucec_most_1se_auc <- auc(actual_bleo_ucec_most, bleo_ucec_most_sensitive_1se)

bleo_ucec_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_ucec_idx]
actual_bleo_ucec_least <- bleomycin_test$least_sensitive[bleo_ucec_idx]
bleo_ucec_least_1se_auc <- auc(actual_bleo_ucec_least, bleo_ucec_least_sensitive_1se)

# UNCLASSIFIED
bleo_unclassified_idx <- bleomycin_test$TCGA_class == 'UNCLASSIFIED' #58

bleo_unclassified_most_sensitive_min  <- new_bleo_most_sensitive_min[bleo_unclassified_idx]
actual_bleo_unclassified_most         <- bleomycin_test$most_sensitive[bleo_unclassified_idx]
bleo_unclassified_most_min_auc        <- auc(actual_bleo_unclassified_most, bleo_unclassified_most_sensitive_min)

bleo_unclassified_least_sensitive_min <- new_bleo_least_sensitive_min[bleo_unclassified_idx]
actual_bleo_unclassified_least <- bleomycin_test$least_sensitive[bleo_unclassified_idx]
bleo_unclassified_least_min_auc <- auc(actual_bleo_unclassified_least, bleo_unclassified_least_sensitive_min)

bleo_unclassified_most_sensitive_1se <- new_bleo_most_sensitive_1se[bleo_unclassified_idx]
actual_bleo_unclassified_most <- bleomycin_test$most_sensitive[bleo_unclassified_idx]
bleo_unclassified_most_1se_auc <- auc(actual_bleo_unclassified_most, bleo_unclassified_most_sensitive_1se)

bleo_unclassified_least_sensitive_1se <- new_bleo_least_sensitive_1se[bleo_unclassified_idx]
actual_bleo_unclassified_least <- bleomycin_test$least_sensitive[bleo_unclassified_idx]
bleo_unclassified_least_1se_auc <- auc(actual_bleo_unclassified_least, bleo_unclassified_least_sensitive_1se)


## CAMPTOTHECAN
table(camptothecan_test$TCGA_class)

# BLCA
campto_blca_idx <- camptothecan_test$TCGA_class == 'BLCA' #5

campto_blca_most_sensitive_min  <- new_campto_most_sensitive_min[campto_blca_idx]
actual_campto_blca_most         <- camptothecan_test$most_sensitive[campto_blca_idx]
campto_blca_most_min_auc        <- auc(actual_campto_blca_most, campto_blca_most_sensitive_min)

campto_blca_least_sensitive_min <- new_campto_least_sensitive_min[campto_blca_idx]
actual_campto_blca_least <- camptothecan_test$least_sensitive[campto_blca_idx]
campto_blca_least_min_auc <- auc(actual_campto_blca_least, campto_blca_least_sensitive_min)

campto_blca_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_blca_idx]
actual_campto_blca_most <- camptothecan_test$most_sensitive[campto_blca_idx]
campto_blca_most_1se_auc <- auc(actual_campto_blca_most, campto_blca_most_sensitive_1se)

campto_blca_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_blca_idx]
actual_campto_blca_least <- camptothecan_test$least_sensitive[campto_blca_idx]
campto_blca_least_1se_auc <- auc(actual_campto_blca_least, campto_blca_least_sensitive_1se)

# BRCA
campto_brca_idx <- camptothecan_test$TCGA_class == 'BRCA' #28

campto_brca_most_sensitive_min  <- new_campto_most_sensitive_min[campto_brca_idx]
actual_campto_brca_most         <- camptothecan_test$most_sensitive[campto_brca_idx]
campto_brca_most_min_auc        <- auc(actual_campto_brca_most, campto_brca_most_sensitive_min)

campto_brca_least_sensitive_min <- new_campto_least_sensitive_min[campto_brca_idx]
actual_campto_brca_least <- camptothecan_test$least_sensitive[campto_brca_idx]
campto_brca_least_min_auc <- auc(actual_campto_brca_least, campto_brca_least_sensitive_min)

campto_brca_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_brca_idx]
actual_campto_brca_most <- camptothecan_test$most_sensitive[campto_brca_idx]
campto_brca_most_1se_auc <- auc(actual_campto_brca_most, campto_brca_most_sensitive_1se)

campto_brca_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_brca_idx]
actual_campto_brca_least <- camptothecan_test$least_sensitive[campto_brca_idx]
campto_brca_least_1se_auc <- auc(actual_campto_brca_least, campto_brca_least_sensitive_1se)

# CESC
campto_cesc_idx <- camptothecan_test$TCGA_class == 'CESC' #4

campto_cesc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_cesc_idx]
actual_campto_cesc_most         <- camptothecan_test$most_sensitive[campto_cesc_idx]
campto_cesc_most_min_auc        <- auc(actual_campto_cesc_most, campto_cesc_most_sensitive_min)

campto_cesc_least_sensitive_min <- new_campto_least_sensitive_min[campto_cesc_idx]
actual_campto_cesc_least <- camptothecan_test$least_sensitive[campto_cesc_idx]
campto_cesc_least_min_auc <- auc(actual_campto_cesc_least, campto_cesc_least_sensitive_min)

campto_cesc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_cesc_idx]
actual_campto_cesc_most <- camptothecan_test$most_sensitive[campto_cesc_idx]
campto_cesc_most_1se_auc <- auc(actual_campto_cesc_most, campto_cesc_most_sensitive_1se)

campto_cesc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_cesc_idx]
actual_campto_cesc_least <- camptothecan_test$least_sensitive[campto_cesc_idx]
campto_cesc_least_1se_auc <- auc(actual_campto_cesc_least, campto_cesc_least_sensitive_1se)

# COREAD
campto_coread_idx <- camptothecan_test$TCGA_class == 'COREAD' #18

campto_coread_most_sensitive_min  <- new_campto_most_sensitive_min[campto_coread_idx]
actual_campto_coread_most         <- camptothecan_test$most_sensitive[campto_coread_idx]
campto_coread_most_min_auc        <- auc(actual_campto_coread_most, campto_coread_most_sensitive_min)

campto_coread_least_sensitive_min <- new_campto_least_sensitive_min[campto_coread_idx]
actual_campto_coread_least <- camptothecan_test$least_sensitive[campto_coread_idx]
campto_coread_least_min_auc <- auc(actual_campto_coread_least, campto_coread_least_sensitive_min)

campto_coread_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_coread_idx]
actual_campto_coread_most <- camptothecan_test$most_sensitive[campto_coread_idx]
campto_coread_most_1se_auc <- auc(actual_campto_coread_most, campto_coread_most_sensitive_1se)

campto_coread_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_coread_idx]
actual_campto_coread_least <- camptothecan_test$least_sensitive[campto_coread_idx]
campto_coread_least_1se_auc <- auc(actual_campto_coread_least, campto_coread_least_sensitive_1se)

# ESCA
campto_esca_idx <- camptothecan_test$TCGA_class == 'ESCA' #13

campto_esca_most_sensitive_min  <- new_campto_most_sensitive_min[campto_esca_idx]
actual_campto_esca_most         <- camptothecan_test$most_sensitive[campto_esca_idx]
campto_esca_most_min_auc        <- auc(actual_campto_esca_most, campto_esca_most_sensitive_min)

campto_esca_least_sensitive_min <- new_campto_least_sensitive_min[campto_esca_idx]
actual_campto_esca_least <- camptothecan_test$least_sensitive[campto_esca_idx]
campto_esca_least_min_auc <- auc(actual_campto_esca_least, campto_esca_least_sensitive_min)

campto_esca_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_esca_idx]
actual_campto_esca_most <- camptothecan_test$most_sensitive[campto_esca_idx]
campto_esca_most_1se_auc <- auc(actual_campto_esca_most, campto_esca_most_sensitive_1se)

campto_esca_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_esca_idx]
actual_campto_esca_least <- camptothecan_test$least_sensitive[campto_esca_idx]
campto_esca_least_1se_auc <- auc(actual_campto_esca_least, campto_esca_least_sensitive_1se)

# GBM
campto_gbm_idx <- camptothecan_test$TCGA_class == 'GBM' #17

campto_gbm_most_sensitive_min  <- new_campto_most_sensitive_min[campto_gbm_idx]
actual_campto_gbm_most         <- camptothecan_test$most_sensitive[campto_gbm_idx]
campto_gbm_most_min_auc        <- auc(actual_campto_gbm_most, campto_gbm_most_sensitive_min)

campto_gbm_least_sensitive_min <- new_campto_least_sensitive_min[campto_gbm_idx]
actual_campto_gbm_least <- camptothecan_test$least_sensitive[campto_gbm_idx]
campto_gbm_least_min_auc <- auc(actual_campto_gbm_least, campto_gbm_least_sensitive_min)

campto_gbm_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_gbm_idx]
actual_campto_gbm_most <- camptothecan_test$most_sensitive[campto_gbm_idx]
campto_gbm_most_1se_auc <- auc(actual_campto_gbm_most, campto_gbm_most_sensitive_1se)

campto_gbm_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_gbm_idx]
actual_campto_gbm_least <- camptothecan_test$least_sensitive[campto_gbm_idx]
campto_gbm_least_1se_auc <- auc(actual_campto_gbm_least, campto_gbm_least_sensitive_1se)

# HNSC
campto_hnsc_idx <- camptothecan_test$TCGA_class == 'HNSC' #19

campto_hnsc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_hnsc_idx]
actual_campto_hnsc_most         <- camptothecan_test$most_sensitive[campto_hnsc_idx]
campto_hnsc_most_min_auc        <- auc(actual_campto_hnsc_most, campto_hnsc_most_sensitive_min)

campto_hnsc_least_sensitive_min <- new_campto_least_sensitive_min[campto_hnsc_idx]
actual_campto_hnsc_least <- camptothecan_test$least_sensitive[campto_hnsc_idx]
campto_hnsc_least_min_auc <- auc(actual_campto_hnsc_least, campto_hnsc_least_sensitive_min)

campto_hnsc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_hnsc_idx]
actual_campto_hnsc_most <- camptothecan_test$most_sensitive[campto_hnsc_idx]
campto_hnsc_most_1se_auc <- auc(actual_campto_hnsc_most, campto_hnsc_most_sensitive_1se)

campto_hnsc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_hnsc_idx]
actual_campto_hnsc_least <- camptothecan_test$least_sensitive[campto_hnsc_idx]
campto_hnsc_least_1se_auc <- auc(actual_campto_hnsc_least, campto_hnsc_least_sensitive_1se)

# KIRC
campto_kirc_idx <- camptothecan_test$TCGA_class == 'KIRC' #15

campto_kirc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_kirc_idx]
actual_campto_kirc_most         <- camptothecan_test$most_sensitive[campto_kirc_idx]
campto_kirc_most_min_auc        <- auc(actual_campto_kirc_most, campto_kirc_most_sensitive_min)

campto_kirc_least_sensitive_min <- new_campto_least_sensitive_min[campto_kirc_idx]
actual_campto_kirc_least <- camptothecan_test$least_sensitive[campto_kirc_idx]
campto_kirc_least_min_auc <- auc(actual_campto_kirc_least, campto_kirc_least_sensitive_min)

campto_kirc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_kirc_idx]
actual_campto_kirc_most <- camptothecan_test$most_sensitive[campto_kirc_idx]
campto_kirc_most_1se_auc <- auc(actual_campto_kirc_most, campto_kirc_most_sensitive_1se)

campto_kirc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_kirc_idx]
actual_campto_kirc_least <- camptothecan_test$least_sensitive[campto_kirc_idx]
campto_kirc_least_1se_auc <- auc(actual_campto_kirc_least, campto_kirc_least_sensitive_1se)

# LGG
campto_lgg_idx <- camptothecan_test$TCGA_class == 'LGG' #8

campto_lgg_most_sensitive_min  <- new_campto_most_sensitive_min[campto_lgg_idx]
actual_campto_lgg_most         <- camptothecan_test$most_sensitive[campto_lgg_idx]
campto_lgg_most_min_auc        <- auc(actual_campto_lgg_most, campto_lgg_most_sensitive_min)

campto_lgg_least_sensitive_min <- new_campto_least_sensitive_min[campto_lgg_idx]
actual_campto_lgg_least <- camptothecan_test$least_sensitive[campto_lgg_idx]
campto_lgg_least_min_auc <- auc(actual_campto_lgg_least, campto_lgg_least_sensitive_min)

campto_lgg_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_lgg_idx]
actual_campto_lgg_most <- camptothecan_test$most_sensitive[campto_lgg_idx]
campto_lgg_most_1se_auc <- auc(actual_campto_lgg_most, campto_lgg_most_sensitive_1se)

campto_lgg_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_lgg_idx]
actual_campto_lgg_least <- camptothecan_test$least_sensitive[campto_lgg_idx]
campto_lgg_least_1se_auc <- auc(actual_campto_lgg_least, campto_lgg_least_sensitive_1se)

# LIHC
campto_lihc_idx <- camptothecan_test$TCGA_class == 'LIHC' #7

campto_lihc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_lihc_idx]
actual_campto_lihc_most         <- camptothecan_test$most_sensitive[campto_lihc_idx]
campto_lihc_most_min_auc        <- auc(actual_campto_lihc_most, campto_lihc_most_sensitive_min)

campto_lihc_least_sensitive_min <- new_campto_least_sensitive_min[campto_lihc_idx]
actual_campto_lihc_least <- camptothecan_test$least_sensitive[campto_lihc_idx]
campto_lihc_least_min_auc <- auc(actual_campto_lihc_least, campto_lihc_least_sensitive_min)

campto_lihc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_lihc_idx]
actual_campto_lihc_most <- camptothecan_test$most_sensitive[campto_lihc_idx]
campto_lihc_most_1se_auc <- auc(actual_campto_lihc_most, campto_lihc_most_sensitive_1se)

campto_lihc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_lihc_idx]
actual_campto_lihc_least <- camptothecan_test$least_sensitive[campto_lihc_idx]
campto_lihc_least_1se_auc <- auc(actual_campto_lihc_least, campto_lihc_least_sensitive_1se)

# LUAD
campto_luad_idx <- camptothecan_test$TCGA_class == 'LUAD' #26

campto_luad_most_sensitive_min  <- new_campto_most_sensitive_min[campto_luad_idx]
actual_campto_luad_most         <- camptothecan_test$most_sensitive[campto_luad_idx]
campto_luad_most_min_auc        <- auc(actual_campto_luad_most, campto_luad_most_sensitive_min)

campto_luad_least_sensitive_min <- new_campto_least_sensitive_min[campto_luad_idx]
actual_campto_luad_least <- camptothecan_test$least_sensitive[campto_luad_idx]
campto_luad_least_min_auc <- auc(actual_campto_luad_least, campto_luad_least_sensitive_min)

campto_luad_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_luad_idx]
actual_campto_luad_most <- camptothecan_test$most_sensitive[campto_luad_idx]
campto_luad_most_1se_auc <- auc(actual_campto_luad_most, campto_luad_most_sensitive_1se)

campto_luad_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_luad_idx]
actual_campto_luad_least <- camptothecan_test$least_sensitive[campto_luad_idx]
campto_luad_least_1se_auc <- auc(actual_campto_luad_least, campto_luad_least_sensitive_1se)

# LUSC
campto_lusc_idx <- camptothecan_test$TCGA_class == 'LUSC' #5

campto_lusc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_lusc_idx]
actual_campto_lusc_most         <- camptothecan_test$most_sensitive[campto_lusc_idx]
campto_lusc_most_min_auc        <- auc(actual_campto_lusc_most, campto_lusc_most_sensitive_min)

campto_lusc_least_sensitive_min <- new_campto_least_sensitive_min[campto_lusc_idx]
actual_campto_lusc_least <- camptothecan_test$least_sensitive[campto_lusc_idx]
campto_lusc_least_min_auc <- auc(actual_campto_lusc_least, campto_lusc_least_sensitive_min)

campto_lusc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_lusc_idx]
actual_campto_lusc_most <- camptothecan_test$most_sensitive[campto_lusc_idx]
campto_lusc_most_1se_auc <- auc(actual_campto_lusc_most, campto_lusc_most_sensitive_1se)

campto_lusc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_lusc_idx]
actual_campto_lusc_least <- camptothecan_test$least_sensitive[campto_lusc_idx]
campto_lusc_least_1se_auc <- auc(actual_campto_lusc_least, campto_lusc_least_sensitive_1se)

# MB
campto_mb_idx <- camptothecan_test$TCGA_class == 'MB' #3

campto_mb_most_sensitive_min  <- new_campto_most_sensitive_min[campto_mb_idx]
actual_campto_mb_most         <- camptothecan_test$most_sensitive[campto_mb_idx]
campto_mb_most_min_auc        <- auc(actual_campto_mb_most, campto_mb_most_sensitive_min)

campto_mb_least_sensitive_min <- new_campto_least_sensitive_min[campto_mb_idx]
actual_campto_mb_least <- camptothecan_test$least_sensitive[campto_mb_idx]
campto_mb_least_min_auc <- auc(actual_campto_mb_least, campto_mb_least_sensitive_min)

campto_mb_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_mb_idx]
actual_campto_mb_most <- camptothecan_test$most_sensitive[campto_mb_idx]
campto_mb_most_1se_auc <- auc(actual_campto_mb_most, campto_mb_most_sensitive_1se)

campto_mb_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_mb_idx]
actual_campto_mb_least <- camptothecan_test$least_sensitive[campto_mb_idx]
campto_mb_least_1se_auc <- auc(actual_campto_mb_least, campto_mb_least_sensitive_1se)

# MESO
campto_meso_idx <- camptothecan_test$TCGA_class == 'MESO' #5

campto_meso_most_sensitive_min  <- new_campto_most_sensitive_min[campto_meso_idx]
actual_campto_meso_most         <- camptothecan_test$most_sensitive[campto_meso_idx]
campto_meso_most_min_auc        <- auc(actual_campto_meso_most, campto_meso_most_sensitive_min)

campto_meso_least_sensitive_min <- new_campto_least_sensitive_min[campto_meso_idx]
actual_campto_meso_least <- camptothecan_test$least_sensitive[campto_meso_idx]
campto_meso_least_min_auc <- auc(actual_campto_meso_least, campto_meso_least_sensitive_min)

campto_meso_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_meso_idx]
actual_campto_meso_most <- camptothecan_test$most_sensitive[campto_meso_idx]
campto_meso_most_1se_auc <- auc(actual_campto_meso_most, campto_meso_most_sensitive_1se)

campto_meso_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_meso_idx]
actual_campto_meso_least <- camptothecan_test$least_sensitive[campto_meso_idx]
campto_meso_least_1se_auc <- auc(actual_campto_meso_least, campto_meso_least_sensitive_1se)

# NB
campto_nb_idx <- camptothecan_test$TCGA_class == 'NB' #13

campto_nb_most_sensitive_min  <- new_campto_most_sensitive_min[campto_nb_idx]
actual_campto_nb_most         <- camptothecan_test$most_sensitive[campto_nb_idx]
campto_nb_most_min_auc        <- auc(actual_campto_nb_most, campto_nb_most_sensitive_min)

campto_nb_least_sensitive_min <- new_campto_least_sensitive_min[campto_nb_idx]
actual_campto_nb_least <- camptothecan_test$least_sensitive[campto_nb_idx]
campto_nb_least_min_auc <- auc(actual_campto_nb_least, campto_nb_least_sensitive_min)

campto_nb_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_nb_idx]
actual_campto_nb_most <- camptothecan_test$most_sensitive[campto_nb_idx]
campto_nb_most_1se_auc <- auc(actual_campto_nb_most, campto_nb_most_sensitive_1se)

campto_nb_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_nb_idx]
actual_campto_nb_least <- camptothecan_test$least_sensitive[campto_nb_idx]
campto_nb_least_1se_auc <- auc(actual_campto_nb_least, campto_nb_least_sensitive_1se)

# OV
campto_ov_idx <- camptothecan_test$TCGA_class == 'OV' #17

campto_ov_most_sensitive_min  <- new_campto_most_sensitive_min[campto_ov_idx]
actual_campto_ov_most         <- camptothecan_test$most_sensitive[campto_ov_idx]
campto_ov_most_min_auc        <- auc(actual_campto_ov_most, campto_ov_most_sensitive_min)

campto_ov_least_sensitive_min <- new_campto_least_sensitive_min[campto_ov_idx]
actual_campto_ov_least <- camptothecan_test$least_sensitive[campto_ov_idx]
campto_ov_least_min_auc <- auc(actual_campto_ov_least, campto_ov_least_sensitive_min)

campto_ov_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_ov_idx]
actual_campto_ov_most <- camptothecan_test$most_sensitive[campto_ov_idx]
campto_ov_most_1se_auc <- auc(actual_campto_ov_most, campto_ov_most_sensitive_1se)

campto_ov_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_ov_idx]
actual_campto_ov_least <- camptothecan_test$least_sensitive[campto_ov_idx]
campto_ov_least_1se_auc <- auc(actual_campto_ov_least, campto_ov_least_sensitive_1se)

# PAAD
campto_paad_idx <- camptothecan_test$TCGA_class == 'PAAD' #13

campto_paad_most_sensitive_min  <- new_campto_most_sensitive_min[campto_paad_idx]
actual_campto_paad_most         <- camptothecan_test$most_sensitive[campto_paad_idx]
campto_paad_most_min_auc        <- auc(actual_campto_paad_most, campto_paad_most_sensitive_min)

campto_paad_least_sensitive_min <- new_campto_least_sensitive_min[campto_paad_idx]
actual_campto_paad_least <- camptothecan_test$least_sensitive[campto_paad_idx]
campto_paad_least_min_auc <- auc(actual_campto_paad_least, campto_paad_least_sensitive_min)

campto_paad_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_paad_idx]
actual_campto_paad_most <- camptothecan_test$most_sensitive[campto_paad_idx]
campto_paad_most_1se_auc <- auc(actual_campto_paad_most, campto_paad_most_sensitive_1se)

campto_paad_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_paad_idx]
actual_campto_paad_least <- camptothecan_test$least_sensitive[campto_paad_idx]
campto_paad_least_1se_auc <- auc(actual_campto_paad_least, campto_paad_least_sensitive_1se)

# SCLC
campto_sclc_idx <- camptothecan_test$TCGA_class == 'SCLC' #31

campto_sclc_most_sensitive_min  <- new_campto_most_sensitive_min[campto_sclc_idx]
actual_campto_sclc_most         <- camptothecan_test$most_sensitive[campto_sclc_idx]
campto_sclc_most_min_auc        <- auc(actual_campto_sclc_most, campto_sclc_most_sensitive_min)

campto_sclc_least_sensitive_min <- new_campto_least_sensitive_min[campto_sclc_idx]
actual_campto_sclc_least <- camptothecan_test$least_sensitive[campto_sclc_idx]
campto_sclc_least_min_auc <- auc(actual_campto_sclc_least, campto_sclc_least_sensitive_min)

campto_sclc_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_sclc_idx]
actual_campto_sclc_most <- camptothecan_test$most_sensitive[campto_sclc_idx]
campto_sclc_most_1se_auc <- auc(actual_campto_sclc_most, campto_sclc_most_sensitive_1se)

campto_sclc_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_sclc_idx]
actual_campto_sclc_least <- camptothecan_test$least_sensitive[campto_sclc_idx]
campto_sclc_least_1se_auc <- auc(actual_campto_sclc_least, campto_sclc_least_sensitive_1se)

# SKCM
campto_skcm_idx <- camptothecan_test$TCGA_class == 'SKCM' #26

campto_skcm_most_sensitive_min  <- new_campto_most_sensitive_min[campto_skcm_idx]
actual_campto_skcm_most         <- camptothecan_test$most_sensitive[campto_skcm_idx]
campto_skcm_most_min_auc        <- auc(actual_campto_skcm_most, campto_skcm_most_sensitive_min)

campto_skcm_least_sensitive_min <- new_campto_least_sensitive_min[campto_skcm_idx]
actual_campto_skcm_least <- camptothecan_test$least_sensitive[campto_skcm_idx]
campto_skcm_least_min_auc <- auc(actual_campto_skcm_least, campto_skcm_least_sensitive_min)

campto_skcm_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_skcm_idx]
actual_campto_skcm_most <- camptothecan_test$most_sensitive[campto_skcm_idx]
campto_skcm_most_1se_auc <- auc(actual_campto_skcm_most, campto_skcm_most_sensitive_1se)

campto_skcm_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_skcm_idx]
actual_campto_skcm_least <- camptothecan_test$least_sensitive[campto_skcm_idx]
campto_skcm_least_1se_auc <- auc(actual_campto_skcm_least, campto_skcm_least_sensitive_1se)

# STAD
campto_stad_idx <- camptothecan_test$TCGA_class == 'STAD' #8

campto_stad_most_sensitive_min  <- new_campto_most_sensitive_min[campto_stad_idx]
actual_campto_stad_most         <- camptothecan_test$most_sensitive[campto_stad_idx]
campto_stad_most_min_auc        <- auc(actual_campto_stad_most, campto_stad_most_sensitive_min)

campto_stad_least_sensitive_min <- new_campto_least_sensitive_min[campto_stad_idx]
actual_campto_stad_least <- camptothecan_test$least_sensitive[campto_stad_idx]
campto_stad_least_min_auc <- auc(actual_campto_stad_least, campto_stad_least_sensitive_min)

campto_stad_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_stad_idx]
actual_campto_stad_most <- camptothecan_test$most_sensitive[campto_stad_idx]
campto_stad_most_1se_auc <- auc(actual_campto_stad_most, campto_stad_most_sensitive_1se)

campto_stad_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_stad_idx]
actual_campto_stad_least <- camptothecan_test$least_sensitive[campto_stad_idx]
campto_stad_least_1se_auc <- auc(actual_campto_stad_least, campto_stad_least_sensitive_1se)

# THCA
campto_thca_idx <- camptothecan_test$TCGA_class == 'THCA' #4

campto_thca_most_sensitive_min  <- new_campto_most_sensitive_min[campto_thca_idx]
actual_campto_thca_most         <- camptothecan_test$most_sensitive[campto_thca_idx]
campto_thca_most_min_auc        <- auc(actual_campto_thca_most, campto_thca_most_sensitive_min)

campto_thca_least_sensitive_min <- new_campto_least_sensitive_min[campto_thca_idx]
actual_campto_thca_least <- camptothecan_test$least_sensitive[campto_thca_idx]
campto_thca_least_min_auc <- auc(actual_campto_thca_least, campto_thca_least_sensitive_min)

campto_thca_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_thca_idx]
actual_campto_thca_most <- camptothecan_test$most_sensitive[campto_thca_idx]
campto_thca_most_1se_auc <- auc(actual_campto_thca_most, campto_thca_most_sensitive_1se)

campto_thca_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_thca_idx]
actual_campto_thca_least <- camptothecan_test$least_sensitive[campto_thca_idx]
campto_thca_least_1se_auc <- auc(actual_campto_thca_least, campto_thca_least_sensitive_1se)

# UCEC
campto_ucec_idx <- camptothecan_test$TCGA_class == 'UCEC' #4

campto_ucec_most_sensitive_min  <- new_campto_most_sensitive_min[campto_ucec_idx]
actual_campto_ucec_most         <- camptothecan_test$most_sensitive[campto_ucec_idx]
campto_ucec_most_min_auc        <- auc(actual_campto_ucec_most, campto_ucec_most_sensitive_min)

campto_ucec_least_sensitive_min <- new_campto_least_sensitive_min[campto_ucec_idx]
actual_campto_ucec_least <- camptothecan_test$least_sensitive[campto_ucec_idx]
campto_ucec_least_min_auc <- auc(actual_campto_ucec_least, campto_ucec_least_sensitive_min)

campto_ucec_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_ucec_idx]
actual_campto_ucec_most <- camptothecan_test$most_sensitive[campto_ucec_idx]
campto_ucec_most_1se_auc <- auc(actual_campto_ucec_most, campto_ucec_most_sensitive_1se)

campto_ucec_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_ucec_idx]
actual_campto_ucec_least <- camptothecan_test$least_sensitive[campto_ucec_idx]
campto_ucec_least_1se_auc <- auc(actual_campto_ucec_least, campto_ucec_least_sensitive_1se)

# UNCLASSIFIED
campto_unclassified_idx <- camptothecan_test$TCGA_class == 'UNCLASSIFIED' #52

campto_unclassified_most_sensitive_min  <- new_campto_most_sensitive_min[campto_unclassified_idx]
actual_campto_unclassified_most         <- camptothecan_test$most_sensitive[campto_unclassified_idx]
campto_unclassified_most_min_auc        <- auc(actual_campto_unclassified_most, campto_unclassified_most_sensitive_min)

campto_unclassified_least_sensitive_min <- new_campto_least_sensitive_min[campto_unclassified_idx]
actual_campto_unclassified_least <- camptothecan_test$least_sensitive[campto_unclassified_idx]
campto_unclassified_least_min_auc <- auc(actual_campto_unclassified_least, campto_unclassified_least_sensitive_min)

campto_unclassified_most_sensitive_1se <- new_campto_most_sensitive_1se[campto_unclassified_idx]
actual_campto_unclassified_most <- camptothecan_test$most_sensitive[campto_unclassified_idx]
campto_unclassified_most_1se_auc <- auc(actual_campto_unclassified_most, campto_unclassified_most_sensitive_1se)

campto_unclassified_least_sensitive_1se <- new_campto_least_sensitive_1se[campto_unclassified_idx]
actual_campto_unclassified_least <- camptothecan_test$least_sensitive[campto_unclassified_idx]
campto_unclassified_least_1se_auc <- auc(actual_campto_unclassified_least, campto_unclassified_least_sensitive_1se)

## CISPLATIN
table(cisplatin_test$TCGA_class)

# BLCA 5
cisplatin_blca_idx <- cisplatin_test$TCGA_class == 'BLCA' #5

cisplatin_blca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_blca_idx]
actual_cisplatin_blca_most         <- cisplatin_test$most_sensitive[cisplatin_blca_idx]
cisplatin_blca_most_min_auc        <- auc(actual_cisplatin_blca_most, cisplatin_blca_most_sensitive_min)

cisplatin_blca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_blca_idx]
actual_cisplatin_blca_least <- cisplatin_test$least_sensitive[cisplatin_blca_idx]
cisplatin_blca_least_min_auc <- auc(actual_cisplatin_blca_least, cisplatin_blca_least_sensitive_min)

cisplatin_blca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_blca_idx]
actual_cisplatin_blca_most <- cisplatin_test$most_sensitive[cisplatin_blca_idx]
cisplatin_blca_most_1se_auc <- auc(actual_cisplatin_blca_most, cisplatin_blca_most_sensitive_1se)

cisplatin_blca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_blca_idx]
actual_cisplatin_blca_least <- cisplatin_test$least_sensitive[cisplatin_blca_idx]
cisplatin_blca_least_1se_auc <- auc(actual_cisplatin_blca_least, cisplatin_blca_least_sensitive_1se)

# BRCA
cisplatin_brca_idx <- cisplatin_test$TCGA_class == 'BRCA' #28

cisplatin_brca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_brca_idx]
actual_cisplatin_brca_most         <- cisplatin_test$most_sensitive[cisplatin_brca_idx]
cisplatin_brca_most_min_auc        <- auc(actual_cisplatin_brca_most, cisplatin_brca_most_sensitive_min)

cisplatin_brca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_brca_idx]
actual_cisplatin_brca_least <- cisplatin_test$least_sensitive[cisplatin_brca_idx]
cisplatin_brca_least_min_auc <- auc(actual_cisplatin_brca_least, cisplatin_brca_least_sensitive_min)

cisplatin_brca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_brca_idx]
actual_cisplatin_brca_most <- cisplatin_test$most_sensitive[cisplatin_brca_idx]
cisplatin_brca_most_1se_auc <- auc(actual_cisplatin_brca_most, cisplatin_brca_most_sensitive_1se)

cisplatin_brca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_brca_idx]
actual_cisplatin_brca_least <- cisplatin_test$least_sensitive[cisplatin_brca_idx]
cisplatin_brca_least_1se_auc <- auc(actual_cisplatin_brca_least, cisplatin_brca_least_sensitive_1se)

# CESC
cisplatin_cesc_idx <- cisplatin_test$TCGA_class == 'CESC' #4

cisplatin_cesc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_cesc_idx]
actual_cisplatin_cesc_most         <- cisplatin_test$most_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_most_min_auc        <- auc(actual_cisplatin_cesc_most, cisplatin_cesc_most_sensitive_min)

cisplatin_cesc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_cesc_idx]
actual_cisplatin_cesc_least <- cisplatin_test$least_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_least_min_auc <- auc(actual_cisplatin_cesc_least, cisplatin_cesc_least_sensitive_min)

cisplatin_cesc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_cesc_idx]
actual_cisplatin_cesc_most <- cisplatin_test$most_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_most_1se_auc <- auc(actual_cisplatin_cesc_most, cisplatin_cesc_most_sensitive_1se)

cisplatin_cesc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_cesc_idx]
actual_cisplatin_cesc_least <- cisplatin_test$least_sensitive[cisplatin_cesc_idx]
cisplatin_cesc_least_1se_auc <- auc(actual_cisplatin_cesc_least, cisplatin_cesc_least_sensitive_1se)

# COREAD
cisplatin_coread_idx <- cisplatin_test$TCGA_class == 'COREAD' #18

cisplatin_coread_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_coread_idx]
actual_cisplatin_coread_most         <- cisplatin_test$most_sensitive[cisplatin_coread_idx]
cisplatin_coread_most_min_auc        <- auc(actual_cisplatin_coread_most, cisplatin_coread_most_sensitive_min)

cisplatin_coread_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_coread_idx]
actual_cisplatin_coread_least <- cisplatin_test$least_sensitive[cisplatin_coread_idx]
cisplatin_coread_least_min_auc <- auc(actual_cisplatin_coread_least, cisplatin_coread_least_sensitive_min)

cisplatin_coread_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_coread_idx]
actual_cisplatin_coread_most <- cisplatin_test$most_sensitive[cisplatin_coread_idx]
cisplatin_coread_most_1se_auc <- auc(actual_cisplatin_coread_most, cisplatin_coread_most_sensitive_1se)

cisplatin_coread_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_coread_idx]
actual_cisplatin_coread_least <- cisplatin_test$least_sensitive[cisplatin_coread_idx]
cisplatin_coread_least_1se_auc <- auc(actual_cisplatin_coread_least, cisplatin_coread_least_sensitive_1se)

# ESCA
cisplatin_esca_idx <- cisplatin_test$TCGA_class == 'ESCA' #13

cisplatin_esca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_esca_idx]
actual_cisplatin_esca_most         <- cisplatin_test$most_sensitive[cisplatin_esca_idx]
cisplatin_esca_most_min_auc        <- auc(actual_cisplatin_esca_most, cisplatin_esca_most_sensitive_min)

cisplatin_esca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_esca_idx]
actual_cisplatin_esca_least <- cisplatin_test$least_sensitive[cisplatin_esca_idx]
cisplatin_esca_least_min_auc <- auc(actual_cisplatin_esca_least, cisplatin_esca_least_sensitive_min)

cisplatin_esca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_esca_idx]
actual_cisplatin_esca_most <- cisplatin_test$most_sensitive[cisplatin_esca_idx]
cisplatin_esca_most_1se_auc <- auc(actual_cisplatin_esca_most, cisplatin_esca_most_sensitive_1se)

cisplatin_esca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_esca_idx]
actual_cisplatin_esca_least <- cisplatin_test$least_sensitive[cisplatin_esca_idx]
cisplatin_esca_least_1se_auc <- auc(actual_cisplatin_esca_least, cisplatin_esca_least_sensitive_1se)

# GBM
cisplatin_gbm_idx <- cisplatin_test$TCGA_class == 'GBM' #17

cisplatin_gbm_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_gbm_idx]
actual_cisplatin_gbm_most         <- cisplatin_test$most_sensitive[cisplatin_gbm_idx]
cisplatin_gbm_most_min_auc        <- auc(actual_cisplatin_gbm_most, cisplatin_gbm_most_sensitive_min)

cisplatin_gbm_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_gbm_idx]
actual_cisplatin_gbm_least <- cisplatin_test$least_sensitive[cisplatin_gbm_idx]
cisplatin_gbm_least_min_auc <- auc(actual_cisplatin_gbm_least, cisplatin_gbm_least_sensitive_min)

cisplatin_gbm_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_gbm_idx]
actual_cisplatin_gbm_most <- cisplatin_test$most_sensitive[cisplatin_gbm_idx]
cisplatin_gbm_most_1se_auc <- auc(actual_cisplatin_gbm_most, cisplatin_gbm_most_sensitive_1se)

cisplatin_gbm_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_gbm_idx]
actual_cisplatin_gbm_least <- cisplatin_test$least_sensitive[cisplatin_gbm_idx]
cisplatin_gbm_least_1se_auc <- auc(actual_cisplatin_gbm_least, cisplatin_gbm_least_sensitive_1se)

# HNSC
cisplatin_hnsc_idx <- cisplatin_test$TCGA_class == 'HNSC' #19

cisplatin_hnsc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_most         <- cisplatin_test$most_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_most_min_auc        <- auc(actual_cisplatin_hnsc_most, cisplatin_hnsc_most_sensitive_min)

cisplatin_hnsc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_least <- cisplatin_test$least_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_least_min_auc <- auc(actual_cisplatin_hnsc_least, cisplatin_hnsc_least_sensitive_min)

cisplatin_hnsc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_most <- cisplatin_test$most_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_most_1se_auc <- auc(actual_cisplatin_hnsc_most, cisplatin_hnsc_most_sensitive_1se)

cisplatin_hnsc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_hnsc_idx]
actual_cisplatin_hnsc_least <- cisplatin_test$least_sensitive[cisplatin_hnsc_idx]
cisplatin_hnsc_least_1se_auc <- auc(actual_cisplatin_hnsc_least, cisplatin_hnsc_least_sensitive_1se)

# KIRC
cisplatin_kirc_idx <- cisplatin_test$TCGA_class == 'KIRC' #15

cisplatin_kirc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_kirc_idx]
actual_cisplatin_kirc_most         <- cisplatin_test$most_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_most_min_auc        <- auc(actual_cisplatin_kirc_most, cisplatin_kirc_most_sensitive_min)

cisplatin_kirc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_kirc_idx]
actual_cisplatin_kirc_least <- cisplatin_test$least_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_least_min_auc <- auc(actual_cisplatin_kirc_least, cisplatin_kirc_least_sensitive_min)

cisplatin_kirc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_kirc_idx]
actual_cisplatin_kirc_most <- cisplatin_test$most_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_most_1se_auc <- auc(actual_cisplatin_kirc_most, cisplatin_kirc_most_sensitive_1se)

cisplatin_kirc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_kirc_idx]
actual_cisplatin_kirc_least <- cisplatin_test$least_sensitive[cisplatin_kirc_idx]
cisplatin_kirc_least_1se_auc <- auc(actual_cisplatin_kirc_least, cisplatin_kirc_least_sensitive_1se)

# LGG
cisplatin_lgg_idx <- cisplatin_test$TCGA_class == 'LGG' #8

cisplatin_lgg_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_lgg_idx]
actual_cisplatin_lgg_most         <- cisplatin_test$most_sensitive[cisplatin_lgg_idx]
cisplatin_lgg_most_min_auc        <- auc(actual_cisplatin_lgg_most, cisplatin_lgg_most_sensitive_min)

cisplatin_lgg_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_lgg_idx]
actual_cisplatin_lgg_least <- cisplatin_test$least_sensitive[cisplatin_lgg_idx]
cisplatin_lgg_least_min_auc <- auc(actual_cisplatin_lgg_least, cisplatin_lgg_least_sensitive_min)

cisplatin_lgg_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_lgg_idx]
actual_cisplatin_lgg_most <- cisplatin_test$most_sensitive[cisplatin_lgg_idx]
cisplatin_lgg_most_1se_auc <- auc(actual_cisplatin_lgg_most, cisplatin_lgg_most_sensitive_1se)

cisplatin_lgg_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_lgg_idx]
actual_cisplatin_lgg_least <- cisplatin_test$least_sensitive[cisplatin_lgg_idx]
cisplatin_lgg_least_1se_auc <- auc(actual_cisplatin_lgg_least, cisplatin_lgg_least_sensitive_1se)

# LIHC
cisplatin_lihc_idx <- cisplatin_test$TCGA_class == 'LIHC' #7

cisplatin_lihc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_lihc_idx]
actual_cisplatin_lihc_most         <- cisplatin_test$most_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_most_min_auc        <- auc(actual_cisplatin_lihc_most, cisplatin_lihc_most_sensitive_min)

cisplatin_lihc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_lihc_idx]
actual_cisplatin_lihc_least <- cisplatin_test$least_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_least_min_auc <- auc(actual_cisplatin_lihc_least, cisplatin_lihc_least_sensitive_min)

cisplatin_lihc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_lihc_idx]
actual_cisplatin_lihc_most <- cisplatin_test$most_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_most_1se_auc <- auc(actual_cisplatin_lihc_most, cisplatin_lihc_most_sensitive_1se)

cisplatin_lihc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_lihc_idx]
actual_cisplatin_lihc_least <- cisplatin_test$least_sensitive[cisplatin_lihc_idx]
cisplatin_lihc_least_1se_auc <- auc(actual_cisplatin_lihc_least, cisplatin_lihc_least_sensitive_1se)

# LUAD
cisplatin_luad_idx <- cisplatin_test$TCGA_class == 'LUAD' #26

cisplatin_luad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_luad_idx]
actual_cisplatin_luad_most         <- cisplatin_test$most_sensitive[cisplatin_luad_idx]
cisplatin_luad_most_min_auc        <- auc(actual_cisplatin_luad_most, cisplatin_luad_most_sensitive_min)

cisplatin_luad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_luad_idx]
actual_cisplatin_luad_least <- cisplatin_test$least_sensitive[cisplatin_luad_idx]
cisplatin_luad_least_min_auc <- auc(actual_cisplatin_luad_least, cisplatin_luad_least_sensitive_min)

cisplatin_luad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_luad_idx]
actual_cisplatin_luad_most <- cisplatin_test$most_sensitive[cisplatin_luad_idx]
cisplatin_luad_most_1se_auc <- auc(actual_cisplatin_luad_most, cisplatin_luad_most_sensitive_1se)

cisplatin_luad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_luad_idx]
actual_cisplatin_luad_least <- cisplatin_test$least_sensitive[cisplatin_luad_idx]
cisplatin_luad_least_1se_auc <- auc(actual_cisplatin_luad_least, cisplatin_luad_least_sensitive_1se)

# LUSC
cisplatin_lusc_idx <- cisplatin_test$TCGA_class == 'LUSC' #5

cisplatin_lusc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_lusc_idx]
actual_cisplatin_lusc_most         <- cisplatin_test$most_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_most_min_auc        <- auc(actual_cisplatin_lusc_most, cisplatin_lusc_most_sensitive_min)

cisplatin_lusc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_lusc_idx]
actual_cisplatin_lusc_least <- cisplatin_test$least_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_least_min_auc <- auc(actual_cisplatin_lusc_least, cisplatin_lusc_least_sensitive_min)

cisplatin_lusc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_lusc_idx]
actual_cisplatin_lusc_most <- cisplatin_test$most_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_most_1se_auc <- auc(actual_cisplatin_lusc_most, cisplatin_lusc_most_sensitive_1se)

cisplatin_lusc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_lusc_idx]
actual_cisplatin_lusc_least <- cisplatin_test$least_sensitive[cisplatin_lusc_idx]
cisplatin_lusc_least_1se_auc <- auc(actual_cisplatin_lusc_least, cisplatin_lusc_least_sensitive_1se)

# MB
cisplatin_mb_idx <- cisplatin_test$TCGA_class == 'MB' #3

cisplatin_mb_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_mb_idx]
actual_cisplatin_mb_most         <- cisplatin_test$most_sensitive[cisplatin_mb_idx]
cisplatin_mb_most_min_auc        <- auc(actual_cisplatin_mb_most, cisplatin_mb_most_sensitive_min)

cisplatin_mb_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_mb_idx]
actual_cisplatin_mb_least <- cisplatin_test$least_sensitive[cisplatin_mb_idx]
cisplatin_mb_least_min_auc <- auc(actual_cisplatin_mb_least, cisplatin_mb_least_sensitive_min)

cisplatin_mb_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_mb_idx]
actual_cisplatin_mb_most <- cisplatin_test$most_sensitive[cisplatin_mb_idx]
cisplatin_mb_most_1se_auc <- auc(actual_cisplatin_mb_most, cisplatin_mb_most_sensitive_1se)

cisplatin_mb_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_mb_idx]
actual_cisplatin_mb_least <- cisplatin_test$least_sensitive[cisplatin_mb_idx]
cisplatin_mb_least_1se_auc <- auc(actual_cisplatin_mb_least, cisplatin_mb_least_sensitive_1se)

# MESO
cisplatin_meso_idx <- cisplatin_test$TCGA_class == 'MESO' #5

cisplatin_meso_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_meso_idx]
actual_cisplatin_meso_most         <- cisplatin_test$most_sensitive[cisplatin_meso_idx]
cisplatin_meso_most_min_auc        <- auc(actual_cisplatin_meso_most, cisplatin_meso_most_sensitive_min)

cisplatin_meso_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_meso_idx]
actual_cisplatin_meso_least <- cisplatin_test$least_sensitive[cisplatin_meso_idx]
cisplatin_meso_least_min_auc <- auc(actual_cisplatin_meso_least, cisplatin_meso_least_sensitive_min)

cisplatin_meso_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_meso_idx]
actual_cisplatin_meso_most <- cisplatin_test$most_sensitive[cisplatin_meso_idx]
cisplatin_meso_most_1se_auc <- auc(actual_cisplatin_meso_most, cisplatin_meso_most_sensitive_1se)

cisplatin_meso_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_meso_idx]
actual_cisplatin_meso_least <- cisplatin_test$least_sensitive[cisplatin_meso_idx]
cisplatin_meso_least_1se_auc <- auc(actual_cisplatin_meso_least, cisplatin_meso_least_sensitive_1se)

# NB
cisplatin_nb_idx <- cisplatin_test$TCGA_class == 'NB' #13

cisplatin_nb_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_nb_idx]
actual_cisplatin_nb_most         <- cisplatin_test$most_sensitive[cisplatin_nb_idx]
cisplatin_nb_most_min_auc        <- auc(actual_cisplatin_nb_most, cisplatin_nb_most_sensitive_min)

cisplatin_nb_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_nb_idx]
actual_cisplatin_nb_least <- cisplatin_test$least_sensitive[cisplatin_nb_idx]
cisplatin_nb_least_min_auc <- auc(actual_cisplatin_nb_least, cisplatin_nb_least_sensitive_min)

cisplatin_nb_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_nb_idx]
actual_cisplatin_nb_most <- cisplatin_test$most_sensitive[cisplatin_nb_idx]
cisplatin_nb_most_1se_auc <- auc(actual_cisplatin_nb_most, cisplatin_nb_most_sensitive_1se)

cisplatin_nb_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_nb_idx]
actual_cisplatin_nb_least <- cisplatin_test$least_sensitive[cisplatin_nb_idx]
cisplatin_nb_least_1se_auc <- auc(actual_cisplatin_nb_least, cisplatin_nb_least_sensitive_1se)

# OV
cisplatin_ov_idx <- cisplatin_test$TCGA_class == 'OV' #17

cisplatin_ov_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_ov_idx]
actual_cisplatin_ov_most         <- cisplatin_test$most_sensitive[cisplatin_ov_idx]
cisplatin_ov_most_min_auc        <- auc(actual_cisplatin_ov_most, cisplatin_ov_most_sensitive_min)

cisplatin_ov_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_ov_idx]
actual_cisplatin_ov_least <- cisplatin_test$least_sensitive[cisplatin_ov_idx]
cisplatin_ov_least_min_auc <- auc(actual_cisplatin_ov_least, cisplatin_ov_least_sensitive_min)

cisplatin_ov_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_ov_idx]
actual_cisplatin_ov_most <- cisplatin_test$most_sensitive[cisplatin_ov_idx]
cisplatin_ov_most_1se_auc <- auc(actual_cisplatin_ov_most, cisplatin_ov_most_sensitive_1se)

cisplatin_ov_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_ov_idx]
actual_cisplatin_ov_least <- cisplatin_test$least_sensitive[cisplatin_ov_idx]
cisplatin_ov_least_1se_auc <- auc(actual_cisplatin_ov_least, cisplatin_ov_least_sensitive_1se)

# PAAD
cisplatin_paad_idx <- cisplatin_test$TCGA_class == 'PAAD' #13

cisplatin_paad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_paad_idx]
actual_cisplatin_paad_most         <- cisplatin_test$most_sensitive[cisplatin_paad_idx]
cisplatin_paad_most_min_auc        <- auc(actual_cisplatin_paad_most, cisplatin_paad_most_sensitive_min)

cisplatin_paad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_paad_idx]
actual_cisplatin_paad_least <- cisplatin_test$least_sensitive[cisplatin_paad_idx]
cisplatin_paad_least_min_auc <- auc(actual_cisplatin_paad_least, cisplatin_paad_least_sensitive_min)

cisplatin_paad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_paad_idx]
actual_cisplatin_paad_most <- cisplatin_test$most_sensitive[cisplatin_paad_idx]
cisplatin_paad_most_1se_auc <- auc(actual_cisplatin_paad_most, cisplatin_paad_most_sensitive_1se)

cisplatin_paad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_paad_idx]
actual_cisplatin_paad_least <- cisplatin_test$least_sensitive[cisplatin_paad_idx]
cisplatin_paad_least_1se_auc <- auc(actual_cisplatin_paad_least, cisplatin_paad_least_sensitive_1se)

# SCLC
cisplatin_sclc_idx <- cisplatin_test$TCGA_class == 'SCLC' #31

cisplatin_sclc_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_sclc_idx]
actual_cisplatin_sclc_most         <- cisplatin_test$most_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_most_min_auc        <- auc(actual_cisplatin_sclc_most, cisplatin_sclc_most_sensitive_min)

cisplatin_sclc_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_sclc_idx]
actual_cisplatin_sclc_least <- cisplatin_test$least_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_least_min_auc <- auc(actual_cisplatin_sclc_least, cisplatin_sclc_least_sensitive_min)

cisplatin_sclc_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_sclc_idx]
actual_cisplatin_sclc_most <- cisplatin_test$most_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_most_1se_auc <- auc(actual_cisplatin_sclc_most, cisplatin_sclc_most_sensitive_1se)

cisplatin_sclc_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_sclc_idx]
actual_cisplatin_sclc_least <- cisplatin_test$least_sensitive[cisplatin_sclc_idx]
cisplatin_sclc_least_1se_auc <- auc(actual_cisplatin_sclc_least, cisplatin_sclc_least_sensitive_1se)

# SKCM
cisplatin_skcm_idx <- cisplatin_test$TCGA_class == 'SKCM' #26

cisplatin_skcm_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_skcm_idx]
actual_cisplatin_skcm_most         <- cisplatin_test$most_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_most_min_auc        <- auc(actual_cisplatin_skcm_most, cisplatin_skcm_most_sensitive_min)

cisplatin_skcm_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_skcm_idx]
actual_cisplatin_skcm_least <- cisplatin_test$least_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_least_min_auc <- auc(actual_cisplatin_skcm_least, cisplatin_skcm_least_sensitive_min)

cisplatin_skcm_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_skcm_idx]
actual_cisplatin_skcm_most <- cisplatin_test$most_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_most_1se_auc <- auc(actual_cisplatin_skcm_most, cisplatin_skcm_most_sensitive_1se)

cisplatin_skcm_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_skcm_idx]
actual_cisplatin_skcm_least <- cisplatin_test$least_sensitive[cisplatin_skcm_idx]
cisplatin_skcm_least_1se_auc <- auc(actual_cisplatin_skcm_least, cisplatin_skcm_least_sensitive_1se)

# STAD
cisplatin_stad_idx <- cisplatin_test$TCGA_class == 'STAD' #8

cisplatin_stad_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_stad_idx]
actual_cisplatin_stad_most         <- cisplatin_test$most_sensitive[cisplatin_stad_idx]
cisplatin_stad_most_min_auc        <- auc(actual_cisplatin_stad_most, cisplatin_stad_most_sensitive_min)

cisplatin_stad_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_stad_idx]
actual_cisplatin_stad_least <- cisplatin_test$least_sensitive[cisplatin_stad_idx]
cisplatin_stad_least_min_auc <- auc(actual_cisplatin_stad_least, cisplatin_stad_least_sensitive_min)

cisplatin_stad_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_stad_idx]
actual_cisplatin_stad_most <- cisplatin_test$most_sensitive[cisplatin_stad_idx]
cisplatin_stad_most_1se_auc <- auc(actual_cisplatin_stad_most, cisplatin_stad_most_sensitive_1se)

cisplatin_stad_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_stad_idx]
actual_cisplatin_stad_least <- cisplatin_test$least_sensitive[cisplatin_stad_idx]
cisplatin_stad_least_1se_auc <- auc(actual_cisplatin_stad_least, cisplatin_stad_least_sensitive_1se)

# THCA
cisplatin_thca_idx <- cisplatin_test$TCGA_class == 'thca' #4

cisplatin_thca_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_thca_idx]
actual_cisplatin_thca_most         <- cisplatin_test$most_sensitive[cisplatin_thca_idx]
cisplatin_thca_most_min_auc        <- auc(actual_cisplatin_thca_most, cisplatin_thca_most_sensitive_min)

cisplatin_thca_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_thca_idx]
actual_cisplatin_thca_least <- cisplatin_test$least_sensitive[cisplatin_thca_idx]
cisplatin_thca_least_min_auc <- auc(actual_cisplatin_thca_least, cisplatin_thca_least_sensitive_min)

cisplatin_thca_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_thca_idx]
actual_cisplatin_thca_most <- cisplatin_test$most_sensitive[cisplatin_thca_idx]
cisplatin_thca_most_1se_auc <- auc(actual_cisplatin_thca_most, cisplatin_thca_most_sensitive_1se)

cisplatin_thca_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_thca_idx]
actual_cisplatin_thca_least <- cisplatin_test$least_sensitive[cisplatin_thca_idx]
cisplatin_thca_least_1se_auc <- auc(actual_cisplatin_thca_least, cisplatin_thca_least_sensitive_1se)

# UCEC
cisplatin_ucec_idx <- cisplatin_test$TCGA_class == 'UCEC' #4

cisplatin_ucec_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_ucec_idx]
actual_cisplatin_ucec_most         <- cisplatin_test$most_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_most_min_auc        <- auc(actual_cisplatin_ucec_most, cisplatin_ucec_most_sensitive_min)

cisplatin_ucec_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_ucec_idx]
actual_cisplatin_ucec_least <- cisplatin_test$least_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_least_min_auc <- auc(actual_cisplatin_ucec_least, cisplatin_ucec_least_sensitive_min)

cisplatin_ucec_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_ucec_idx]
actual_cisplatin_ucec_most <- cisplatin_test$most_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_most_1se_auc <- auc(actual_cisplatin_ucec_most, cisplatin_ucec_most_sensitive_1se)

cisplatin_ucec_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_ucec_idx]
actual_cisplatin_ucec_least <- cisplatin_test$least_sensitive[cisplatin_ucec_idx]
cisplatin_ucec_least_1se_auc <- auc(actual_cisplatin_ucec_least, cisplatin_ucec_least_sensitive_1se)

# UNCLASSIFIED
cisplatin_unclassified_idx <- cisplatin_test$TCGA_class == 'UNCLASSIFIED' #52

cisplatin_unclassified_most_sensitive_min  <- new_cisplatin_most_sensitive_min[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_most         <- cisplatin_test$most_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_most_min_auc        <- auc(actual_cisplatin_unclassified_most, cisplatin_unclassified_most_sensitive_min)

cisplatin_unclassified_least_sensitive_min <- new_cisplatin_least_sensitive_min[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_least <- cisplatin_test$least_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_least_min_auc <- auc(actual_cisplatin_unclassified_least, cisplatin_unclassified_least_sensitive_min)

cisplatin_unclassified_most_sensitive_1se <- new_cisplatin_most_sensitive_1se[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_most <- cisplatin_test$most_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_most_1se_auc <- auc(actual_cisplatin_unclassified_most, cisplatin_unclassified_most_sensitive_1se)

cisplatin_unclassified_least_sensitive_1se <- new_cisplatin_least_sensitive_1se[cisplatin_unclassified_idx]
actual_cisplatin_unclassified_least <- cisplatin_test$least_sensitive[cisplatin_unclassified_idx]
cisplatin_unclassified_least_1se_auc <- auc(actual_cisplatin_unclassified_least, cisplatin_unclassified_least_sensitive_1se)

#### graphing test auc values for overall and by cancer type --------
#graphing auc of each drug model overall and by cancer type ##
#have this model come from GLM_models folder
model_stats <- read.csv('model_stats_final.csv', header = TRUE, stringsAsFactors = FALSE, row.names = 1)

model_stats <- as.matrix(model_stats)

colors <- colorRampPalette(c("blue", "white", "red"))(100)

png(filename = 'Images/something.png', width = 480, height = 480, units = px)

test_auc_heat <- heatmap.2(model_stats, trace = 'none', Rowv = FALSE, Colv = FALSE, col = colors, key = FALSE, cexRow = 0.9, cexCol = 0.8, cellnote = model_stats, notecol = 'black', rowsep = 1, sepwidth = c(0.1,0.1), margins = c(8,8))
plot(test_auc_heat)
dev.off()

### graphing auc curves for gdsc testing results ---------
pred <- prediction(new_bleo_least_sensitive_min[bleo_blca_idx], bleomycin_test$least_sensitive[bleo_blca_idx])
perf <- performance(pred, 'tpr', 'fpr')
png(file = '')
plot(perf, main = '')

auc <- performance(pred, 'auc')
auc <- unlist(slot(auc, 'y.values'))
minauc <- min(round(auc, digits = 2))
maxauc <- max(round(auc, digits = 2))
minauct <- paste(c('minAUC = '), minauc, sep = '')
maxauct <- paste(c("max(AUC) = "),maxauc,sep="")
legend(0.3,0.6,c(minauct,maxauct,"\n"),border="white",cex=1.7,box.col = "white")
dev.off()