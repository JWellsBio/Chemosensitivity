# THIS SCRIPT GETS THE RAW GENE EXPRESSION FILES AND CONVERTS THEM TO USABLE FORMATS FOR MODEL BUILDING

### load (and install) necessary packages -----
if (!require ('openxlsx')) install.packages('openxlsx')
library(openxlsx) # to read in excel files

if (!require ('biomaRt')) BiocInstaller::biocLite('biomaRt')
library(biomaRt) # for conversion between gene names, Entrez IDs, and ENSG IDs


## functions needed ----
#process TCGA expression files
tcga_gene_exp <- function(tcga_class) {
  tcga_rna_seq <- read.delim(paste0('TCGA_files/TCGA-', tcga_class, '.htseq_counts.tsv'), sep = '\t', row.names = 1, stringsAsFactors = FALSE)
  rownames(tcga_rna_seq) <- gsub('\\..*$', '', rownames(tcga_rna_seq))
  return(tcga_rna_seq)
}


### load raw expression files ------
### CELL LINES
# GDSC
gdsc_rna_seq <- read.csv('Gene_Expression_Files/gdsc_rnaseq.csv', header = TRUE, stringsAsFactors = FALSE)
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$?.., unique = TRUE) # ENSG IDs as row names
gdsc_rna_seq <- gdsc_rna_seq[, -1] # remove gene names
gdsc_rna_seq <- gdsc_rna_seq[-c(17738,17739), ] # removes NAs that were there due to previous use in excel
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))
# dimensions now: 17737 genes (by ENSG ID) x 962 cell lines (by COSMIC ID)


### TUMOR SETS (TCGA)
# BLCA
blca_tcga_rnaseq <- tcga_gene_exp('BLCA')

#BRCA
brca_tcga_rnaseq <- tcga_gene_exp('BRCA')

#CESC
cesc_tcga_rnaseq <- tcga_gene_exp('CESC')

#CHOL
chol_tcga_rnaseq <- tcga_gene_exp('CHOL')

#COAD
coad_tcga_rnaseq <- tcga_gene_exp('COAD')

#DLBC
dlbc_tcga_rnaseq <- tcga_gene_exp('DLBC')

#ESCA
esca_tcga_rnaseq <- tcga_gene_exp('ESCA')

#HNSC
hnsc_tcga_rnaseq <- tcga_gene_exp('HNSC')

#KIRC
kirc_tcga_rnaseq <- tcga_gene_exp('KIRC')

#LIHC
lihc_tcga_rnaseq <- tcga_gene_exp('LIHC')

#LUAD
luad_tcga_rnaseq <- tcga_gene_exp('LUAD')

#LUSC
lusc_tcga_rnaseq <- tcga_gene_exp('LUSC')

#MESO
meso_tcga_rnaseq <- tcga_gene_exp('MESO')

#OV
ov_tcga_rnaseq   <- tcga_gene_exp('OV')

#PAAD
paad_tcga_rnaseq <- tcga_gene_exp('PAAD')

#READ
read_tcga_rnaseq <- tcga_gene_exp('READ')

#SARC
sarc_tcga_rnaseq <- tcga_gene_exp('SARC')

#SKCM
skcm_tcga_rnaseq <- tcga_gene_exp('SKCM')

#STAD
stad_tcga_rnaseq <- tcga_gene_exp('STAD')

#TGCT
tgct_tcga_rnaseq <- tcga_gene_exp("TGCT")

#UCEC
ucec_tcga_rnaseq <- tcga_gene_exp('UCEC')

#UCS
ucs_tcga_rnaseq  <- tcga_gene_exp('UCS')

## make sure all gene expression sets have the same features ------
# makes sure model doesn't include something not in another set
# gdsc is by ENSG ID
# ccle is by gene name
# TCGA is by ENSG ID

# get gene names and entrez IDs
ccle_genes <- rownames(ccle_microarray)

# convert them to ENSG IDs
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
ccle_annots <- getBM(mart = mart, attributes = c('hgnc_symbol', 'ensembl_gene_id'), filters = 'hgnc_symbol', values = ccle_genes, uniqueRows = TRUE)
ccle_annots <- ccle_annots[!duplicated(ccle_annots[, 2]), ]
ccle_microarray <- ccle_microarray[which(rownames(ccle_microarray) %in% ccle_annots[, 1]), ]
ccle_annots <- ccle_annots[which(ccle_annots[, 1] %in% rownames(ccle_microarray)), ]
ccle_microarray <- ccle_microarray[match(ccle_annots[, 1], rownames(ccle_microarray)), ]
rownames(ccle_microarray) <- make.names(ccle_annots[, 2], unique = TRUE)
# dimensions now: 17846 genes (by ENSG ID) x 855 cell lines (name and type)


# limit each gene expression set to common ENSG IDs
ensg_i_need <- Reduce(intersect, list(rownames(gdsc_rna_seq), rownames(ccle_microarray)))
#14209 ENSG ID in common

gdsc_rna_seq <- gdsc_rna_seq[which(rownames(gdsc_rna_seq) %in% ensg_i_need), ]
#14209 x 962

ccle_microarray <- ccle_microarray[which(rownames(ccle_microarray) %in% ensg_i_need), ]
#14209 x 855


blca_tcga_rnaseq <- blca_tcga_rnaseq[which(rownames(blca_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 430

brca_tcga_rnaseq <- brca_tcga_rnaseq[which(rownames(brca_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 1217

cesc_tcga_rnaseq <- cesc_tcga_rnaseq[which(rownames(cesc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 309

chol_tcga_rnaseq <- chol_tcga_rnaseq[which(rownames(chol_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 45

coad_tcga_rnaseq <- coad_tcga_rnaseq[which(rownames(coad_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 512

dlbc_tcga_rnaseq <- dlbc_tcga_rnaseq[which(rownames(dlbc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 48

esca_tcga_rnaseq <- esca_tcga_rnaseq[which(rownames(esca_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 173

hnsc_tcga_rnaseq <- hnsc_tcga_rnaseq[which(rownames(hnsc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 546

kirc_tcga_rnaseq <- kirc_tcga_rnaseq[which(rownames(kirc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 607

lihc_tcga_rnaseq <- lihc_tcga_rnaseq[which(rownames(lihc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 424

luad_tcga_rnaseq <- luad_tcga_rnaseq[which(rownames(luad_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 585

lusc_tcga_rnaseq <- lusc_tcga_rnaseq[which(rownames(lusc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 550

meso_tcga_rnaseq <- meso_tcga_rnaseq[which(rownames(meso_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 86

ov_tcga_rnaseq <- ov_tcga_rnaseq[which(rownames(ov_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 379

paad_tcga_rnaseq <- paad_tcga_rnaseq[which(rownames(paad_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 182

read_tcga_rnaseq <- read_tcga_rnaseq[which(rownames(read_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 177

sarc_tcga_rnaseq <- sarc_tcga_rnaseq[which(rownames(sarc_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 265

skcm_tcga_rnaseq <- skcm_tcga_rnaseq[which(rownames(skcm_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 472

stad_tcga_rnaseq <- stad_tcga_rnaseq[which(rownames(stad_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 407

tgct_tcga_rnaseq <- tgct_tcga_rnaseq[which(rownames(tgct_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 156

ucec_tcga_rnaseq <- ucec_tcga_rnaseq[which(rownames(ucec_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 583

ucs_tcga_rnaseq <- ucs_tcga_rnaseq[which(rownames(ucs_tcga_rnaseq) %in% ensg_i_need), ]
#14209 x 56

## write these to files for later ------
write.csv(gdsc_rna_seq, file = 'Processed_Gene_Expression/gdsc_rna_seq_processed.csv')
write.csv(ccle_microarray, file = 'Processed_Gene_Expression/ccle_microarray_processed.csv')
write.csv(blca_tcga_rnaseq, file = 'Processed_Gene_Expression/blca_tcga_rna_seq_processed.csv')
write.csv(brca_tcga_rnaseq, file = 'Processed_Gene_Expression/brca_tcga_rna_seq_processed.csv')
write.csv(cesc_tcga_rnaseq, file = 'Processed_Gene_Expression/cesc_tcga_rna_seq_processed.csv')
write.csv(chol_tcga_rnaseq, file = 'Processed_Gene_Expression/chol_tcga_rna_seq_processed.csv')
write.csv(coad_tcga_rnaseq, file = 'Processed_Gene_Expression/coad_tcga_rna_seq_processed.csv')
write.csv(dlbc_tcga_rnaseq, file = 'Processed_Gene_Expression/dlbc_tcga_rna_seq_processed.csv')
write.csv(esca_tcga_rnaseq, file = 'Processed_Gene_Expression/esca_tcga_rna_seq_processed.csv')
write.csv(hnsc_tcga_rnaseq, file = 'Processed_Gene_Expression/hnsc_tcga_rna_seq_processed.csv')
write.csv(kirc_tcga_rnaseq, file = 'Processed_Gene_Expression/kirc_tcga_rna_seq_processed.csv')
write.csv(lihc_tcga_rnaseq, file = 'Processed_Gene_Expression/lihc_tcga_rna_seq_processed.csv')
write.csv(luad_tcga_rnaseq, file = 'Processed_Gene_Expression/luad_tcga_rna_seq_processed.csv')
write.csv(lusc_tcga_rnaseq, file = 'Processed_Gene_Expression/lusc_tcga_rna_seq_processed.csv')
write.csv(meso_tcga_rnaseq, file = 'Processed_Gene_Expression/meso_tcga_rna_seq_processed.csv')
write.csv(ov_tcga_rnaseq, file = 'Processed_Gene_Expression/ov_tcga_rna_seq_processed.csv')
write.csv(paad_tcga_rnaseq, file = 'Processed_Gene_Expression/paad_tcga_rna_seq_processed.csv')
write.csv(read_tcga_rnaseq, file = 'Processed_Gene_Expression/read_tcga_rna_seq_processed.csv')
write.csv(sarc_tcga_rnaseq, file = 'Processed_Gene_Expression/sarc_tcga_rna_seq_processed.csv')
write.csv(skcm_tcga_rnaseq, file = 'Processed_Gene_Expression/skcm_tcga_rna_seq_processed.csv')
write.csv(stad_tcga_rnaseq, file = 'Processed_Gene_Expression/stad_tcga_rna_seq_processed.csv')
write.csv(tgct_tcga_rnaseq, file = 'Processed_Gene_Expression/tgct_tcga_rna_seq_processed.csv')
write.csv(ucec_tcga_rnaseq, file = 'Processed_Gene_Expression/ucec_tcga_rna_seq_processed.csv')
write.csv(ucs_tcga_rnaseq, file = 'Processed_Gene_Expression/ucs_tcga_rna_seq_processed.csv')
