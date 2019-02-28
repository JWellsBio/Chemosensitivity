# THIS SCRIPT GETS THE RAW GENE EXPRESSION FILES AND CONVERTS THEM TO USABLE FORMATS FOR MODEL BUILDING

### load (and install) necessary packages -----
if (!require ('openxlsx')) install.packages('openxlsx')
library(openxlsx) # to read in excel files

if (!require ('biomaRt')) BiocInstaller::biocLite('biomaRt')
library(biomaRt) # for conversion between gene names, Entrez IDs, and ENSG IDs

### load raw expression files ------
### CELL LINES
# GDSC
gdsc_rna_seq <- read.csv('Gene_Expression_Files/gdsc_rnaseq.csv', header = TRUE, stringsAsFactors = FALSE)
rownames(gdsc_rna_seq) <- make.names(gdsc_rna_seq$ï.., unique = TRUE) # ENSG IDs as row names
gdsc_rna_seq <- gdsc_rna_seq[, -1] # remove gene names
gdsc_rna_seq <- gdsc_rna_seq[-c(17738,17739), ] # removes NAs that were there for some reason
colnames(gdsc_rna_seq) <- gsub('X', '', colnames(gdsc_rna_seq))
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

#HNSC
hnsc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-HNSC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(hnsc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(hnsc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 546 cases (TCGA ID)

#KIRC
kirc_tcga_rnaseq <- read.delim('TCGA_files/TCGA-KIRC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(kirc_tcga_rnaseq) <- gsub('\\..*$', '', rownames(kirc_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 607 cases (TCGA ID)

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

#UCEC
ucec_tcga_rnaseq <- read.delim('TCGA_files/TCGA-UCEC.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(ucec_tcga_rnaseq) <- gsub('\\..*$', '', rownames(ucec_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 583 cases (TCGA ID)

#UCS
ucs_tcga_rnaseq <- read.delim('TCGA_files/TCGA-UCS.htseq_counts.tsv', sep = '\t', row.names = 1, stringsAsFactors = FALSE)
rownames(ucs_tcga_rnaseq) <- gsub('\\..*$', '', rownames(ucs_tcga_rnaseq))
# dimensions now: 60483 genes (by ENSG ID) x 56 cases (TCGA ID)

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
#14209 ENSG ID in common

gdsc_rna_seq <- gdsc_rna_seq[which(rownames(gdsc_rna_seq) %in% ensg_i_need), ]
#14209 x 962

ccle_microarray <- ccle_microarray[which(rownames(ccle_microarray) %in% ensg_i_need), ]
#14209 x 855

nci60_microarray <- nci60_microarray[which(rownames(nci60_microarray) %in% ensg_i_need), ]
#14209 x 52

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
write.csv(nci60_microarray, file = 'Processed_Gene_Expression/nci60_microarray_processed.csv')
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
