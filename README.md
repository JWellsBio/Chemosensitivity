# Chemosensitivity
Chemosensitivity Project Scripts and Files

**GENE EXPRESSION FILES SHOULD BE DOWNLOADED FROM LINKS BELOW**

## Idea for Project
There are plenty of projects that predict sensitivity on a small scale, using only one drug and/or one cancer type. We sought to go as big as possible and build pan-cancer models to predict sensitivity for as many chemotherapy drugs as possible. One factor we considered was that it had to be able to be sequenced fairly inexpensively. The standard here would be Nanostring, which limits us to ~800 genes.

## Study Design
The flowchart for the study design can be found [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/study_design_synth.png). 
+ Briefly, we sought to use only publicly available gene expression data to build our predictive models. The best sources of gene expression data paired with sensitivity data are the [Genomics of Drug Sensitivity in Cancer](https://www.cancerrxgene.org/) database and the [Cancer Cell Line Encyclopedia](https://portals.broadinstitute.org/ccle).
+ Generalized linear models (GLMs) were developed per drug and tested against testing data. [Figure](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/Fig_2_testing.png)
+ Finalized pan-cancer models were then applied to withheld testing data and AUC was measured.
+ Further, pan-cancer models were tested against specific tissue types in the testing data to see how well pan-cancer models predicted sensitivity to individual tissue types.
+ Finally, drug models were tested on human tumor datasets from [The Cancer Genome Atlas](https://cancergenome.nih.gov/).


## TCGA results
+ Pan-cancer models were tested against available TCGA class/drug combinations where n >= 10. Sensitivity measure used was RFS.
+ Recurrence-free survival curves split based on predicted labels were generated and measured for significance by multivariate Cox regression analysis, adjusting for age, sex, and tumor stage. An example can be found [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/blca_cisplatin_surv.tiff)


