# Chemosensitivity
Chemosensitivity Project Scripts and Files
**NO GENE EXPRESSION FILES ARE LOADED AT THIS TIME**

## Idea for Project
There plenty of projects that predict sensitivity on a small scale, using only one drug and/or one cancer type. We sought to go as big as possible and build pan-cancer models to predict sensitivity for as many chemotherapy drugs as possible. One factor we considered was that it had to be able to be sequenced fairly inexpensively. The standard here would be Nanostring, which limits us to ~800 genes.

## Study Design
The flowchart for the study design can be found in the [Images](https://github.com/JWellsBio/Chemosensitivity/tree/master/Images) folder. 
+ Briefly, we sought to use only publicly available gene expression data to build our predictive models. The best sources of gene expression data paired with sensitivity data are the [Genomics of Drug Sensitivity](https://www.cancerrxgene.org/) database and the [Cancer Cell Line Encyclopedia](https://portals.broadinstitute.org/ccle).
+ After splitting into testing and training data, generalized linear models (GLMs) were developed per drug. Because trying to predict a continuous variable (such as IC50 or AUC) is difficult, sensitivity measures were binned into the top 20% (most sensitive) and the bottom 20% (least sensitive).
+ Models were evaluated based on AUC (10x-cross-validated within the training data) and were tuned to have the fewest genes possible for the highest AUC.
+ Finalized pan-cancer models were then applied to withheld testing data and AUC was measured.
+ Further, pan-cancer models were tested against specific tissue types in the testing data to see how well pan-cancer models predicted sensitivity to individual tissue types.
+ Finally, drug models were tested on human tumor datasets from [The Cancer Genome Atlas] (https://cancergenome.nih.gov/).
