# Chemosensitivity
Chemosensitivity Project Scripts and Files
**NO GENE EXPRESSION FILES ARE LOADED AT THIS TIME**

## Idea for Project
There are plenty of projects that predict sensitivity on a small scale, using only one drug and/or one cancer type. We sought to go as big as possible and build pan-cancer models to predict sensitivity for as many chemotherapy drugs as possible. One factor we considered was that it had to be able to be sequenced fairly inexpensively. The standard here would be Nanostring, which limits us to ~800 genes.

## Study Design
The flowchart for the study design can be found in the [Images](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/study_design.png) folder. 
+ Briefly, we sought to use only publicly available gene expression data to build our predictive models. The best sources of gene expression data paired with sensitivity data are the [Genomics of Drug Sensitivity](https://www.cancerrxgene.org/) database and the [Cancer Cell Line Encyclopedia](https://portals.broadinstitute.org/ccle).
+ After splitting into testing and training data, generalized linear models (GLMs) were developed per drug. Because trying to predict a continuous variable (such as IC50 or AUC) is difficult, sensitivity measures were binned into the top 20% (most sensitive) and the bottom 20% (least sensitive).
+ Models were evaluated based on AUC (10x-cross-validated within the training data) and were tuned to have the fewest genes possible for the highest AUC. [Cisplatin example](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/cisplatin_most_auc.png)
+ Finalized pan-cancer models were then applied to withheld testing data and AUC was measured.
+ Further, pan-cancer models were tested against specific tissue types in the testing data to see how well pan-cancer models predicted sensitivity to individual tissue types.
+ Finally, drug models were tested on human tumor datasets from [The Cancer Genome Atlas](https://cancergenome.nih.gov/).

## Results for pan-cancer models applied to testing data
+ Pan-cancer models built from GDSC data can be seen [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/GDSC_AUC_heatmap.png)
+ Pan-cancer models built from CCLE data can be seen [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/CCLE_AUC_heatmap.png)

## Total gene numbers
+ The total gene set for all drugs considered comes out to 248. A break-down of number of genes needed per drug model can be seen [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/geneset_treemap.png)
+ A table detailing genes needed in GDSC models can be seen [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/gdsc_gene_table.PNG)
+ And for CCLE [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/ccle_gene_table.PNG)

## TCGA results **STILL SLIGHTLY IN PROGRESS**
+ Pan-cancer models were tested against available TCGA class/drug combinations where n > 5. Sensitivity measure used was RFS.
+ ROC curves were generated for every combination tested. An exmaple can be found [here](https://github.com/JWellsBio/Chemosensitivity/blob/master/Images/ov_cisplatin_tcga_gdsc_auc.png)
+ A massive table with all combinaitons tested is still in progress.

## Steps to repeat analysis
+ Run this script
+ Run this script
+ Run this script
+ Blah blah blah
