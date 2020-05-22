# Gene-network-analysis

This repository has been created to present the work done on computational analysis of fused co-expression networks for the identification of candidate cancer gene biomarkers. 

The jupyther notebook Pipeline.ipynb contains the complete pipeline for the construction of the fused co-expression networks and for the extraction of relevant gene biomarkers.
The jupyter notebook Knowledge-based evaluation of the results.ipynb contains systematic and statistic evalution the the extracted genes (here called IC genes). 

# What can you find in this repository?
This repository contains all data, scripts and results related to LIHC tumor.
In particular, you will find:
- 4 .py files which contain the main steps of the process,
- two folders Datasets and Datasets_created that contains the data retrieved from GMQL and PyGMQL,
- one folder Graphs that contains the graph of the network in the format required by Gephi,
- one folder Extracted that contains the data regarding all the genes extracted with the pipeline,
- one folder Auc_acc_f1 that contains the tables with the performances values after having executed extraction_classification.py,
- one folder Boxplots that contains the boxplots with the comparisons of the performances of the classification,
- one file matrix_GTypeGName.xlsx that contains the correspondances between the gene symbols and the gene type.
- one folder Results containing the gene symbols of the LIHC fused network, the IC gene symbols and the PubMed evaluation results presentd in Knowledge-based evaluation of the results.ipynb.
- one folder Supplementary data containing the DB_pharmacologically_active.csv file downloaded from DrugBank to be used in the DrugBank evaluation part of Knowledge-based evaluation of the results.ipynb. 

# How to run the notebook
pip install -r requirements.txt

Execute the jupyter notebook Pipeline.ipynb until the part 'After extraction of communities with Gephi'.
Use Gephi in order to extract the relevant communities and save the genes in the folder Extracted with the name 'comm_'+str(tumor)+'.csv'.
