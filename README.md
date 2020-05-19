# Gene-network-analysis

This repository has been created to present the work done on computational analysis of fused co-expression networks for the identification of candidate cancer gene biomarkers. 

The jupyther notebook Pipeline.ipynb presents these results.

# What can you find in this repository?
This repository contains all data, scripts and results related to LIHC tumor.
In particular, you will find:
- 4 .py files which contain the main steps of the process,
- two folders Datasets and Datasets_created that contains the data retrieved from GMQL and PyGMQL,
- one folder Matrices that contains the matrices computed after having executed preprocessing.py,
- one folder Networks that contains the networks built with creation_matrices.py and fusion_networks.py,
- one folder Graphs that contains the graph of the network in the format required by Gephi,
- one folder Extracted that contains the data regarding all the genes extracted with the pipeline,
- one file matrix_GTypeGName.xlsx that contains the correspondances between the gene symbols and the gene type.

#How to run the notebook
pip install -r requirements.txt

Execute the jupyter notebook Pipeline.ipynb until the part 'After extraction of communities with Gephi'.
Use Gephi in order to extract the relevant communities and save the genes in the folder Extracted with the name 'comm_'+str(tumor)+'.csv'.
