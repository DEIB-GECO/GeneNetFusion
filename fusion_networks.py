
import pandas as pd
import numpy as np
import warnings
from scipy.stats import binom
import snf
import os
#from creation_matrices import adj_to_list
os.environ["MKL_NUM_THREADS"] = "10"
os.environ["NUMEXPR_NUM_THREADS"] = "10"
os.environ["OMP_NUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["VECLIB_MAXIMUM_THREADS"] = "10"

warnings.filterwarnings('ignore')

#Takes the adjacency matrix and transform it into a list of edges and saves it into output_filename
def adj_to_list(A):
    List=[]
    A=pd.DataFrame(A)
    for source in A.index.values:
        for target in A.index.values:
            if source>target:
                if A[source][target]!=0:
                    List.append((target,source,A[source][target]))
    List1=pd.DataFrame(List)
    return List1
   
#It allows to read a file (as a Dataframe) that has 3 or more columns (sep: tab) nd one line header with column names
# ('src' and 'trg' that indicates the source and destination, the column of interest indicate the weights)
#It return the parsed network data, the number of nodes in the network and the number of edges.
def read(filename, column_of_interest, consider_self_loops = True, drop_zeroes = True):
   table = pd.read_csv(filename, sep = '\t')
   table = table[["src", "trg", column_of_interest]]
   table.rename(columns = {column_of_interest: "nij"}, inplace = True)

   if drop_zeroes:
      table = table[table["nij"] > 0]
   if not consider_self_loops:
      table = table[table["src"] != table["trg"]]

   original_nodes = len(set(table["src"]) | set(table["trg"]))
   original_edges = table.shape[0]
   return table, original_nodes, original_edges

#Reads a preprocessed edge table and returns only the edges supassing a significance threshold.
def thresholding(table, threshold):
   table = table.copy()
   if "sdev_cij" in table:
      return table[(table["score"] - (threshold * table["sdev_cij"])) > 0][["src", "trg", "nij", "score"]]
   else:
      return table[table["score"] > threshold][["src", "trg", "nij", "score"]]

#It allows to correct the noise of the  table
def noise_corrected(table, return_self_loops = False, calculate_p_value = False):
   table = table.copy()
   src_sum = table.groupby(by = "src").sum()[["nij"]]
   table = table.merge(src_sum, left_on = "src", right_index = True, suffixes = ("", "_src_sum"))
   trg_sum = table.groupby(by = "trg").sum()[["nij"]]
   table = table.merge(trg_sum, left_on = "trg", right_index = True, suffixes = ("", "_trg_sum"))
   table.rename(columns = {"nij_src_sum": "ni.", "nij_trg_sum": "n.j"}, inplace = True)
   table["n.."] = table["nij"].sum()
   table["mean_prior_probability"] = ((table["ni."] * table["n.j"]) / table["n.."]) * (1 / table["n.."])

   if calculate_p_value:
      table["score"] = binom.cdf(table["nij"], table["n.."], table["mean_prior_probability"])
      return table[["src", "trg", "nij", "score"]]

   table["kappa"] = table["n.."] / (table["ni."] * table["n.j"])
   table["score"] = ((table["kappa"] * table["nij"]) - 1) / ((table["kappa"] * table["nij"]) + 1)
   table["var_prior_probability"] = (1 / (table["n.."] ** 2)) * (table["ni."] * table["n.j"] * (table["n.."] - table["ni."]) * (table["n.."] - table["n.j"])) / ((table["n.."] ** 2) * ((table["n.."] - 1)))
   table["alpha_prior"] = (((table["mean_prior_probability"] ** 2) / table["var_prior_probability"]) * (1 - table["mean_prior_probability"])) - table["mean_prior_probability"]
   table["beta_prior"] = (table["mean_prior_probability"] / table["var_prior_probability"]) * (1 - (table["mean_prior_probability"] ** 2)) - (1 - table["mean_prior_probability"])
   table["alpha_post"] = table["alpha_prior"] + table["nij"]
   table["beta_post"] = table["n.."] - table["nij"] + table["beta_prior"]
   table["expected_pij"] = table["alpha_post"] / (table["alpha_post"] + table["beta_post"])
   table["variance_nij"] = table["expected_pij"] * (1 - table["expected_pij"]) * table["n.."]
   table["d"] = (1.0 / (table["ni."] * table["n.j"])) - (table["n.."] * ((table["ni."] + table["n.j"]) / ((table["ni."] * table["n.j"]) ** 2)))
   table["variance_cij"] = table["variance_nij"] * (((2 * (table["kappa"] + (table["nij"] * table["d"]))) / (((table["kappa"] * table["nij"]) + 1) ** 2)) ** 2) 
   table["sdev_cij"] = table["variance_cij"] ** .5

   if not return_self_loops:
      table = table[table["src"] != table["trg"]]

   return table[["src", "trg", "nij", "score", "sdev_cij"]]

#It computes the disparity filter for a matrix
def disparity_filter(table, return_self_loops = False):

   table = table.copy()
   table_sum = table.groupby(table["src"]).sum().reset_index()
   table_deg = table.groupby(table["src"]).count()["trg"].reset_index()
   table = table.merge(table_sum, on = "src", how = "left", suffixes = ("", "_sum"))
   table = table.merge(table_deg, on = "src", how = "left", suffixes = ("", "_count"))
   table["score"] = 1.0 - ((1.0 - (table["nij"] / table["nij_sum"])) ** (table["trg_count"] - 1))
   table["variance"] = (table["trg_count"] ** 2) * (((20 + (4.0 * table["trg_count"])) / ((table["trg_count"] + 1.0) * (table["trg_count"] + 2) * (table["trg_count"] + 3))) - ((4.0) / ((table["trg_count"] + 1.0) ** 2)))

   if not return_self_loops:
      table = table[table["src"] != table["trg"]]

   return table[["src", "trg", "nij", "score", "variance"]]

#To pass from edgelist to adjacency matrix
def from_edgelist_to_matrix(A, length):
    List=np.zeros((length,length))
    for source in A.index:
            List[A['src'][source],A['trg'][source]]=A['nij'][source]
            List[A['trg'][source],A['src'][source]]=A['nij'][source]
    return List

class fusion_snf:
    def __init__(self, normal_sum, cancer_sum, tumor):
        self.tumor=tumor
        self.normal=normal_sum
        self.cancer=cancer_sum
        self.both=[self.normal,self.cancer]
        self.fused=pd.DataFrame(snf.snf(self.both,  K=20))
        self.fused_edgelist=adj_to_list(self.fused)
        self.fused_edgelist.columns=['src', 'trg', 'weight']
        self.fused_edgelist.to_csv('./Networks/fused_edgelist_'+ str(self.tumor)+'.csv', sep = '\t')
        
        self.fused_thr=self.threshold_snf()
        np.save('./Networks/fused_'+ str(self.tumor)+'.npy',self.fused_thr)

    #Threshold the fused matrix
    def threshold_snf(self):
        [fused, nodes, edges]=read('./Networks/fused_edgelist_'+ str(self.tumor)+'.csv', 'weight')
        fused=noise_corrected(fused)
        fused_df=disparity_filter(fused)
        quantiles = [99,99.5,99.8]
        for q in quantiles:
            thr=np.percentile(fused_df['score'].values, q)
            fused_dfThr=thresholding(fused_df, thr)

            #200000 is the threshold for Gephi to represent relevant edges
            if str(len(fused_dfThr))< 200000:
                break

        print('fused: ' + str(len(fused_dfThr)))
        fused=from_edgelist_to_matrix(fused_dfThr, len(self.normal))
        return fused
       
