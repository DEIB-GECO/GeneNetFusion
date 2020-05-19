import pandas as pd
import numpy as np
import warnings
from scipy.special import betainc
from tqdm import tqdm
import os
os.environ["MKL_NUM_THREADS"] = "10"
os.environ["NUMEXPR_NUM_THREADS"] = "10"
os.environ["OMP_NUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["VECLIB_MAXIMUM_THREADS"] = "10"

warnings.filterwarnings('ignore')

#To calculate the strengths in an adjacency matrix
def strength(X):
    strength=[]
    for i in range(0, len(X)):
        strength.append(X[i].sum())
    return strength       

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
           

class Net:
    def __init__(self, normal, cancer, indices, tumor, normal_perc, cancer_perc):
        self.tumor =tumor
        self.normal=normal.set_index(np.arange(0,len(normal)))
        self.index_mirna_null_normal=self.normal[self.normal.duplicated()].index
        self.cancer=cancer.set_index(np.arange(0,len(cancer)))
        self.index_mirna_null_cancer=self.cancer[self.cancer.duplicated()].index
        self.indices=indices
        self.normal_eucl = self.eucl(self.normal)
        self.cancer_eucl = self.eucl(self.cancer)
        self.normal_pearson=self.pearson(self.normal, self.index_mirna_null_normal)
        self.cancer_pearson=self.pearson(self.cancer, self.index_mirna_null_cancer)
        self.normal_normalized=(self.normal-self.normal.mean())/self.normal.std()
        self.cancer_normalized=(self.cancer-self.cancer.mean())/self.cancer.std()
        self.normal_pearson_p_all, self.cancer_pearson_p_all = self.permutation(N=10)
        self.normal_pearson_thr =  self.threshold_permuted(self.normal_pearson, self.normal_pearson_p_all)
        self.cancer_pearson_thr =  self.threshold_permuted(self.cancer_pearson, self.cancer_pearson_p_all)
        np.save('./Networks/normal_pearson_permutation_'+str(tumor)+'.npy', self.normal_pearson_thr)
        np.save('./Networks/cancer_pearson_permutation_'+str(tumor)+'.npy', self.cancer_pearson_thr)
        self.normal_eucl_thr, self.cancer_eucl_thr = self.normalize_eucl() 
        self.normal_pearson_thr, self.cancer_pearson_thr = self.threshold_pvalue(normal_perc, cancer_perc)
        
        self.n_sum, self.c_sum = self.sum_pears_eucl()

    # It computes the Pearson's correlation among the rows of a matrix A, keeping the correlation of the submatrix B to 1
    # In our case the submatrix represents the mirnas whose expression in normal or cancer conditions has a mean equal to zero
    def pearson(self, A, B):
        Pearson = np.corrcoef(A)
        for i in B:
            for j in B:
                Pearson[i, j] = 1
        Pearson = np.nan_to_num(Pearson)
        return Pearson

    # It computes the Euclidean network starting from a matrix with variables as rows and observations as columns
    def eucl(self, X):
        values = X.values
        ds = []
        for i in tqdm(range(values.shape[0])):
            di = np.sqrt(np.power(values[i] - values, 2).sum(1))
            ds.append(di)

        eucl = np.vstack(ds)
        return eucl

    # It allows to permutate a matrix along both dimensions
    def permutate(self,X):
        dim = X.shape
        X_flat = X.flatten()
        X_flat = np.random.permutation(X_flat)
        return np.reshape(X_flat, dim)

    # It fills some elements of the matrix to one and put to 0 the nans
    def fill_null_mirna(self,index, matrix):
        for i in index:
            for j in index:
                matrix[i, j] = 1

        matrix = np.nan_to_num(matrix)
        return matrix

    # It computes the Pearson correlation and the relative pvalues and keeps only the triangular matrix of pvalues
    def corrcoef(self, matrix):
        r = np.corrcoef(matrix)
        rf = r[np.triu_indices(r.shape[0], 1)]
        df = matrix.shape[1] - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = betainc(0.5 * df, 0.5, df / (df + ts))
        p = np.zeros(shape=r.shape)
        p[np.triu_indices(p.shape[0], 1)] = pf
        p[np.tril_indices(p.shape[0], -1)] = pf
        p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])

        r = np.nan_to_num(r)
        p = np.nan_to_num(p)

        # make pvalues triangular
        for i in range(len(matrix)):
            for j in range(i, len(matrix)):
                p[j][i] = p[i][j]

        return r, p

    # It is used to find the 10th percentile value without considering zeros
    def find_first_perc(self, A):
        for i in range(0, len(A)):
            if A[i] != 0:
                first_indices = i
                one_perc = int((10 * (len(A) - first_indices)) / 100)
                return A[i + one_perc], A[i:i + one_perc]

    #Threshold the euclidean matrix
    def threshold_eucl(self, X):

        p99 = np.percentile(X,99)
        X[X<p99] = 0

        return X

    #Normalize it to have values from 0 to 1
    def normalize_eucl(self):

        N = self.threshold_eucl(self.normal_eucl)
        C = self.threshold_eucl(self.cancer_eucl)

        N_log = np.log10(N)
        N_log_norm = N_log/N_log.max()
        N_log_norm[ N_log_norm == -np.inf] = 0

        C_log = np.log10(C)
        C_log_norm = C_log/C_log.max()
        C_log_norm[C_log_norm == -np.inf] = 0

        np.save('./Networks/normal_eucl_'+ str(self.tumor)+'.npy', N_log_norm)
        np.save('./Networks/cancer_eucl_'+ str(self.tumor)+'.npy', C_log_norm)

        print('Eucl ' + str(self.tumor)+ ' normalized')

        return N_log_norm, C_log_norm

    # It allows to permutate the Pearson matrices ten times
    def permutation(self, N=10):
        if N < 2:
            raise ValueError('N too small')
        normal_pearson_p_all = np.zeros((len(self.normal), len(self.normal)))
        cancer_pearson_p_all = np.zeros((len(self.cancer), len(self.cancer)))

        for _ in range(N):
            normal_p = pd.DataFrame(self.permutate(self.normal.values))
            normal_p.columns = self.normal.columns
            cancer_p = pd.DataFrame(self.permutate(self.cancer.values))
            cancer_p.columns = self.cancer.columns
            index_mirna_null_normal_p = normal_p[normal_p.duplicated()].index
            index_mirna_null_cancer_p = cancer_p[cancer_p.duplicated()].index
            normal_pearson_p = np.corrcoef(normal_p)
            cancer_pearson_p = np.corrcoef(cancer_p)
            normal_pearson_p = self.fill_null_mirna(index_mirna_null_normal_p, normal_pearson_p)
            cancer_pearson_p = self.fill_null_mirna(index_mirna_null_cancer_p, cancer_pearson_p)
            normal_norm_p = pd.DataFrame(self.permutate(self.normal_normalized.values))
            normal_norm_p.columns = self.normal.columns
            cancer_norm_p = pd.DataFrame(self.permutate(self.cancer_normalized.values))
            cancer_norm_p.columns = self.cancer.columns
            normal_pearson_p_all = normal_pearson_p_all + normal_pearson_p
            cancer_pearson_p_all = cancer_pearson_p_all + cancer_pearson_p

        normal_pearson_p_all = normal_pearson_p_all / 10
        cancer_pearson_p_all = cancer_pearson_p_all / 10
        print('permutation calculated')
        return normal_pearson_p_all, cancer_pearson_p_all

    #Threshold Pearson networks using the permutation
    def threshold_permuted(self, matrix_original, matrix_permutated):
        
        linspace=np.linspace(matrix_original.min(), matrix_original.max(), num=20)
        count=np.zeros(20)
        #Low and high random values
        val_low=-99999999999999999
        val_high=99999999999999999
        for i in range(20):
            count[i]=np.count_nonzero(matrix_permutated[(matrix_permutated<(linspace[i]))])/2
            if (i!=0):
                count[i]=count[i]-count[i-1]
            if (count[i]==0):
                if val_low<linspace[i]:
                    val_low=linspace[i]
            if (count[i]!=0):
                if (i<19):
                    if val_high>linspace[i+1]:
                        val_high=linspace[i+1]
       
        matrix_original[((matrix_original>(val_low)) & (matrix_original>(val_high)))]=0
        return matrix_original

    # Threshold Pearson networks using the p-values
    # We pass two arrays of percentiles for threshold on pvalues such that the Pearson networks have a number of edges similar to the Euclidean one
    def threshold_pvalue(self, normal_perc, cancer_perc):
        val_eucl=(self.normal_eucl_thr!=0).sum()/2
        normal_pearson, pvalues_normal=self.corrcoef(self.normal)
        cancer_pearson, pvalues_cancer=self.corrcoef(self.cancer)
        pvalues_normal_df=pd.DataFrame(pvalues_normal)
        pvalues_cancer_df=pd.DataFrame(pvalues_cancer)
        pvalues_normal=np.zeros((len(self.normal),len(self.normal)))
        pvalues_cancer=np.zeros((len(self.normal),len(self.normal)))
        
        pvalues_normal[self.normal_pearson_thr!=0]=pvalues_normal_df.values[self.normal_pearson_thr!=0]
        pvalues_cancer[self.cancer_pearson_thr!=0]=pvalues_cancer_df.values[self.cancer_pearson_thr!=0]
        pvalues_normal=np.nan_to_num(pvalues_normal)
        pvalues_cancer=np.nan_to_num(pvalues_cancer)
        pvalues_normal_df=pd.DataFrame(pvalues_normal)
        pvalues_cancer_df=pd.DataFrame(pvalues_cancer)
        triang_sup_noDiag_normal=list(pvalues_normal_df.get_values()[np.triu_indices(len(self.normal), k=1)])
        triang_sup_noDiag_normal_df=pd.DataFrame(triang_sup_noDiag_normal)
        pvalues_normal_array=triang_sup_noDiag_normal_df.get_values().flatten()
        triang_sup_noDiag_cancer=list(pvalues_cancer_df.get_values()[np.triu_indices(len(self.cancer), k=1)])
        triang_sup_noDiag_cancer_df=pd.DataFrame(triang_sup_noDiag_cancer)
        pvalues_cancer_array=triang_sup_noDiag_cancer_df.get_values().flatten()
        critical_pval_normal, ten_perc_pvalues_normal=self.find_first_perc(pvalues_normal_array)
        critical_pval_cancer, ten_perc_pvalues_cancer=self.find_first_perc(pvalues_cancer_array)

        for perc in normal_perc:
            percentile_normal=np.percentile(ten_perc_pvalues_normal, perc)
            n_p = self.normal_pearson_thr.copy()
            n_p[pvalues_normal_df.values>=percentile_normal]=0
            n_p = np.nan_to_num(n_p)
            sum_n =(n_p!=0).sum()
            print('pval: ' +str(perc))
            print('normal_pears: '+str(sum_n))
         
            if (sum_n<5000000):
                if(abs(sum_n/2 - val_eucl)<1000000):
                    print('Here!Perc: '+str(perc))
                    break
                elif(sum_n/2 - val_eucl)<0:
                    for perc1 in range(perc+1, perc+4):  
                        percentile_normal=np.percentile(ten_perc_pvalues_normal, perc1)
                        n_p = self.normal_pearson_thr
                        n_p[pvalues_normal_df.values>=percentile_normal]=0
                        sum_n =(n_p!=0).sum()
                        if(abs(sum_n/2 - val_eucl)<1000000):
                            break
        print('thresh pval: ' +str(perc))
        print('normal_pears: '+str((n_p!=0).sum()))

        for perc in cancer_perc:
            
            percentile_cancer=np.percentile(ten_perc_pvalues_cancer, perc)
            c_p = self.cancer_pearson_thr
            c_p[pvalues_cancer_df.values>=percentile_cancer]=0
            sum_c =(c_p!=0).sum()
        
            if (sum_c<5000000):
                if(abs(sum_c/2 - val_eucl)<1000000):
                    print('Here!Perc: '+str(perc))
                    break
                elif(sum_c/2 - val_eucl)<0:
                    for perc1 in range(perc+1, perc+4):  
                        percentile_cancer=np.percentile(ten_perc_pvalues_cancer, perc1)
                        c_p = self.cancer_pearson_thr
                        c_p[pvalues_cancer_df.values>=percentile_cancer]=0
                        sum_c =(c_p!=0).sum()
                        if(abs(sum_c/2 - val_eucl)<1000000):
                            break
                        
        print('thresh pval: ' +str(perc))
        print('cancer_pears: '+str((c_p!=0).sum()))
        
        np.save('./Networks/normal_pearson_'+ str(self.tumor)+'.npy', n_p)
        np.save('./Networks/cancer_pearson_'+ str(self.tumor)+'.npy', c_p)
        return n_p, c_p
   
    
    #Put together the two measures
    def sum_pears_eucl(self):

        n_e = self.normal_eucl_thr
        c_e = self.cancer_eucl_thr
        

        n_abs = abs(self.normal_pearson_thr)
        c_abs = abs(self.cancer_pearson_thr)

        n = n_e + n_abs
        c = c_e + c_abs

        np.save('./Networks/normal_sum_pears_eucl_' + str(self.tumor) +'.npy', n)
        np.save('./Networks/cancer_sum_pears_eucl_' + str(self.tumor) +'.npy', c)

        print('Normal and cancer ' + str(self.tumor)+ ' done')

        return n,c   
