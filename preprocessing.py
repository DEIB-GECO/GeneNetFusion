import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

class Preprocessing:
    def __init__(self,miRNACancer, miRNACancerlength, miRNANormal, miRNANormallength, geneCancer, geneNormal, tumor):
        self.miRNACancer=miRNACancer
        self.miRNACancerlength=miRNACancerlength
        self.miRNANormal=miRNANormal
        self.miRNANormallength=miRNANormallength
        self.geneCancer=geneCancer
        self.geneNormal=geneNormal
        self.tumor=tumor
        self.fix_matrices()
        self.geneCancer, self.geneNormal= self.filter_genes()
        self.filter_noise()
        self.Normal, self.Cancer, self.Indices=self.matrices()
        
        
    #It allows to have matrices that contains fpkm with ensembl as index and patients in columns
    def fix_matrices(self):
        #Remove indices not necessary
        self.miRNACancer.index=self.miRNACancer.index.droplevel(0)
        self.miRNANormal.index=self.miRNANormal.index.droplevel(0)
        self.miRNANormallength.index=self.miRNANormallength.index.droplevel(0)
        self.miRNACancerlength.index=self.miRNACancerlength.index.droplevel(0)
        self.miRNANormal.columns=self.miRNANormal.columns.droplevel(0)
        self.miRNACancer.columns=self.miRNACancer.columns.droplevel(0)
        self.miRNANormal.columns.name=None
        self.miRNACancer.columns.name=None

        self.geneNormal.columns=self.geneNormal.columns.droplevel(0)
        self.geneCancer.columns=self.geneCancer.columns.droplevel(0)
        self.geneNormal.columns.name=None
        self.geneCancer.columns.name=None
        print('geneNorm ' + str(self.geneNormal.shape))
        print('geneCanc ' + str(self.geneCancer.shape))
        l=len(self.miRNANormal)

        #Transform rpm into fpkm
        for i in range (0, l-1):
            self.miRNANormal.get_values()[i]=(self.miRNANormal.get_values()[i]*(10**3))/(2*abs(self.miRNANormallength.get_values()[i,1]-self.miRNANormallength.get_values()[i,0]))
            self.miRNACancer.get_values()[i]=(self.miRNACancer.get_values()[i]*(10**3))/(2*abs(self.miRNACancerlength.get_values()[i,1]-self.miRNACancerlength.get_values()[i,0]))
            
        self.miRNANormal.index.name='ensembl'
        self.miRNACancer.index.name='ensembl'

    #Remove genes that are not protein coding, long non coding or mirna.
    #Remove also duplicated genes
    def filter_genes(self):
        matrix_GTypeGName= pd.read_excel('matrix_GTypeGName.xls')
        symbolAndEnsembl=pd.DataFrame(self.geneNormal.index.get_level_values(0), self.geneNormal.index.get_level_values(1))
        self.geneCancer['ensembl']=self.geneCancer.index.get_level_values(0)
        self.geneNormal['ensembl']=self.geneNormal.index.get_level_values(0)
        GeneCancerTypes=pd.merge(self.geneCancer,matrix_GTypeGName, left_on="gene_symbol", right_on="Gene name")
        GeneNormalTypes=pd.merge(self.geneNormal,matrix_GTypeGName, left_on="gene_symbol", right_on="Gene name")
        GeneCancerTypes.index=GeneCancerTypes['ensembl']
        GeneNormalTypes.index=GeneNormalTypes['ensembl']
        GeneCancerTypes=GeneCancerTypes.drop_duplicates()
        GeneNormalTypes=GeneNormalTypes.drop_duplicates()
        
        #Remove all gene types that are not interesting
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_J_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='polymorphic_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='unitary_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_D_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='scaRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_V_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='sense_overlapping')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='rRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_J_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_C_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='bidirectional_promoter_lncRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_C_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='3prime_overlapping_ncRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_C_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_J_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='ribozyme')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_J_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_D_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='Mt_rRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='antisense')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='misc_RNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='rRNA_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='transcribed_unprocessed_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='processed_transcript')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='transcribed_processed_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_V_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='sense_intronic')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='IG_V_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TR_V_gene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='TEC')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='transcribed_unitary_pseudogene')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='vaultRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='macro_lncRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='scRNA')]
        GeneCancerTypes=GeneCancerTypes[(GeneCancerTypes['Gene type']!='non_coding')]

        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_J_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='polymorphic_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='unitary_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_D_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='scaRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_V_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='sense_overlapping')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='rRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_J_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_C_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='bidirectional_promoter_lncRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_C_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='3prime_overlapping_ncRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_C_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_J_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='ribozyme')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_J_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_D_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='Mt_rRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='antisense')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='misc_RNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='rRNA_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='transcribed_unprocessed_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='processed_transcript')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='transcribed_processed_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_V_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='sense_intronic')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='IG_V_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TR_V_gene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='TEC')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='transcribed_unitary_pseudogene')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='vaultRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='macro_lncRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='scRNA')]
        GeneNormalTypes=GeneNormalTypes[(GeneNormalTypes['Gene type']!='non_coding')]

        
        #Remove the shortest gene since data are random
        GeneCancer_noMiRNA=GeneCancerTypes[(GeneCancerTypes['Gene type']!='miRNA')]
        GeneNormal_noMiRNA=GeneNormalTypes[(GeneNormalTypes['Gene type']!='miRNA')]
        GeneCancer_noSnRNA=GeneCancer_noMiRNA[(GeneCancer_noMiRNA['Gene type']!='snRNA')]
        GeneNormal_noSnRNA=GeneNormal_noMiRNA[(GeneNormal_noMiRNA['Gene type']!='snRNA')]
        GeneCancer_noSnoRNA=GeneCancer_noSnRNA[(GeneCancer_noSnRNA['Gene type']!='snoRNA')]
        GeneNormal_noSnoRNA=GeneNormal_noSnRNA[(GeneNormal_noSnRNA['Gene type']!='snoRNA')]
        GeneCancer_noTRNA=GeneCancer_noSnoRNA[(GeneCancer_noSnoRNA['Gene type']!='Mt_tRNA')]
        GeneNormal_noTRNA=GeneNormal_noSnoRNA[(GeneNormal_noSnoRNA['Gene type']!='Mt_tRNA')]


        #Drop the type, keeping only fpkm and  as index gene symbol
        GeneCancerLong=GeneCancer_noTRNA.drop(['Gene type','Gene name', 'ensembl'], axis=1)
        GeneNormalLong=GeneNormal_noTRNA.drop(['Gene type','Gene name', 'ensembl'], axis=1)

        
        GeneCancerLong.index.name='ensembl'
        GeneNormalLong.index.name='ensembl'
        
        indexDuplicated_cancer=GeneCancerLong[GeneCancerLong.index.duplicated()].index
        indexDuplicated_normal=GeneNormalLong[GeneNormalLong.index.duplicated()].index
        indexDuplicated=indexDuplicated_cancer.intersection(indexDuplicated_normal)
        noDuplCancer=GeneCancerLong.drop(indexDuplicated)
        noDuplNormal=GeneNormalLong.drop(indexDuplicated)
        DuplCancer=GeneCancerLong.drop(noDuplCancer.index)
        DuplNormal=GeneNormalLong.drop(noDuplNormal.index)
        geneCancer_noDupl=noDuplCancer.append(DuplCancer.drop_duplicates())
        geneNormal_noDupl=noDuplNormal.append(DuplNormal.drop_duplicates())
        diff=geneCancer_noDupl.index.difference(geneNormal_noDupl.index)
        geneCancer_noDupl=geneCancer_noDupl.drop(diff, axis=0 )
        diff_norm=geneNormal_noDupl.index.difference(geneCancer_noDupl.index)
        geneNormal_noDupl=geneNormal_noDupl.drop(diff_norm, axis=0 )
        print('normal ' + str(geneNormal_noDupl.shape))
        print('cancer ' + str(geneCancer_noDupl.shape))
        return geneCancer_noDupl, geneNormal_noDupl
    
    #Remove noisy genes using z-score of fpkm
    def filter_noise(self):
        #Calculate the mean of each gene
        miRNANormal_mean=self.miRNANormal.mean(axis=1)
        miRNACancer_mean=self.miRNACancer.mean(axis=1)
        geneNormal_mean=self.geneNormal.mean(axis=1)
        geneCancer_mean=self.geneCancer.mean(axis=1)
        
        #Extract indices of all the genes/miRNA that have mean equal to 0
        miRNANormalIndex_null=miRNANormal_mean[miRNANormal_mean.get_values()==0].index
        miRNACancerIndex_null=miRNACancer_mean[miRNACancer_mean.get_values()==0].index
        geneNormalIndex_null=geneNormal_mean[geneNormal_mean.get_values()==0].index
        geneCancerIndex_null=geneCancer_mean[geneCancer_mean.get_values()==0].index
        #Find only those that are in both cancer and normal equal to 0
        miRNA_null=miRNANormalIndex_null.intersection(miRNACancerIndex_null)
        gene_null=geneNormalIndex_null.intersection(geneCancerIndex_null)
        #Drop them
        self.miRNANormal=self.miRNANormal.drop(miRNA_null, axis=0)
        self.miRNACancer=self.miRNACancer.drop(miRNA_null, axis=0)
        self.geneNormal=self.geneNormal.drop(gene_null, axis=0)
        self.geneCancer=self.geneCancer.drop(gene_null, axis=0)

        #Recalculate the mean
        geneNormal_mean=self.geneNormal.mean(axis=1)
        geneCancer_mean=self.geneCancer.mean(axis=1)

        #Calculate log2
        geneNormal_log2=np.log2(self.geneNormal)
        geneCancer_log2=np.log2(self.geneCancer)

        #Calculate the standard deviation
        geneNormal_std=self.geneNormal.std(axis=1)
        geneCancer_std=self.geneCancer.std(axis=1)

        #Calculate zfpkm
        geneNormal_zfpkm=(geneNormal_log2.sub(geneNormal_mean, axis=0)).divide(geneNormal_std, axis=0)
        geneCancer_zfpkm=(geneCancer_log2.sub(geneCancer_mean, axis=0)).divide(geneCancer_std, axis=0)

        #Calculate the mean for each gene
        geneNormal_mean_zfpkm=geneNormal_zfpkm.mean(axis=1)
        geneCancer_mean_zfpkm=geneCancer_zfpkm.mean(axis=1)
    
        geneNormalIndex_less=geneNormal_mean_zfpkm[geneNormal_mean_zfpkm.get_values()<(-3)].index
        geneCancerIndex_less=geneCancer_mean_zfpkm[geneCancer_mean_zfpkm.get_values()<(-3)].index
        gene_less=geneNormalIndex_less.intersection(geneCancerIndex_less)

        #Drop them
        self.geneNormal=self.geneNormal.drop(gene_less, axis=0)
        self.geneCancer=self.geneCancer.drop(gene_less, axis=0)
        
        
    #Retrieve common patients between genes and mirnas
    def common_patients(self):
        #Find common patients for cancer and common for normal
        a=pd.DataFrame(self.miRNACancer.columns)
        b=pd.DataFrame(self.geneCancer.columns)
        Cancerpatients=pd.merge(a,b)

        c=pd.DataFrame(self.miRNANormal.columns)
        d=pd.DataFrame(self.geneNormal.columns)
        Normalpatients=pd.merge(c,d)

        #Remove patients differents
        MiRNACancerDiff=self.miRNACancer.columns.difference(Cancerpatients[0])
        MiRNACancerCommon=self.miRNACancer.columns.difference(self.miRNACancer[MiRNACancerDiff].columns)
        self.miRNACancer=self.miRNACancer[MiRNACancerCommon]
        MiRNANormalDiff=self.miRNANormal.columns.difference(Normalpatients[0])
        MiRNANormalCommon=self.miRNANormal.columns.difference(self.miRNANormal[MiRNANormalDiff].columns)
        self.miRNANormal=self.miRNANormal[MiRNANormalCommon]
        GeneCancerDiff=self.geneCancer.columns.difference(Cancerpatients[0])
        GeneCancerCommon=self.geneCancer.columns.difference(self.geneCancer[GeneCancerDiff].columns)
        self.geneCancer=self.geneCancer[GeneCancerCommon]
        GeneNormalDiff=self.geneNormal.columns.difference(Normalpatients[0])
        GeneNormalCommon=self.geneNormal.columns.difference(self.geneNormal[GeneNormalDiff].columns)
        self.geneNormal=self.geneNormal[GeneNormalCommon]

    #It create the two matrices and it saves the indices with the respective ensembl
    def matrices(self):
        #Create the two matrices
        Normal=self.geneNormal.append(self.miRNANormal)
        Cancer=self.miRNACancer.append(self.geneCancer)


        #Order the index
        Normal=Normal.sort_index()
        Cancer=Cancer.sort_index()
        #Cancer_stage=Cancer_stage.sort_index()

        #Fill nan with 0
        Normal=Normal.fillna(0)
        Cancer=Cancer.fillna(0)
        
        Normal.to_csv('./Matrices/normal_'+str(self.tumor),sep=';')
        Cancer.to_csv('./Matrices/cancer_'+str(self.tumor),sep=';')
        
        #Save indices
        Indices=pd.DataFrame(np.arange(0,len(Normal)), Normal.index)
        Indices.to_csv('./Matrices/indices_'+str(self.tumor), sep=';')

        return Normal, Cancer, Indices


