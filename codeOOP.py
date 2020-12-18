#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import umap.umap_ as umap
import glob
from sklearn import preprocessing
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import utils
from collections import Counter

class DataProcesss:
#initializing the program, importing sampleNeeded (sample no. to patient info)
    def __init__(self, data1, path2):
        self.sampleNeeded = pd.read_csv(data1, sep='\t')
        self.sampleNeeded = self.sampleNeeded.loc[self.sampleNeeded['Characteristics[treatment]'].isin(['none','placebo_none','placebo_placebo'])]
        self.sampleNeeded = self.sampleNeeded[self.sampleNeeded['Characteristics[organism part]'] == 'blood']

        all_files = glob.glob(path2 + "/*.txt")

        allFiles = []
        wantedList = self.sampleNeeded['Derived Array Data File'].tolist()
        fileNameList=[]
#go through each patient microarray data and compile dataframe
        for filename in all_files:
            if filename[-27:] in wantedList:
                fileNameList.append(filename[-27:])
                df = pd.read_csv(filename, sep='\t')
                allFiles.append(df)

        masterlist = pd.DataFrame(columns = allFiles[0]['ID_REF'].tolist())
        for i in range(len(allFiles)):
            masterlist.loc[i]= allFiles[i]['VALUE'].tolist()

        self.sampleNeeded = self.sampleNeeded.rename(columns = {'Derived Array Data File':'ID'})

        masterlist=masterlist.values
#scale the expression data with minmax method
        min_max_scaler = preprocessing.MinMaxScaler()
        self.masterlistScaled_before= min_max_scaler.fit_transform(masterlist)
        self.masterlistScaled_before = pd.DataFrame(self.masterlistScaled_before)
        self.masterlistScaled = self.masterlistScaled_before.copy()
        self.masterlistScaled.insert(loc = 0, column = "ID", value = fileNameList)

#only looking at diseased pop data
        finalList = self.masterlistScaled.merge(self.sampleNeeded, on = 'ID', how = 'inner')
        self.finalListNoDup = finalList.drop_duplicates(subset=['ID'])
        self.DiseasedPop = self.finalListNoDup.copy()

        self.DiseasedPop  = self.DiseasedPop.set_index(self.masterlistScaled_before.index)
        normalIndex = self.DiseasedPop[self.DiseasedPop["Factor Value[disease]"] == 'normal'].index

        self.DiseasedPop = self.DiseasedPop.drop(normalIndex)

        self.masterlistScaledDiseased = self.masterlistScaled_before.copy()
        self.masterlistScaledDiseased = self.masterlistScaledDiseased.drop(normalIndex)

#umap for dimensionaity reduction
    def umap(self,n_neigh,min_d,n_comp):
        self.reducerDiseased = umap.UMAP(n_neighbors=15, min_dist=0.15, n_components=2).fit(self.masterlistScaledDiseased)
        self.reducedDataDiseased = pd.DataFrame(self.reducerDiseased.transform(self.masterlistScaledDiseased))
        return self.reducedDataDiseased

#using DBSCAN to cluster data based on density
    def densCluster(self,EPS,min_samp):
        db = DBSCAN(eps=EPS, min_samples=min_samp).fit(self.reducedDataDiseased)
        self.labelsDiseased=db.labels_
        self.finalTestDiseased =self.DiseasedPop.copy()
        self.finalTestDiseased = self.finalTestDiseased.set_index(self.reducedDataDiseased.index)


        self.finalTestDiseased.insert(loc = 1, column = "UMAPx", value = self.reducedDataDiseased[0])
        self.finalTestDiseased.insert(loc = 1, column = "UMAPy", value = self.reducedDataDiseased[1])
        self.finalTestDiseased.insert(loc = 1, column = "clusterID", value = self.labelsDiseased)

#The dataset is too large and took very long to run. The code is what I used and I saved it for users to directly import and use
    # def rankExpression(self,data):
    #     list_string = list(map(str, list(data.columns)))
    #     self.masterlistScaledDiseasedRank = pd.DataFrame(columns = list_string)
    #
    #     for i in range(len(list_string)):
    #        self.masterlistScaledDiseasedRank[list_string[i]] =self. masterlistScaledDiseased[i].rank(pct=True)
    #
    #     self.masterlistScaledDiseasedRank_transpose = masterlistScaledDiseased.T
    #
    #     masterlistScaledDiseasedRank_transpose.to_csv('rankedData.csv', index=True)

#ranking the expressio of each gene from different patient using
    def rankedDataProcessing(self,path):
        self.rankedData = pd.read_csv(path)
        self.rankedData.insert(loc = 0, column = "clusterID", value = self.finalTestDiseased['clusterID'].tolist())
        clustercopy = self.rankedData.copy()
        n_clusters_ = len(set(self.labelsDiseased)) - (1 if -1 in self.labelsDiseased else 0)

        listofClusters=[]
        for i in range(n_clusters_):
            clusterTemp = clustercopy[clustercopy['clusterID'] == i]
            listofClusters.append(clusterTemp)

        for i in range(len(listofClusters)):
            listofClusters[i] = listofClusters[i].drop(columns='clusterID')
            listofClusters[i] = listofClusters[i].T
        listofDysGenebyCluster=[]

        for i in range(len(listofClusters)):
            listofDysGeneinOneCluster=[]
            for j in range(len(listofClusters[i].index)):
                if all((k<0.05 for k in listofClusters[0].iloc[j,:])) or all((k>0.95 for k in listofClusters[i].iloc[j,:])):
                    listofDysGeneinOneCluster.append(j)
            listofDysGenebyCluster.append(listofDysGeneinOneCluster)

        totalgene=[]
        for i in range(len(listofDysGenebyCluster)):
            totalgene =  totalgene+listofDysGenebyCluster[i]

        #30 genes upregulated in total of all clusters
        totalgeneUnique = list(set(totalgene))
        clusterNoandUniqueGene=[]

        for item, count in Counter(totalgene).most_common():
            for i in range(10):
                if count == i:
                    clusterNoandUniqueGene.append([i,item])
        # len(clusterNoandUniqueGene)
        clusterNoandUniqueGene
        clusterNoandUniqueGene_pd=pd.DataFrame(clusterNoandUniqueGene,columns=['count','index'])

        codetoGene = pd.read_csv('codeToGene.txt', sep='\t')

        matchingGene=[]
        for i in range(len(clusterNoandUniqueGene_pd.index)):
            matchingGene.append(codetoGene.iloc[clusterNoandUniqueGene_pd.iloc[i,1],1])
        matchingGene
        clusterNoandUniqueGene_pd.insert(loc = 2, column = "gene", value = matchingGene)
        clusterGeneMasterlist=[]

        for geneIndex in clusterNoandUniqueGene_pd['index']:
            print(geneIndex)
            clusterswithGene=[]
            for i in range(len(listofDysGenebyCluster)):
                if geneIndex in listofDysGenebyCluster[i]:
                    clusterswithGene.append(i)
            clusterGeneMasterlist.append(clusterswithGene)

        clusterNoandUniqueGene_pd.insert(loc = 3, column = "clusterWithGene", value = clusterGeneMasterlist)

        return clusterNoandUniqueGene_pd


if __name__ == '__main__':
    data1='sampleInfo.txt'
    path2=r"C:\\Users\\jeffy\\OneDrive\\Desktop\\IBDdataset\\combined\\"
    result=DataProcesss(data1,path2)
    result.umap(15, 0.15, 2)
    result.densCluster(0.5,8)
    final=result.rankedDataProcessing('rankedDataNoIndexNoTransp.csv')
    print(final)
