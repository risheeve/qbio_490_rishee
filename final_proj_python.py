#!/usr/bin/env python
# coding: utf-8

# In[18]:


import os

os.chdir('/PATH/TO/final_project')

# import necessary libraries
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats

import umap
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
import warnings
import mygene
warnings.filterwarnings('ignore')

# In[17]:


# Import cptac
import cptac
# Examine the data sets available with list_datasets()
cptac.list_datasets()


# In[19]:


# Download the colon cancer data set
cptac.download(dataset="Colon")
colon = cptac.Colon()


# In[5]:


colon.list_data()


# In[20]:


protein_data = colon.get_proteomics() # get the proteomics data
protein_data = protein_data.drop_duplicates()
protein_data.columns = protein_data.columns.get_level_values(0) 
protein_data.index = protein_data.index.get_level_values(0) 

rna_data = colon.get_transcriptomics() # get the transcriptomics data
rna_data = rna_data.drop_duplicates()
rna_data.columns = rna_data.columns.get_level_values(0) 
rna_data.index = rna_data.index.get_level_values(0) 

clinical_data = colon.get_clinical() # get the clinical data
mutation_data = colon.get_somatic_mutation() # get the mutation data


# In[23]:


# load data from csv
protein_data = pd.read_csv('CPTAC_COAD_Proteomics.csv', index_col=0) # get the proteomics data
protein_data = protein_data.drop_duplicates()
protein_data.columns = protein_data.columns.get_level_values(0) 
protein_data.index = protein_data.index.get_level_values(0) 

rna_data = pd.read_csv('CPTAC_COAD_RNA.csv', index_col=0) # get the transcriptomics data
rna_data = rna_data.drop_duplicates()
rna_data.columns = rna_data.columns.get_level_values(0) 
rna_data.index = rna_data.index.get_level_values(0) 

clinical_data = pd.read_csv('CPTAC_COAD_Clinical.csv', index_col=0) # get the clinical data
mutation_data = pd.read_csv('CPTAC_COAD_Mutation.csv', index_col=0) # get the mutation data


# In[24]:


# remove na values ()
rna_data = rna_data.loc[:, rna_data.count(axis=0) >= 3]
protein_data = protein_data.loc[:, protein_data.count(axis=0) >= 3]


# In[25]:


# adding columns with helpful information
clinical_data["Age_in_year"] = clinical_data["Age"]/12
clinical_data["diagnos_status"] = np.where(clinical_data["Age_in_year"]<=50, "Early Onset", "Late Onset")


# In[26]:


# get data from patients with tumor
tumor_mask = clinical_data["Sample_Tumor_Normal"]=="Tumor"
clinical_data = clinical_data[tumor_mask]


# In[27]:


early_mask = np.where(clinical_data["diagnos_status"]=="Early Onset")
early_data = clinical_data.iloc[early_mask]
late_mask = np.where(clinical_data["diagnos_status"]=="Late Onset")
late_data = clinical_data.iloc[late_mask]


# In[28]:


early_rna_mask = np.intersect1d(early_data.index, rna_data.index)
early_rna = rna_data.loc[early_rna_mask, :]
early_prot_mask = np.intersect1d(early_data.index, protein_data.index)
early_prot = protein_data.loc[early_prot_mask, :]
# make sure patients are share in rna and protein
early_rna_prot = np.intersect1d(early_rna.index, early_prot.index)
early_rna = rna_data.loc[early_rna_prot, :]
early_prot = protein_data.loc[early_rna_prot, :]

late_rna_mask = np.intersect1d(late_data.index, rna_data.index)
late_rna = rna_data.loc[late_rna_mask, :]
late_prot_mask = np.intersect1d(late_data.index, protein_data.index)
late_prot = protein_data.loc[late_prot_mask, :]

# make sure patients are share in rna and protein
late_rna_prot = np.intersect1d(late_rna.index, late_prot.index)
late_rna = rna_data.loc[late_rna_prot, :]
late_prot = protein_data.loc[late_rna_prot, :]


# In[29]:


# select only genes share in rna and protein datasets
shared_rna_prot = np.intersect1d(rna_data.columns, protein_data.columns)

early_rna_shared = early_rna.loc[:, shared_rna_prot]
early_prot_shared = early_prot.loc[:, shared_rna_prot]
late_rna_shared = late_rna.loc[:, shared_rna_prot]
late_prot_shared = late_prot.loc[:, shared_rna_prot]


# In[30]:


# find the most mutated genes with available in both omics
share_mutation = mutation_data.loc[mutation_data['Gene'].isin(shared_rna_prot)] 
mutation_counts = share_mutation['Gene'].value_counts()
high_mutation = mutation_counts.index.tolist()


# In[81]:


early_rna_shared['MXRA8']


# In[82]:


early_prot_shared['MXRA8']


# In[61]:


gene_names


# In[76]:


# becuase of problems with CPTAC, we will not be using this graph
# set top20 mutated genes as the target
gene_names = high_mutation[:16]
ncomparisons = 16

# create heatmap for female early onset
fig, ax = plt.subplots()
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        corr, pval = stats.spearmanr(early_rna_shared[g1], early_prot_shared[g2], nan_policy = 'omit')
        corr_df.loc[g1, g2] = corr
plot = sns.heatmap(
    corr_df,
    cmap='mako',
    ax = ax
)
fig.suptitle('Early Onset')
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig('early.png', bbox_inches='tight')


# In[11]:


# create separate dataset for both genders
female_mask = np.where(clinical_data["Gender"]=="Female")
female_clinical_data = clinical_data.iloc[female_mask]
male_mask = np.where(clinical_data["Gender"]=="Male")
male_clinical_data = clinical_data.iloc[male_mask]


# In[12]:


# create separate dataset for early and late onset
female_early_mask = np.where(female_clinical_data["diagnos_status"]=="Early Onset")
female_early_data = clinical_data.iloc[female_early_mask]
female_late_mask = np.where(female_clinical_data["diagnos_status"]=="Late Onset")
female_late_data = clinical_data.iloc[female_late_mask]


# In[13]:


# create separate dataset for early and late onset
male_early_mask = np.where(male_clinical_data["diagnos_status"]=="Early Onset")
male_early_data = clinical_data.iloc[male_early_mask]
male_late_mask = np.where(male_clinical_data["diagnos_status"]=="Late Onset")
male_late_data = clinical_data.iloc[male_late_mask]


# In[14]:


# create separate rna and protein dataset for female early, female late, male early, and male late
fe_rna_mask = np.intersect1d(female_early_data.index, rna_data.index)
fe_rna = rna_data.loc[fe_rna_mask, :]
fe_prot_mask = np.intersect1d(female_early_data.index, protein_data.index)
fe_prot = protein_data.loc[fe_prot_mask, :]
# patients in fe_rna and fe_prot overlap completely, so no adjustments needed

fl_rna_mask = np.intersect1d(female_late_data.index, rna_data.index)
fl_rna = rna_data.loc[fl_rna_mask, :]
fl_prot_mask = np.intersect1d(female_late_data.index, protein_data.index)
fl_prot = protein_data.loc[fl_prot_mask, :]
# make sure patients are share in rna and protein
fl_rna_prot = np.intersect1d(fl_rna.index, fl_prot.index)
fl_rna = rna_data.loc[fl_rna_prot, :]
fl_prot = protein_data.loc[fl_rna_prot, :]

me_rna_mask = np.intersect1d(male_early_data.index, rna_data.index)
me_rna = rna_data.loc[me_rna_mask, :]
me_prot_mask = np.intersect1d(male_early_data.index, protein_data.index)
me_prot = protein_data.loc[me_prot_mask, :]
# patients in me_rna and me_prot overlap completely, so no adjustments needed

ml_rna_mask = np.intersect1d(male_late_data.index, rna_data.index)
ml_rna = rna_data.loc[ml_rna_mask, :]
ml_prot_mask = np.intersect1d(male_late_data.index, protein_data.index)
ml_prot = protein_data.loc[ml_prot_mask, :]
# make sure patients are share in rna and protein
ml_rna_prot = np.intersect1d(ml_rna.index, ml_prot.index)
ml_rna = rna_data.loc[ml_rna_prot, :]
ml_prot = protein_data.loc[ml_rna_prot, :]


# In[15]:


# select only genes share in rna and protein datasets
shared_rna_prot = np.intersect1d(rna_data.columns, protein_data.columns)

fe_rna_shared = fe_rna.loc[:, shared_rna_prot]
fe_prot_shared = fe_prot.loc[:, shared_rna_prot]
fl_rna_shared = fl_rna.loc[:, shared_rna_prot]
fl_prot_shared = fl_prot.loc[:, shared_rna_prot]

me_rna_shared = me_rna.loc[:, shared_rna_prot]
me_prot_shared = me_prot.loc[:, shared_rna_prot]
ml_rna_shared = ml_rna.loc[:, shared_rna_prot]
ml_prot_shared = ml_prot.loc[:, shared_rna_prot]


# In[16]:


# find the most mutated genes with available in both omics
share_mutation = mutation_data.loc[mutation_data['Gene'].isin(shared_rna_prot)] 
mutation_counts = share_mutation['Gene'].value_counts()
high_mutation = mutation_counts.index.tolist()


# In[21]:


ml_prot_shared


# In[19]:


# set top20 mutated genes as the target
gene_names = high_mutation[:20]
ncomparisons = 20

# create heatmap for female early onset
fig, ax = plt.subplots()
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        corr, pval = stats.spearmanr(fe_rna_shared[g1], fe_prot_shared[g2], nan_policy = 'omit')
        corr_df.loc[g1, g2] = corr
plot = sns.heatmap(
    corr_df,
    cmap='mako',
    ax = ax
)
fig.suptitle('Female Early Onset')
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig('female_early.png', bbox_inches='tight')


# In[20]:


gene_names = high_mutation[:20]
ncomparisons = 20

# create heatmap for female late onset
fig, ax = plt.subplots()
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        corr, pval = stats.spearmanr(fl_rna_shared[g1], fl_prot_shared[g2], nan_policy = 'omit')
        corr_df.loc[g1, g2] = corr
plot = sns.heatmap(
    corr_df,
    cmap='mako',
    ax = ax
)
fig.suptitle('Female Late Onset')
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig('female_late.png', bbox_inches='tight')


# In[21]:


gene_names = high_mutation[:20]
ncomparisons = 20

# create heatmap for male early onset
fig, ax = plt.subplots()
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        corr, pval = stats.spearmanr(me_rna_shared[g1], me_prot_shared[g2], nan_policy = 'omit')
        corr_df.loc[g1, g2] = corr
plot = sns.heatmap(
    corr_df,
    cmap='mako',
    ax = ax
)
fig.suptitle('Male Early Onset')
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig('male_early.png', bbox_inches='tight')


# In[22]:


gene_names = high_mutation[:20]
ncomparisons = 20

# create heatmap for male late onset
fig, ax = plt.subplots()
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        corr, pval = stats.spearmanr(ml_rna_shared[g1], ml_prot_shared[g2], nan_policy = 'omit')
        corr_df.loc[g1, g2] = corr
plot = sns.heatmap(
    corr_df,
    cmap='mako',
    ax = ax
)
fig.suptitle('Male Late Onset')
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)
plt.savefig('male_late.png', bbox_inches='tight')


# In[ ]:


#--------------------------------------------------------------------------
#machine learning
mg = mygene.MyGeneInfo()

# load in data
clinical = pd.read_csv('tcga_coad_clinical_data.csv', index_col=0)
rna = pd.read_csv('tcga_rna_count_data.csv', index_col=0).T
rna.columns = rna.columns.str.split('.').str[0]
clinical = clinical.set_index('Tumor_Sample_Barcode')
ids = []
labels = []
for i in gene_names:
    if i['query'] in rna.columns:
        try:
            labels.append(i['symbol'])
            ids.append(i['query'])
        except:
            pass
rna.rename(columns=dict(zip(ids, labels)), inplace=True)
# select genes to include in ML model features
genes_of_interest = ['SCARNA5', 'SNORA73B','SNORA23','RNU4-2','SNORA74A','RNU4-1','PIWIL1','RNVU1-7','MIR3609','SLC14A1','PLEKHG4B','GTSF1','THPO'] # change this to whichever genes you are using

rna = rna.loc[clinical.index, :]
clinical.loc[:, 'onset'] = np.where(clinical.loc[:, 'age_at_initial_pathologic_diagnosis'] >= 50, 'late', 'early')

rna = rna[~rna.index.duplicated(keep='first')]
ML_data = pd.DataFrame(index = clinical.index)
target = clinical.loc[:, 'onset']

for gene in genes_of_interest:
    try: ML_data.loc[:, f'{gene}_rna'] = rna.loc[:, gene]
    except: print('rna data not found for ' + gene)
        
ML_data = ML_data.dropna(axis=1) # remove columns with NA values, check which columns 
# are present and make sure you have enough data for the tests

# Step 0. Scaling
scaler = StandardScaler()

ML_data = scaler.fit_transform(ML_data)
# Step 1. Dimensionality Reduction
reducer = umap.UMAP()
embedding = reducer.fit_transform(ML_data)


# Plot DR
fig, ax = plt.subplots(constrained_layout=True)
scatters = {}

stages = target.unique()
stages.sort()

stage = 'early'
stage_index = target == 'early'
scatters[stage] = ax.scatter(embedding[stage_index, 0], embedding[stage_index, 1], s=22, alpha = 1)

stage = 'late'
stage_index = target == 'late'
scatters[stage] = ax.scatter(embedding[stage_index, 0], embedding[stage_index, 1], s=12, alpha = 0.4)

fig.legend(scatters, title='Onset', fancybox=True, bbox_to_anchor=[1.15, 0.75]) # May need to play
# around with the bbox_to_anchor args to get legend in right place
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')

# Step 2. Classification
classifier = MLPClassifier()
encoder = LabelEncoder()

encoder.fit(np.sort(target.unique()))
encoded_target = encoder.transform(target)

n_tests = 20
accuracies = []

for n in range(n_tests):
    X_train, X_test, y_train, y_test = train_test_split(ML_data, encoded_target, train_size=0.75)

    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    
    acc = sum(y_pred == y_test)/len(y_pred)
    accuracies.append(acc)

print('mean accuracy:', np.mean(accuracies))
