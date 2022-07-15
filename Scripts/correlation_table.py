# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:23:34 2020

@author: danie
"""

from pylab import *
from pathlib import Path
import scipy.stats as sts
import scipy.spatial as sps
import pickle
import numpy as np
import os
import umap
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-bright')
from statsmodels.stats.multitest import multipletests as mt

def corr_sig(df=None):
    p_matrix = np.zeros(shape=(df.shape[1],df.shape[1]))
    for col in df.columns:
        for col2 in df.drop(col,axis=1).columns:
            _ , p = sts.pearsonr(df[col],df[col2])
            p_matrix[df.columns.to_list().index(col),df.columns.to_list().index(col2)] = p
    return p_matrix


def plot_cor_matrix(corr, mask=None):
    f, ax = plt.subplots(figsize=(14, 11))
    sns.heatmap(corr, ax=ax,
                mask=mask,
                # cosmetics
                annot=True, vmin=-1, vmax=1, center=0,
                cmap='coolwarm_r', linewidths=.5, square =1, cbar_kws={"shrink": .8, "spacing": "uniform"},  annot_kws={'fontsize':10, 'color':'k'})



data_folder = os.path.join(Path(os.getcwd()).parents[0], 'Files', 'data')

store =  pickle.load(open(data_folder + '/pickles/all_fams.tertiaryDS.pkl', 'rb'))
store.pop('Bacillaceae_plus_Anaplasmataceae', None)


families = np.array(list(store.keys()))


nb={}


with open(os.path.join(data_folder, 'Family_Niche12082020.txt')) as f:
    header = f.readline().strip().split('\t')
    for line in f:
        a=line.strip().split('\t')
        fm = a[0].split('family.')[1]
        if fm in families:
            nb[fm] = float(a[9])
            print(a)




for i in families:
    if i not in nb:
        store.pop(i, None)

families = np.array(list(store.keys()))

#hamming d acros
hamming_d_across={}

for i in families:
    hamming_d_across[i] =  np.std(1-store[i]['hamming_d_within'])/np.std(1-store[i]['hamming_d_across'])

#fir diversity
fir_diversity={}

for i in families:
    print(i)
    mat = np.array(list(store[i]['freq_reactions'].values()))
    
    mat=mat.T
    fir_diversity[i] = np.mean(sps.distance.pdist(mat, metric='sqeuclidean'))    

#fir pan size
fir_pan_size={i: np.mean(store[i]['size_of_pan_firs']) for i in families}

#fir average size
fir_av_size={i: np.mean(store[i]['size_of_firs']['av']) for i in families}

#fir core size
fir_core_size={i: np.mean(store[i]['size_of_core_firs']) for i in families}

#fir shell size
fir_shell_size={i: np.mean(store[i]['size_of_shell_firs']) for i in families}

# fir cloud size
fir_cloud_size={i: np.mean(store[i]['size_of_cloud_firs']) for i in families}

# pan reactome size
reactome_pan_size={i: len(store[i]['freq_reactions']) for i in families}

# Average reactome size
reactome_av_size={i: np.mean(store[i]['size_of_models']) for i in families}
    
# Core reactome size
reactome_core_size={i: store[i]['size_of_core_model_reactome'] for i in families}

# Shell reactome size
reactome_shell_size={i: store[i]['size_of_shell_model_reactome'] for i in families}

#Cloud reactome size
reactome_cloud_size={i: store[i]['size_of_cloud_model_reactome'] for i in families}


#fir fluidity
fir_fluidity={i: np.mean(store[i]['fluidity_firs_across']) for i in families}

# Reactome fluidity
reactome_fluidity={i: np.mean(store[i]['fluidity_models']) for i in families}

# Niche size
niche_size={i: np.mean(store[i]['size_of_niches']['av']) for i in families}


# Number of niche driven reactions
niche_d_r={i: len(store[i]['x_reactome']) for i in families}


# Number of niche driven metabolites
niche_d_m={i: len(store[i]['y_metabolome']) for i in families}

# Number of reactions statistically assciated to metabolites
r_stat_ass={i: len(store[i]['reac_met_freq_association']) for i in families}

#Number of metabolites statistically associated to reactions
met_stat_ass = {}

for i in families:
    met=[]
    for r in store[i]['reac_met_freq_association']:
        for m in store[i]['reac_met_freq_association'][r]['pos']:
            if m not in met:
                met.append(m)
        for m in store[i]['reac_met_freq_association'][r]['neg']:
            if m not in met:
                met.append(m)
    met_stat_ass[i] = len(met)
        
#diference in niches
niche_size={}

for i in families:
    t=np.array(list(store[i]['freq_metabolite_use'].values())).T
    niche_size[i] =  mean(sps.distance.pdist(t, metric='hamming'))
    

# 0 niche breadth (nb)
# 1 fir diversity (fir_diversity)
# 2 fir pan size (fir_pan_size)
# 3 fir average size (fir_av_size)
# 4 fir core size (fir_core_size)
# 5 fir shell size (fir_shell_size)
# 6 fir cloud size (fir_cloud_size)
# 7 pan reactome size (reactome_pan_size)
# 8 Average reactome size (reactome_av_size)
# 9 Core reactome size (reactome_core_size)
# 10 Shell reactome size (reactome_shell_size)
# 11 Cloud reactome size (reactome_cloud_size)
# 12 fir fluidity (fir_fluidity)
# 13 Reactome fluidity (reactome_fluidity)
# 14 Niche size (niche_size)
# 15 Number of niche driven reactions (niche_d_r)
# 16 Number of niche driven metabolites (niche_d_m)

variables = np.zeros((len(families), 17))

for i, f in enumerate(families):
    variables[i][0] = nb[f]
    variables[i][1] = fir_diversity[f]
    variables[i][2] = fir_fluidity[f]
    variables[i][3] = reactome_fluidity[f]
    
    variables[i][4] = fir_pan_size[f]
    variables[i][5] = reactome_pan_size[f]
    
    variables[i][6] = fir_av_size[f]
    variables[i][7] = reactome_av_size[f]
    
    variables[i][8] = fir_core_size[f]
    variables[i][9] = reactome_core_size[f]
    
    variables[i][10] = fir_shell_size[f]
    variables[i][11] = reactome_shell_size[f]
    
    variables[i][12] = fir_cloud_size[f]
    variables[i][13] = reactome_cloud_size[f]
    
    
    
   
    
    variables[i][14] = niche_size[f]
    variables[i][15] = niche_d_r[f]
    variables[i][16] = niche_d_m[f]



names=np.array(['NicheBreadth', 'diversity(pEFMs)','fluidity(pEFMs)','fluidity(Reactomes)', 'pan(pEFMs)', 'pan(Reactomes)','size(pEFMs)', 'size(Reactomes)', 'core(pEFMs)', 'core(Reactomes)', 'shell(pEFMs)', 'shell(Reactomes)', 'cloud(pEFMs)','cloud(Reactomes)','diversity(Metabs)', 'EnvDReacs','EnvDMetabs'])

correl = np.zeros((17,17))
pv = np.zeros((17,17))


for idx, v in enumerate(variables.T):
    for idx2, v2 in enumerate(variables.T):
        pc = sts.pearsonr(v,v2)
        correl[idx][idx2] = np.round(pc[0],3)
        pv[idx][idx2] = pc[1]

#sorter=np.argsort(np.mean(correl, axis=0))[::-1]
#variables = variables.T[sorter].T

#names = names[sorter]

correl = np.zeros((17,17))
pv = np.zeros((17,17))


for idx, v in enumerate(variables.T):
    for idx2, v2 in enumerate(variables.T):
        pc = sts.pearsonr(v,v2)
        correl[idx][idx2] = np.round(pc[0],3)
        pv[idx][idx2] = pc[1]


pv1 = sps.distance.squareform(pv)

adjpv = mt(pv1, method='fdr_bh')[1]
pv = sps.distance.squareform(adjpv)


# with open('C:/Users/danie/Documents/random_work_stuff/home/home/Tables/production_tables/pv.tsv', 'w') as f:
#     for i,v in enumerate(correl):
#         f.write(names[i] + '\t')
#         for z,v2 in enumerate(v):
#             f.write('%.2f (' %v2)
#             f.write("{:.1e}".format(pv[i][z]) +')\t')
            
#         f.write('\n')






vd = {}
for i,v in enumerate(names):
    vd[v] = {}
    for i2,v2 in enumerate(families):
        vd[v][v2] = variables[i2][i]

df = pd.DataFrame.from_dict(vd)


corr = df.corr()                          
p_values = corr_sig(df) 
ti = np.triu_indices_from(p_values)
pvs = p_values[ti].flatten()
corrected = mt(pvs, alpha=0.05, method='fdr_bh')
p_values[ti] = corrected[1]
np.fill_diagonal(p_values, 1)
mask = np.invert(np.tril(p_values<0.05))    
plot_cor_matrix(corr,mask)

savefig(os.path.join(Path(os.getcwd()).parents[0], 'Files', 'figures', 'Figure4.pdf'), dpi=600, tight_layot=1)
