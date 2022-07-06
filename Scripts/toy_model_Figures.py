#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 10:53:12 2020

@author: daniel
"""

from pylab import *
import os
from pathlib import Path
import numpy as np
import pickle
import seaborn as sns
import scipy.stats as sts
import matplotlib.pyplot as plt

def get_residual_scores(matrix):
    n=len(matrix)
    av_mat = np.round(np.mean(matrix, axis=0)*n)
    mat = np.round(matrix*n)
    r1  = (mat - av_mat)
    score = np.std(r1,axis=0)
    
    stand_r1 = np.zeros(r1.shape)
    
    return score/max(score)


#get the data
data_folder = os.path.join(Path(os.getcwd()).parents[0], 'Files', 'data')
store = pickle.load(open(data_folder + '/pickles/toy_model.pkl', 'rb'))


####labels###
reactions = np.array(['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15'])
mets = np.array(['M1_e','M2_e','M3_e','M4_e','M5_e','M6_e','M7_e','M8_e','M9_e','M10_e'])

####sorters: reaction frequency and metabolite usage ####

reaction_freq = store['full_freq_m']
s_reaction_freq = np.sum(reaction_freq, axis=1)
s0_reaction_freq = np.sum(reaction_freq, axis=0)
sorter_reac = np.argsort(s_reaction_freq)
sorter_reac = sorter_reac[::-1]


used_mets = store['used_environment']
s_mets = np.sum(used_mets, axis=1)
s0_mets = np.sum(used_mets, axis=0)
sorter_mets = np.argsort(s_mets)
sorter_mets = sorter_mets[::-1]



####Classifying environment_drivers###

residual_met_u = store['diff_used_env']

resid_mets_rank = np.zeros(residual_met_u.shape)

for i,v in enumerate(residual_met_u):
    rk = sts.rankdata(np.abs(np.round(v,2)), method='min')-1
    resid_mets_rank[i] = rk/max(rk)
    
unormed_met_score=np.zeros(len(mets))

for i,v in enumerate(resid_mets_rank.T):
    unormed_met_score[i] = sum(resid_mets_rank.T[i]==1)
    
env_driv_met_score = unormed_met_score/len(resid_mets_rank)
env_driv_met_score = get_residual_scores(used_mets)
sorter_mets0 = np.argsort(env_driv_met_score)

#####classifying driven reactions ######

residual_reac_f = store['diff_freq_m']

resid_reac_rank = np.zeros(residual_reac_f.shape)

for i,v in enumerate(residual_reac_f):
    rk = sts.rankdata(np.abs(np.round(v,2)), method='min')-1
    resid_reac_rank[i] = rk/max(rk)
    
unormed_reac_score=np.zeros(len(reactions))

for i,v in enumerate(resid_reac_rank.T):
    unormed_reac_score[i] = sum(resid_reac_rank.T[i]==1)
    

env_driv_reac_score = get_residual_scores(reaction_freq)
sorter_reac0 = np.argsort(env_driv_reac_score)

######residual reaction frequencies###


sns.heatmap(residual_reac_f[sorter_mets].T[sorter_reac0].T, cmap=cm.coolwarm); 


xticks(np.arange(0.5, len(reactions)+0.5), reactions[sorter_reac0], rotation=85, fontsize=12)
yticks([0, 100,208],['208', '100','0'], fontsize=12)

for i in range(15):
    vlines(i,0,208,color='w', lw =.5)
#savefig(pathToSave, dpi=600)
show()




####Association Matrix #####

ass_d = store['cosine_dict']

u_reacts = reactions[store['reac_filter1']]
u_reacts = u_reacts[store['reac_filter2']]



u_mets = mets[store['met_filter1']]
u_mets = u_mets[store['met_filter2']]

sorter_u_mets=np.zeros(len(u_mets), np.int)

c=-1
for i,v in enumerate(mets[sorter_mets0]):
    if v in u_mets:
        c+=1
        sorter_u_mets[c] = np.arange(len(u_mets))[u_mets==v]

sorter_u_reac=np.zeros(len(u_reacts), np.int)

c=-1
for i,v in enumerate(reactions[sorter_reac0]):
    if v in u_reacts:
        c+=1
        sorter_u_reac[c] = np.arange(len(u_reacts))[u_reacts==v]



ass_m = np.array([ass_d[m] for m in store['envd_reactions']])


sns.heatmap(ass_m[sorter_u_reac].T[sorter_u_mets].T, cmap=cm.coolwarm, linewidths=0.1, linecolor='w')

xticks(np.arange(0.5, len(u_mets)+0.5), u_mets[sorter_u_mets], rotation=85, fontsize=12)
yticks(np.arange(0.5, len(u_reacts)+0.5), u_reacts[sorter_u_reac], rotation=0, fontsize=12)


#savefig(pathToSave, dpi=600)
show()



# #####evolved metabolite usage #########
met_usage = store['true_used_env'].copy()


# #####predicted metabolite usage #########
pred_met_usage = store['predicted_environment'].copy()


####predict evolved####

sns.regplot(met_usage.flatten(), pred_met_usage.flatten(), scatter_kws={'s':5}, line_kws={'color':'g', 'lw':.5})
#savefig(pathToSave, dpi=600)
show()
