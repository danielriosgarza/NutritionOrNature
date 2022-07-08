#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:37:21 2020

@author: daniel
"""

from pathlib import Path
import os
import pickle


from pylab import *
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-darkgrid')
import scipy.stats as sts
import scipy.spatial as sps
import seaborn as sns
import pandas as pd



def betafit(data):
    mean=np.mean(data)
    var=np.var(data,ddof=1)
    alpha1=mean**2*(1-mean)/var-mean
    beta1=alpha1*(1-mean)/mean
    return alpha1, beta1

#import the basic ds
data_folder = os.path.join(Path(os.getcwd()).parents[0], 'Files', 'data')
fam_TDS=pickle.load(open(data_folder + '/pickles/all_fams.tertiaryDS.pkl', 'rb'))
families = np.array(list(fam_TDS.keys()))

##Characteristic FIR size####

fir_size= np.array([fam_TDS[i]['size_of_firs']['av'] for i in families])

fir_size_m = np.mean(fir_size, axis=1)

sorter = np.argsort(fir_size_m)

boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)
boxplot(fir_size[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)

show()


###average_model_size#####

mod_sizes=[]
av_mod_size=np.zeros(len(families))

for i in range(len(families)):
    av_mod_size[i] = median(fam_TDS[families[i]]['size_of_models'])
    mod_sizes.append(fam_TDS[families[i]]['size_of_models'])
mod_sizes = np.array(mod_sizes)   

boxplot(mod_sizes[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops, showmeans=True)
xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)

show()


#####pan reactome size #############
pan_reactome_size=np.zeros(len(families))

for i,v in enumerate(families):
    pan_reactome_size[i] = len(fam_TDS[v]['freq_reactions'])


###environment associated reaction size ######
env_associate_reac_size = np.zeros(len(families))
for i,v in enumerate(families):
    l = np.round(np.array([fam_TDS[v]['niche_score_reactions'][i] for i in fam_TDS[families[i]]['niche_score_reactions']]),2)
    env_associate_reac_size[i] = sum(l>0)


##Characteristic fluidity within environments####

fluidity_within= np.array([fam_TDS[i]['fluidity_firs_within'] for i in families])

fluidty_within_m = np.median(fluidity_within, axis=1)


boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)
boxplot(fluidity_within[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)
show()

    

##Characteristic Niche size####

niche_size= np.array([fam_TDS[i]['size_of_niches']['av'] for i in families])

niche_size_m = np.median(niche_size, axis=1)

sorter=np.argsort(niche_size_m)
boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)

boxplot(niche_size[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)

show()

##Characteristic fluidity across environments####

fluidity_across= np.array([fam_TDS[i]['fluidity_firs_across'] for i in families])

fluidty_across_m = np.median(fluidity_across, axis=1)


boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)
boxplot(fluidity_across[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)

show()


##Characteristic dn across environments####

dn= np.array([fam_TDS[i]['dn_reactions'] for i in families])

dn_m = np.median(dn, axis=1)


boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)
boxplot(dn[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)
show()


##Characteristic ds across environments####

ds= np.array([fam_TDS[i]['ds_reactions'] for i in families])

ds_m = np.median(ds, axis=1)


boxprops = dict(linestyle='-', linewidth=.5, color='k')
medianprops = dict(linestyle='-', linewidth=1.0, color='green')
flierprops = dict(marker='o', markerfacecolor='green', markersize=8,
                  markeredgecolor='w', alpha=.9)
boxplot(ds[sorter].T, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops)

xticks(ticks=np.arange(1,len(families)+1), labels=families[sorter], rotation=90)
show()



niche_s_r={}


for i in families:
    for z in fam_TDS[i]['niche_score_reactions']:
        if z not in niche_s_r:
            niche_s_r[z]=[]
        niche_s_r[z].append(fam_TDS[i]['niche_score_reactions'][z])




r_n_ass={}

for i in families:
    for z in fam_TDS[i]['reac_met_flux_association']:
        if z not in r_n_ass:
            r_n_ass[z]={'n':0}
        r_n_ass[z]['n'] +=1
        for met in fam_TDS[i]['reac_met_flux_association'][z]:
            if met not in r_n_ass[z]:
                r_n_ass[z][met]=0
        for met in fam_TDS[i]['reac_met_flux_association'][z]:
            r_n_ass[z][met]+=1
            
        
            

