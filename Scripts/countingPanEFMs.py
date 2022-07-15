# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 21:27:11 2022

@author: u0139894
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

core_panEFMs = []
for i in families:
    core_panEFMs.append(np.mean(store[i]['size_of_core_firs']/len(store[i]['x_reactome'])))


shell_panEFMs = []
for i in families:
    shell_panEFMs.append(np.mean(store[i]['size_of_shell_firs']/len(store[i]['x_reactome'])))


cloud_panEFMs = []
for i in families:
    cloud_panEFMs.append(np.mean(store[i]['size_of_cloud_firs']/len(store[i]['x_reactome'])))



