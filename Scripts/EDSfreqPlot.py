# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:08:39 2020

@author: danie
"""


# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 14:36:06 2020

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

def reaction_violin_plot(store_3, family):
    d = store_3[family]
    reacs = np.array(list(d['freq_reactions'].keys()))
    efms = d['freq_reactions'].copy()
    eds = d['niche_score_reactions'].copy()
    
    x_efm = []
    y_eds = []
    x_efm_0 = []
    y_eds_0 = []
    
    for i in reacs:
        #if eds[i]>0.001:
        eds_i = [eds[i]]*1000
        y_eds += eds_i
        x_efm += list(efms[i])
        # else:
        #     eds_i = [eds[i]]*1000
        #     y_eds_0 += eds_i
        #     x_efm_0 += list(efms[i])
        
    sample = np.random.choice(np.arange(len(x_efm)), replace=0, size = 10000)
    #sample0 = np.random.choice(np.arange(len(x_efm_0)), replace=0, size = 10000)

    
    x = np.array(x_efm)[sample]
    y = np.array(y_eds)[sample]
    #x0 = np.array(x_efm_0)[sample0]
    #y0 = np.array(y_eds_0)[sample0]
    
    
    return x,y#,x0, y0

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


    

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Generate the data
data = families.copy()
 
# The PDF document
pdf = PdfPages(os.path.join(Path(os.getcwd()).parents[0], 'Files', 'figures', 'outEDSFreq.pdf'))
 
# Generate the pages
nb_plots = len(data)
nb_plots_per_page = 8
nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
grid_size = (nb_plots_per_page, 2)
 
for i, samples in enumerate(data):
    x,y = reaction_violin_plot(store, samples)
    # Create a figure instance (ie. a new page) if needed
    if i % nb_plots_per_page == 0:
        fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    
    # Plot stuffs !
    if i%2==0:
        plt.subplot2grid(grid_size, (i % nb_plots_per_page, 0))
    else:
        plt.subplot2grid(grid_size, (i % nb_plots_per_page-1, 1))
    
    
    sns.kdeplot(x, y, color='b', shade=True, cmap="Blues", shade_lowest=False)
    #sns.kdeplot(x0, y0, color='r', shade=True, cmap="Reds", shade_lowest=False)
    plt.title(str(i+1)+ '. ' + samples, fontsize=10)
    plt.xlabel('Frequency in panEFMs', fontsize=8)
    plt.ylabel('EDS', fontsize=8)
    
    # Close the page if needed
    if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
        
        pdf.savefig(fig)
 
# Write the PDF document to the disk
pdf.close()



    