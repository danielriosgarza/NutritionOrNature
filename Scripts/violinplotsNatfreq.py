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
    panr = d['freq_mod_reactions'].copy()
    
    x_efm = []
    y_mod = []
    
    for i in reacs:
        mod_i = [panr[i]]*1000
        y_mod += mod_i
        x_efm += list(efms[i])
        
    sample = np.random.choice(np.arange(len(x_efm)), replace=0, size = 100000)
    
    x = np.array(x_efm)[sample]
    y = np.array(y_mod)[sample]
    corr=sts.pearsonr(x,y)
    
    tx = ' (Pearsonr = ' + '{:#.2g}'.format(corr[0]) + '; adjP =' + '{:.1e}'.format(corr[1]*46) + ')'
    print(tx)

    
    return x,y,tx

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
pdf = PdfPages(os.path.join(Path(os.getcwd()).parents[0], 'Files', 'figures', 'outNatFreq.pdf'))
 
# Generate the pages
nb_plots = len(data)
nb_plots_per_page = 8
nb_pages = int(np.ceil(nb_plots / float(nb_plots_per_page)))
grid_size = (nb_plots_per_page, 2)
 
for i, samples in enumerate(data):
    x,y, tx = reaction_violin_plot(store, samples)
    # Create a figure instance (ie. a new page) if needed
    if i % nb_plots_per_page == 0:
        fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
    
    # Plot stuffs !
    if i%2==0:
        plt.subplot2grid(grid_size, (i % nb_plots_per_page, 0))
    else:
        plt.subplot2grid(grid_size, (i % nb_plots_per_page-1, 1))
    
    
    sns.kdeplot(x, y, color='r', shade=True, cmap="Reds", shade_lowest=False)
    plt.title(str(i+1)+ '. ' + samples + '\n' + tx, fontsize=10)
    plt.xlabel('Frequency in panEFMs', fontsize=8)
    plt.ylabel('Frequency in pan-reactomes', fontsize=8)
    
    # Close the page if needed
    if (i + 1) % nb_plots_per_page == 0 or (i + 1) == nb_plots:
        
        pdf.savefig(fig)
 
# Write the PDF document to the disk
pdf.close()



    