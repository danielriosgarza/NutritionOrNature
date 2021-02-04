# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:14:03 2020

@author: danie
"""

from pathlib import Path
import os
import pickle
import gc

import scipy.stats as sts
import scipy.spatial as sps
import numpy as np








data_folder = os.path.join(Path(os.getcwd()))

#file_folder =  os.path.join(Path(os.getcwd()).parents[0], 'files')
fam_panEFM=pickle.load(open(data_folder + '/all_fams.tertiaryDS.pkl', 'rb'))


families = np.array(list(fam_panEFM.keys()))



metab_d={}

for i,v in enumerate(families):
    
    x, y, z = [], [], []
    mets=[]
    
    for r in fam_panEFM[v]['reac_met_freq_association']:
        x.append(fam_panEFM[v]['freq_reactions'][r])
    
    sel =np.random.choice(np.arange(1000), replace=False, size=200)
    x = np.array(x).T
    
    for m in fam_panEFM[v]['y_metabolome']:
        y.append(fam_panEFM[v]['freq_metabolite_use'][m])
        mets.append(m)
        
    
    y = np.array(y).T
    
    for r in fam_panEFM[v]['reac_met_freq_association']:
        z.append(fam_panEFM[v]['freq_mod_reactions'][r])
    z=np.array(z)
    from sklearn.linear_model import MultiTaskElasticNetCV as en
    enet=en(cv=2, n_jobs=5,verbose=1, max_iter=10000)
    EN=enet.fit(x[sel],y[sel])
    p=EN.predict(z.reshape(1,-1)).flatten()
    d={mets[i]:p[i] for i in range(len(mets))}
    
    metab_d[v] = d.copy()
    

en_res={}

for fam in metab_d:
    for m in metab_d[fam]:
        if m not in en_res:
            en_res[m]=np.zeros(len(families))

for i,fam in enumerate(families):
    for m in metab_d[fam]:
        en_res[m][i] = metab_d[fam][m]

with open('C:\\Users\\danie\\Documents\\random_work_stuff\EN_pred.tsv', 'w') as f:
    f.write('metabolite\t')
    f.write('\t'.join(families))
    f.write('\n')
    
    for i in en_res:
        f.write(i+'\t')
        for z in en_res[i]:
            f.write('%.6f\t' %z)
        f.write('\n')
        
    
