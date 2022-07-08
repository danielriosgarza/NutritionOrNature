from pylab import *
import scipy.stats as sts
import scipy.spatial as sps
from pathlib import Path
import pickle
import numpy as np
import os
import umap
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')


data_folder = os.path.join(Path(os.getcwd()).parents[0], 'Files', 'data')


phylum={}
order={}
with open(os.path.join(data_folder, 'strainData.txt')) as f:
    h=f.readline().strip().split('\t')
    for line in f:
        a=line.strip().split('\t')
        phylum[a[9]]=a[6]
        order[a[9]] = a[8]


def make_umap(data_m, decomposition=None, n_neighbors=20, min_dist=.1, metric='cosine'):
    if decomposition is None:
        pc = data_m.copy()
    else:
        from sklearn.decomposition import PCA
        pca=PCA(n_components=decomposition)
        pc=pca.fit_transform(data_m)
    
    import umap
    trans = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric).fit(pc)
    
    return trans


store = pickle.load(open(data_folder + '/pickles/all_fams.tertiaryDS.pkl', 'rb'))
store.pop('Bacillaceae_plus_Anaplasmataceae', None)




families = np.array(list(store.keys()))



nb={}


with open(os.path.join(data_folder, 'family.niche_breadth_score.tsv')) as f:
    header = f.readline().strip().split('\t')
    for line in f:
        a=line.strip().split('\t')
        fm = a[0].split('family.')[1]
        if fm in families:
            nb[fm] = float(a[-1])




for i in families:
    if i not in nb:
        store.pop(i, None)

families = np.array(list(store.keys()))


reac_dim = (1000*len(families))+len(families)
all_reactions={}

for i in families:
    
    for reaction in store[i]['reac_met_freq_association']:
        if reaction not in all_reactions:
            all_reactions[reaction] = np.zeros(reac_dim)

for i, fam in enumerate(families):
    for reaction in all_reactions:
             
        if reaction in store[fam]['reac_met_freq_association']:
            all_reactions[reaction][i*1000: (i*1000+1000)] = store[fam]['freq_reactions'][reaction]
            if sum(store[fam]['freq_reactions'][reaction])>1:
                all_reactions[reaction][len(families)*1000 + i] = store[fam]['freq_mod_samp_reactions'][reaction]
            else:
                all_reactions[reaction][len(families)*1000 + i] = np.mean(store[fam]['freq_reactions'][reaction])
            
            
reactions = np.array(list(all_reactions.keys()))

data = np.zeros((len(reactions), reac_dim))

for i,v in enumerate(reactions):
    data[i] = all_reactions[v]

data_reactions=data.T

# for i,v in enumerate(data_reactions):
#     c = np.random.uniform(size = sum(v==0))
#     data_reactions[i][data_reactions[i]==0] = c
    

###UMAP Reactions#####


reactions_umap =  make_umap(data_reactions, n_neighbors=200, min_dist=.5, metric='minkowski')




c2=np.arange(len(families))
np.random.shuffle(c2)
c2=c2/(max(c2) +1)
cmap = matplotlib.cm.get_cmap('tab20')

for i,fam in enumerate(families):
    scatter(reactions_umap.embedding_.T[0][i*1000:(i*1000 + 1000)], reactions_umap.embedding_.T[1][i*1000:(i*1000 + 1000)], c=np.array([cmap(c2[i])]), alpha=1, s=20, lw=.2,edgecolors='w')
    scatter(reactions_umap.embedding_.T[0][1000*len(families)+i], reactions_umap.embedding_.T[1][1000*len(families)+i], c=np.array([cmap(c2[i])]), alpha=.7, s=100, zorder=50, lw=1,edgecolors='w',label=fam+'(' + str(i+1) +')')
    text(reactions_umap.embedding_.T[0][1000*len(families)+i], reactions_umap.embedding_.T[1][1000*len(families)+i], str(i+1), zorder=100)

legend(loc='lower right', ncol=4,bbox_to_anchor=(1.1, -1.0))


#savefig(pathToFig, dpi =300,bbox_inches="tight")
show()

############################################

with open(os.path.join(data_folder, 'EN_pred.tsv')) as f:
    fms = f.readline().strip().split('\t')[1::]
    
    en_preds={i:{} for i in fms}
    
    for line in f:
        a=line.strip().split('\t')
        mt = a[0]
        vs=list(map(float, a[1::]))
        for z,v in enumerate(fms):
            en_preds[v][mt] = vs[z]

        




all_metabolites={i: np.zeros(reac_dim) for i in en_preds[families[0]]}


for i, fam in enumerate(families):
    for met in all_metabolites:
        if met in store[fam]['freq_metabolite_use']:
            all_metabolites[met][i*1000: (i*1000+1000)] = store[fam]['freq_metabolite_use'][met]
            all_metabolites[met][len(families)*1000 + i] = en_preds[fam][met]

metabolites = np.array(list(all_metabolites.keys()))

data = np.zeros((len(metabolites), reac_dim))

for i,v in enumerate(metabolites):
    data[i] = all_metabolites[v]

data_metabolites=data.T

######Metabolite UMAP####################################

metabolites_umap =  make_umap(data_metabolites,n_neighbors=20, min_dist=.5, metric='minkowski')

for i,fam in enumerate(families):
    scatter(metabolites_umap.embedding_.T[0][i*1000:(i*1000 + 1000)], metabolites_umap.embedding_.T[1][i*1000:(i*1000 + 1000)], c=np.array([cmap(c2[i])]), alpha=1, s=20, lw=.2,edgecolors='w' )
    scatter(metabolites_umap.embedding_.T[0][1000*len(families)+i], metabolites_umap.embedding_.T[1][1000*len(families)+i], c=np.array([cmap(c2[i])]), alpha=1, s=100, zorder=50, lw=1,edgecolors='w',label=fam+'(' + str(i+1) +')')
    text(metabolites_umap.embedding_.T[0][1000*len(families)+i], metabolites_umap.embedding_.T[1][1000*len(families)+i], str(i+1), zorder=100)

legend(loc='lower right', ncol=4,bbox_to_anchor=(1.1, -1.0))

#savefig(pathToFig, dpi =300,bbox_inches="tight")

show()

######################################################

# d={i:[] for i in range(25)}

# for i,v in enumerate(families):
#     r = np.array(list(store[v]['freq_reactions'].keys()))
#     freqM = np.array([store[v]['freq_mod_reactions'][z] for z in r])
    
#     freqF = np.array([mean(store[v]['freq_reactions'][z]) for z in r])
    
    
#     count=-1
#     for mr in linspace(0,1,25):
#         count+=1
#         l=list(freqM[(freqF>=mr) & (freqF<(mr+0.15))])
#         d[count]+=l
    
# boxprops = dict(linestyle='-', linewidth=.5, color='k')
# medianprops = dict(linestyle='-', linewidth=2.0, color='green')
# flierprops = dict(marker='o', markerfacecolor='w', markersize=3, markeredgecolor='k',  markeredgewidth=.1)
# meanprops = dict(markerfacecolor='r',markeredgecolor='r' )
    
    
# b=np.array([d[i] for i in range(25)]) 
# boxplot(b, boxprops =boxprops, medianprops=medianprops, flierprops=flierprops,showmeans=True, meanprops=meanprops)

# ts =[]
# for i in np.round(linspace(0,1,25),2):
#     ts.append('%.2f'%i)

# ts2 =[]
# for i in np.round(linspace(0,1,8),2):
#     ts2.append('%.2f'%i)


# xticks(np.arange(1,26), ts, fontsize=10, rotation=90)
# yticks(linspace(0,1,8), ts2, fontsize=10)

# #savefig(pathToFig, dpi =300,bbox_inches="tight")

# show()

