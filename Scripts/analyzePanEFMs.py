#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 06:08:37 2020

@author: daniel
"""

from pathlib import Path
import os
import pickle
import gc
import cobra

import scipy.stats as sts
import scipy.spatial as sps
import numpy as np


from parse_MFS_class import MFS_family
from parse_MFS_class import apply_environment
from Env_ball_class import Env_ball

def cosine(a,b):
    return np.inner(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))



def get_fluidity_index(binary_mat, niter=1000):
    b=np.arange(len(binary_mat))
    sel_ind_1 = np.random.choice(b, replace=True, size=niter)
    sel_ind_2 = np.random.choice(b, replace=True, size=niter)
    
    s1 = np.sum(np.clip(binary_mat[sel_ind_1]-binary_mat[sel_ind_2], 0,1), axis=1)
    s2 = np.sum(np.clip(binary_mat[sel_ind_2]-binary_mat[sel_ind_1], 0,1), axis=1)
    
    fluidity = (s1+s2)/np.sum((binary_mat[sel_ind_1] + binary_mat[sel_ind_2])>0,axis=1)
    return np.mean(fluidity)

def flip_p(m1,m2, niter):
    m = np.arange(m1.shape[1])
    choice_idx = np.random.choice(m, replace=True, size=niter)
    
    k =(m1.T[choice_idx].T - m2.T[choice_idx].T).flatten()
    
    return sum(k!=0)/len(k)

def one_to_all_distance(inde, matrix, metric='hamming'):
    '''
    Returns
    -------
    distance of one vector in a matrix to all other vectores in the matrix.
    '''
    return sps.distance.cdist(matrix[inde].reshape(1,-1),matrix[np.arange(len(matrix))!=inde], metric=metric)

def distance_to_neighbour(matrix, metric='hamming'):
    '''
    Returns
    -------
    distance to closest neigbour
    '''
    v = np.zeros(len(matrix))
    for i in range(len(matrix)):
        v[i] = min(one_to_all_distance(i, matrix, metric)[0])
    return v

def get_even_distances(matrix, metric='hamming'):
    '''
    Buld the even_distance matrix (S)
    1) pick a random vector and add to S
    2) define the walking distance
    3) iterate:
      a) select a random vector that has distance greater then the walking distance to all vectors in S
      b) add selected vector to set
      c) repeat until no vectors are left

    Returns
    -------
    indices of vectors in the matrix that consist of a random sample that 
    don't violate the min walking distance (redundancy)

    '''
    dim = len(matrix)
    selected = [np.random.choice(np.arange(dim))]
    #redundancy threshold
    threshold = max(distance_to_neighbour(matrix))/2.
    a =1
    while a==1:
        k = np.array([i for i in np.arange(dim) if i not in selected])
        if len(k)==dim:
            a=0
        else:
            distance_vs = sps.distance.cdist(matrix[selected], matrix[k], metric=metric)
            m = np.ones(len(k))
            for i in distance_vs:
                m*=i>threshold
            if sum(m)==0:
                a=0
            else:
                selected.append(np.random.choice(k[m.astype(np.bool)]))
    print('done')
    return selected


def main(family):
    #reference from the git script
    data_folder = os.path.join(Path(os.getcwd()).parents[1], 'data')

    #load pickles
    fam_mfs=pickle.load(open(data_folder + '/pickles/' + family + '.mfs.pkl', 'rb'))
    
    fam_associate =pickle.load(open(data_folder + '/pickles/' + family + '.associate.pkl', 'rb'))

    #obtain all the primary data structures (see: Measures & data structures)
    
    #Reactome:[r]
    
    reactome = fam_associate['reactome'].copy()
    
    # Metabolome: [M]
    
    metabolome = fam_associate['transporter'].copy()
    
    # Environment Ball: [N1][M]
    
    environment_ball_class =Env_ball(1000) 
    
    #exclude oxygen and water
    met_idx = np.array([environment_ball_class.transporters.index(i) for i in metabolome])
    environment_ball = environment_ball_class.matrix.T[met_idx].T
    
    #FIRS: [E][N1][r]
           
    #get only the included reactions
    
    firs = {i:fam_mfs.mfs[i][fam_mfs.include_reactome].T for i in fam_mfs.mfs}
    
    # Niche: [E][N1][M]
    
    niche = fam_associate['used_env'].copy()
    
    
    
    #shuffle sample orders
    for i in range(1000):
        r=np.arange(1000)
        np.random.shuffle(r)
        firs[str(i)] = firs[str(i)][r]
        niche[i] = niche[i][r]
    
    # Niche binary: [E][N1][M]
    niche_binary ={}
    
    for i in niche:
        niche_binary[i] = niche[i].copy()
        niche_binary[i][np.round(niche_binary[i], 10) !=0] =1.0
    
    #Models: [s]
    
    models = np.array([i.replace('.sbml','') for i in fam_mfs.model_files])
    
    # Model reactomes: [s][r]
    model_reactomes = fam_mfs.model_reactomes.copy()
    model_reactomes = model_reactomes.T[fam_mfs.include_reactome].T
    
    #Model sample: [d][r]
    
    model_sample_idx = np.array(get_even_distances(model_reactomes, metric='hamming'))
    model_sample = model_reactomes[model_sample_idx]
    
    # FIRS growth rate: [E][N1]
    
    firs_growth_rate = np.zeros((1000, 1000))
    
    for i in range(1000):
        firs_growth_rate[i] = np.sum(niche[i], axis=1)
    
    
    
    #remove CO2 and H+
        
    met_idx =(metabolome!='EX_cpd00011_e') & (metabolome!='EX_cpd00067_e')
    metabolome = metabolome[met_idx]
    environment_ball = environment_ball.T[met_idx].T
    
    for i in niche:
        niche[i] = niche[i].T[met_idx].T
    
    for i in niche_binary:
        niche_binary[i] = niche_binary[i].T[met_idx].T
    
    ######Secondary Data Structures###
    
    
    
    
    #Size of FIRS: [E][N1]
    size_of_firs = np.zeros((1000, 1000))
    for i in range(1000):
        size_of_firs[i] = np.sum(firs[str(i)], axis=1)
    
    #Size of Niche: [E][N1]
    size_of_niches = np.zeros((1000,1000))
    for i in range(1000):
        size_of_niches[i] = np.sum(niche_binary[i], axis=1)
    
    #Size of models: [s]
    size_of_models = np.zeros(len(models))
    for i,v in enumerate(model_reactomes):
        size_of_models[i] = sum(v>0)
    
    #Fluidity of FIRS within environments: [E]
    
    fluidity_firs_within = np.zeros(1000)
    for i in range(1000):
        fluidity_firs_within[i]=get_fluidity_index(firs[str(i)], 1000)
    
    #Fluidity of FIRS across environments: [N2]
    
    fluidity_firs_across = np.zeros(10000)
    
    for i in range(10000):
        
        rintA = np.random.randint(0,1000, size=2)
        rintB = np.random.randint(0,1000, size=2)
        s1 = sum(np.clip(firs[str(rintA[0])][rintB[0]]-firs[str(rintA[1])][rintB[1]], 0,1))
        s2 = sum(np.clip(firs[str(rintA[1])][rintB[1]]-firs[str(rintA[0])][rintB[0]], 0,1))
        fluidity_firs_across[i] = (s1+s2)/sum((firs[str(rintA[0])][rintB[0]]+firs[str(rintA[1])][rintB[1]])>0)
    
    #Fluidity of niches: [E]
    fluidity_niche_within = np.zeros(1000)
    for i in range(1000):
        
        fluidity_niche_within[i]=get_fluidity_index(niche_binary[i], 1000)
    
    #Fluidity across niches: [N2]
    fluidity_niche_across = np.zeros(10000)
    
    for i in range(10000):
        
        rintA = np.random.randint(0,1000, size=2)
        rintB = np.random.randint(0,1000, size=2)
        s1 = sum(np.clip(niche_binary[rintA[0]][rintB[0]]-niche_binary[rintA[1]][rintB[1]], 0,1))
        s2 = sum(np.clip(niche_binary[rintA[1]][rintB[1]]-niche_binary[rintA[0]][rintB[0]], 0,1))
        fluidity_niche_across[i] = (s1+s2)/sum((niche_binary[rintA[0]][rintB[0]]+niche_binary[rintA[1]][rintB[1]])>0)
    
    #Fluidity of models: [N2]
    
    fluidity_models = np.zeros(10000)
    for i in range(10000):
        print(i)
        fluidity_models[i]=get_fluidity_index(model_reactomes, 2)
        
        
    #Fluidity of model samples: [N2]
    fluidity_model_samples = np.zeros(10000)
    for i in range(10000):
        print(i)
        fluidity_model_samples[i]=get_fluidity_index(model_sample, 2)
    
    #Frequency of reactions: [E][r]
    
    freq_reactions = np.zeros((1000,len(reactome)))
    for i in range(1000):
        freq_reactions[i] = np.sum(firs[str(i)], axis=0)/1000
    
    #Residual reactions frequency: [E][r]
    freq_reactions_m = np.mean(freq_reactions, axis=0)    
    
    residual_reaction_freq = freq_reactions - freq_reactions_m
    
    #niche driven score for reactions: [r]
    
    niche_score_reactions = np.round(np.std(residual_reaction_freq, axis=0),5)
    
    #Reaction frequency in models: [r]
    freq_mod_reactions = np.sum(model_reactomes, axis=0)/len(models)
    
    #Reaction frequency in model sample[r]
    freq_mod_samp_reactions = np.sum(model_sample, axis=0)/ len(model_sample)
    
    #Metabolite usage frequency: [E][M]
    
    freq_metabolite_use = np.zeros((1000, len(metabolome)))
    for i in range(1000):
        freq_metabolite_use[i] = np.sum(niche_binary[i], axis=0)/1000
    freq_metabolite_use_m = np.mean(freq_metabolite_use, axis=0)
    
    residual_metabolite_freq = freq_metabolite_use - freq_metabolite_use_m
    
    #metabolite usage flux [E][M]
    
    metabolite_usage_flux = np.zeros((1000, len(metabolome)))
    for i in range(1000):
        metabolite_usage_flux[i] = np.sum(niche[i], axis=0)/1000
    
    #residual metabolite usage flux: [E][M]
    
    metabolite_usage_flux_m = np.mean(metabolite_usage_flux, axis=0)
    residual_metabolite_usage_flux = metabolite_usage_flux-metabolite_usage_flux_m
    
    #niche driven score for metabolites: [M]
    niche_score_metabolites = np.round(np.std(residual_metabolite_usage_flux, axis=0),5)
    
    #####x, y: non zero reactions frequencies and metabolites usage flux####
    
    
    x_reactome = reactome[niche_score_reactions!=0]
    x_reac_freq = freq_reactions.T[niche_score_reactions!=0].T
    
    
    y_metabolome = metabolome[niche_score_metabolites!=0]
    y_met_usage_flux = metabolite_usage_flux.T[niche_score_metabolites!=0].T
    
    y_met_freq =freq_metabolite_use.T[niche_score_metabolites!=0].T
    
    #correlation: [r][M]
    
    correlation=np.zeros((len(x_reactome), len(y_metabolome)))
    
    for i, reac in enumerate(x_reac_freq.T):
        correlation[i] = np.array([sts.pearsonr(reac.flatten(), metab.flatten())[0] for metab in y_met_usage_flux.T])
    
    #correlation metabolite frequency
        
    correlation_met_freq=np.zeros((len(x_reactome), len(y_metabolome)))
    
    for i, reac in enumerate(x_reac_freq.T):
        correlation_met_freq[i] = np.array([sts.pearsonr(reac.flatten(), metab.flatten())[0] for metab in y_met_freq.T])
    
    
    #Reaction pairwise distance: [E][E]
    reaction_pairwise_distance = sps.distance.squareform(sps.distance.pdist(freq_reactions))
    
    #FIRS pairwise distance: [E]
    firs_pairwise_distance = np.zeros(1000)
    
    for i in range(1000):
        firs_pairwise_distance[i] = np.mean(sps.distance.pdist(firs[str(i)], metric='hamming'))
    
    #Niche binary pairwise distance: [E]
    niche_binary_pairwise_distance = np.zeros(1000)
    
    for i in range(1000):
        niche_binary_pairwise_distance[i] = np.mean(sps.distance.pdist(niche_binary[i], metric='hamming'))
        
    #Niche pairwise distance: [E]
    niche_pairwise_distance = np.zeros(1000)
    
    for i in range(1000):
        niche_pairwise_distance[i] = np.mean(sps.distance.pdist(niche[i]))
    
    #Niche distance: [E][E]
    niche_distance = sps.distance.squareform(sps.distance.pdist(metabolite_usage_flux))
    
    #DNDS_reaction: [N1]
    dn_reactions = np.zeros(1000)
    ds_reactions = np.zeros(1000)
    rand_idx =np.arange(1000)
    for i in range(1000):
        f1 = np.random.randint(0,1000,size=2)
        np.random.shuffle(rand_idx)
        ds_reactions[i] = flip_p(firs[str(f1[0])], firs[str(f1[0])][rand_idx],1000)
        dn_reactions[i] = flip_p(firs[str(f1[0])], firs[str(f1[1])],1000)
    
    store ={'size_of_firs': size_of_firs, 'size_of_niches': size_of_niches,\
            'size_of_models': size_of_models, 'fluidity_firs_within': fluidity_firs_within,\
            'fluidity_firs_across':fluidity_firs_across, 'fluidity_niche_within':fluidity_niche_within,\
            'fluidity_niche_across': fluidity_niche_across, 'fluidity_models': fluidity_models,\
            'fluidity_model_samples': fluidity_model_samples, 'freq_reactions': freq_reactions, \
            'residual_reaction_freq': residual_reaction_freq, 'niche_score_reactions': niche_score_reactions,\
            'freq_mod_reactions': freq_mod_reactions, 'freq_mod_samp_reactions': freq_mod_samp_reactions,\
            'freq_metabolite_use': freq_metabolite_use, 'metabolite_usage_flux': metabolite_usage_flux,\
            'metabolite_usage_flux_m': metabolite_usage_flux_m, 'niche_score_metabolites':niche_score_metabolites,\
            'x_reactome': x_reactome, 'x_reac_freq':x_reac_freq,\
            'y_metabolome' : y_metabolome, 'y_met_usage_flux':y_met_usage_flux,\
            'y_met_freq': y_met_freq, 'correlation': correlation, \
            'correlation_met_freq':correlation_met_freq, 'reaction_pairwise_distance':reaction_pairwise_distance,\
            'firs_pairwise_distance': firs_pairwise_distance, 'niche_binary_pairwise_distance': niche_binary_pairwise_distance,\
            'niche_pairwise_distance': niche_pairwise_distance, 'niche_distance': niche_distance,\
            'dn_reactions': dn_reactions, 'ds_reactions':ds_reactions}
    pickle.dump(store, open(data_folder + '/pickles/' + family + '.secondaryDS.pkl', 'wb'))
    
fams = os.listdir('/mnt/comics/danielg/metabolic_models/git/data')

fams.remove('pickles')
import gc
for fam in fams:
    main(fam)
    gc.collect()
