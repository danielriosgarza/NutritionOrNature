#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:32:34 2019

@author: danielg
"""
#import dependencies
import numpy as np
import os
import umap
import cobra
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



class Env_ball:
    '''
    Environmental ball. Manipulate a predifined uniform random environment.
    
    '''
    def __init__(self, number_of_realizations):
        
        
        
        self.rn = int(number_of_realizations)
        
        
        path=os.getcwd().replace('Scripts', 'files')
        
        self.path=os.path.join(path,'env_ball_' + str(self.rn) +'.tsv') 
        
        with open(self.path) as f:
            self.transporters=f.readline().strip().split('\t')
            
            self.metabolites_d={i.split(':')[1]:i.split(':')[0] for i in f.readline().strip().split('\t')}
            self.dict, self.matrix=self.__get_ds(f)
    
    def __get_ds(self, file_object):
        '''
        build env ball matrix and dict referencing env number and env vector.
        Order is the same as the self.transporters object
        
        '''
        
        counter=0
        d={}
        m=np.zeros((self.rn, len(self.transporters)))
        for line in file_object:
            v=np.array([i for i in map(float, line.strip().split('\t'))])
            m[counter]=v
            d[counter]=v
            counter+=1
        return d, m
    
    def get_main_components(self,n=3):
        '''
        use uMAP to generate a representation of the space in 3 dimensions.
        
        '''
        
        water = self.transporters.index('EX_cpd00001_e')
        oxygen = self.transporters.index('EX_cpd00007_e')
        i=np.arange(len(self.transporters))
        mat = self.matrix.T[(i!=water) & (i!=oxygen)].T
        trans = umap.UMAP(n_neighbors=30, random_state=666, n_components=n, min_dist=.8).fit(mat)
        return trans
    
    def plot(self, color='green', edge_c = 'k', size= 20, line_w = .25, cmap=plt.cm.Greens):
        '''
        plot the three uMAP components in a scatter plot.
        
        '''
        
        trans = self.get_main_components()
        
        fig = plt.figure()
        ax = fig.gca()
        plt.style.use('ggplot')
        ax.scatter(trans.embedding_.T[0], trans.embedding_.T[1], s=size, alpha=.5, c=color, edgecolor=edge_c, lw=line_w, cmap=cmap)
        #fig.colorbar(f)
        plt.setp(ax.get_xticklabels(), fontsize=5)
        plt.setp(ax.get_yticklabels(), fontsize=5)
        
        ax.set_aspect('equal')
        ax.set_xlabel('UMAP_1', fontsize=8)
        ax.set_ylabel('UMAP_2', fontsize=8)
               
               
        
        #ax.set_xticklabels("")
        #ax.set_yticklabels("")
        
    
    def apply_environment(self, mdl, env_vec,transporters):
        for i in range(len(transporters)):
            try:
                mdl.reactions.get_by_id(transporters[i]).lower_bound=-env_vec[i]
                mdl.reactions.get_by_id(transporters[i]).upper_bound=1000.
            except KeyError:
                pass
        sol_m=mdl.optimize()
        return sol_m.objective_value
    
    def get_env_growth_profile(self, mdl_path):
        gv=np.zeros(self.rn)
        mod=cobra.io.read_sbml_model(mdl_path)
        for i in range(self.rn):
            print (i)
            gv[i] = self.apply_environment(mod, self.matrix[i], self.transporters)
        return gv
        
        
        
    
        
    
