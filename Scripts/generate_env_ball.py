#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:11:38 2019

@author: danielg
"""
import cobra
import numpy as np
import os
import os
from pathlib import Path

def main(simul_n, transporters, outdir):
    '''
    simul_n: the number of random environments to generate.
    transporters: file containing list of trasporters (with modelSeed ids)
    outdir: directory to store results
    '''
    
    
    with open(transporters) as f:
        transporters=f.readline().strip().split('\t')  
        names = f.readline().strip().split('\t')
    water_index=transporters.index('EX_cpd00001_e')
    
    oxygen_index=transporters.index('EX_cpd00007_e')

    
    random_environments=np.zeros((simul_n,len(transporters)))
    
    
    
    for i in range(simul_n):

        #Draw from a Dirichlet distribution
        a=np.random.dirichlet(np.ones(len(transporters)))

        #Set water to a constant value
        a[water_index]=1.0

        #Set a 50% chance of being anoxic
        if np.random.uniform()<0.5:
            a[oxygen_index]=0.0
        else:
            a[oxygen_index]=1.0

        #stretch the distribution
        a=a*1000
        random_environments[i]=a
    
    with open(os.path.join(outdir,   'env_ball_'+str(simul_n)+'.tsv'),'w') as f:
        
        f.write('\t'.join(transporters) + '\n')
        f.write('\t'.join(names) + '\n')
        
        
        for i in random_environments:
            
            for num in i:
                f.write('%.20f\t'%num)
            f.write('\n')
    
#Example of usage:
#main(1000, pathToFileWithListOfTransporters, pathToWriteOutput)
