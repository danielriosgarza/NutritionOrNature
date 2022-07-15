# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:43:03 2022

@author: u0139894
"""


from pathlib import Path
import pickle
import os


data_folder = os.path.join(Path(os.getcwd()).parents[0], 'Files', 'data')

store =  pickle.load(open(data_folder + '/pickles/all_fams.tertiaryDS.pkl', 'rb'))

store.pop('Bacillaceae_plus_Anaplasmataceae', None)

print(store.keys())

print(store['Enterobacteriaceae'].keys())

print(store['Enterobacteriaceae']['y_metabolome'])
print(store['Enterobacteriaceae']['x_reactome'])