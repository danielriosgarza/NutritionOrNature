# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 03:32:21 2020

@author: danie
"""


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import cobra
import numpy as np
from Env_ball_class import Env_ball
import gc
import pickle

#family=sys.argv[1]

def apply_environment(mdl, env_vec,transporters):
    for i in range(len(transporters)):
        mdl.reactions.get_by_id(transporters[i]).lower_bound=-env_vec[i]
        mdl.reactions.get_by_id(transporters[i]).upper_bound=1000.

    sol_m=mdl.optimize()
    return sol_m.objective_value


def remove_reaction(orig_obj, model, reaction, up_low_b_dict):

    model.reactions.get_by_id(reaction).upper_bound=0
    model.reactions.get_by_id(reaction).lower_bound=0



    new_obj=np.round(float(model.slim_optimize()),decimals=12)


    if new_obj<orig_obj:
        model.reactions.get_by_id(reaction).upper_bound=up_low_b_dict[reaction][0]
        model.reactions.get_by_id(reaction).lower_bound=up_low_b_dict[reaction][1]


def add_reaction(orig_obj, model, reaction, up_low_b_dict):

    model.reactions.get_by_id(reaction).upper_bound=up_low_b_dict[reaction][0]
    model.reactions.get_by_id(reaction).lower_bound=up_low_b_dict[reaction][1]



def get_model_sample (niter, del_p, model, up_low_b_dict, reactions, transporters, environment):
    apply_environment(model, environment, transporters)
    samples = np.random.randint(0, len(reactions), niter)
    orig_obj = .1*np.round(float(model.slim_optimize()),decimals=12)
    for i in range(niter):
        if np.random.binomial(1, del_p):
            remove_reaction(orig_obj, model, reactions[samples[i]], up_low_b_dict)
        else:
            add_reaction(orig_obj, model, reactions[samples[i]], up_low_b_dict)
    print(orig_obj, np.round(float(model.slim_optimize()),decimals=12))
    
    prof_r = np.ones(len(reactions))
    for i,v in enumerate(reactions):
        if (model.reactions.get_by_id(v).upper_bound==0) and (model.reactions.get_by_id(v).lower_bound==0):
            prof_r[i] = 0.0
    prof_t = np.ones(len(transporters))
    for i,v in enumerate(transporters):
        if model.reactions.get_by_id(v).flux>=0:
            prof_t[i] = 0.0
    return prof_r, prof_t
            

eb=Env_ball(1000)
transporters=eb.transporters[:]
random_environments=eb.matrix.copy()

model = cobra.io.read_sbml_model('C:/Users/danie/Documents/random_work_stuff/home/home/Files/Aeromonadaceae.ensembl.sbml')

model.solver='gurobi'
reactions=[i.id for i in model.reactions if 'rxn' in i.id]

up_low_b_dict = {}
for reaction in reactions:
    up_low_b_dict[reaction] = ( model.reactions.get_by_id(reaction).upper_bound,  model.reactions.get_by_id(reaction).lower_bound)


prof_r = np.zeros(len(reactions))
prof_t  = np.zeros(len(transporters))


for i in range(1000):
    env_choice = np.random.randint(0,1000)
    print(i)
    mc=model.copy()
    p_r, p_t = get_model_sample(10000, .8, mc, up_low_b_dict, reactions, transporters, random_environments[env_choice])
    prof_r += (p_r/1000)
    prof_t += (p_t/1000)
    gc.collect()

prof_rb, prof_tb = get_model_sample(10000, .8, model, up_low_b_dict, reactions, transporters, random_environments[10])

#print(prof_r)
#print(prof_t)

store = {'reactions':prof_r.copy(), 'metabolites':prof_t.copy(), 'r':reactions, 'm':transporters}
pickle.dump(store, open('C:/Users/danie/Documents/random_work_stuff/home/home/Files/AeromonadaceaeRand.pkl', 'wb'))
