# code B > src > initialization
# Elsa Abs, Feb 2020
"""
This module, with only one function "initialized_data(), initializes data related to 
substrate, monomer, enzyme, and microbe, as well as their distribution on the spatial grid,
preceding the actual decompostion-related computations.
"""

import os
import sys
import pandas as pd
import numpy as np
import csv

from utility   import isFloat
from utility   import isFloatorInteger



def preinitialize_data():
    """
    This class deals with importing the Data dictionary from codes A (residents) and B (invaders)
    """
    
    residents      = sys.argv[5]
    invaders       = sys.argv[6]
    enzmax         = sys.argv[7]
    seed_index     = sys.argv[8]

    # Resident dictionary
    # Import all the csv files and convert to dictionary
    with open(r'./dictionaries/dic'+residents+seed_index+'/keys.txt', mode="r",) as f:       #read data from keys
        keys = eval(f.read())
    B = {}
    for key in keys:
        val = pd.read_csv('./dictionaries/dic'+residents+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=0)
        if len(val) == 0 :                                                   # it's an integer or a float (len(val)=0)
            with open('./dictionaries/dic'+residents+seed_index+'/data_{}.csv'.format(str(key)),'r') as f:
                for line in f:
                    if isFloat(line) == True:                               # it's a float
                        B[key] = np.float32(line)
                    else :                                                  # it's an integer
                        B[key] = np.int16(line)
        elif len(val.shape) == 1 :                                          # it's a series (len(val.shape)=1)
            B[key] = val.astype('float32')
        elif val.shape[1] == 0 :                                            # it's an array (len(val.shape)=2 but val.shape[1]=0)
            B[key] = np.loadtxt('./dictionaries/dic'+residents+seed_index+'/data_{}.csv'.format(str(key)), dtype='float32')
        elif isFloatorInteger(val.iloc[0,0]) == False:                      # it's a multiindex dataframe (val.iloc[0,0]=not a float)
            #print(key, 'is a multiindex dataframe')
            val = pd.read_csv('./dictionaries/dic'+residents+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=[0,1])
            B[key] = val.astype('float32')
        else :                                                              # it's a 1 index dataframe
            #print(key, 'is a 1 index dataframe')
            val = pd.read_csv('./dictionaries/dic'+residents+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=0)
            B[key] = val.astype('float32')
    Data_Dictionary_Residents = B
     
     
    # Invader dictionary
    if invaders == '0':
        Data_Dictionary_Invaders = Data_Dictionary_Residents
    else:
        # Import all the csv files and convert to dictionary
        with open(r'./dictionaries/dic'+invaders+seed_index+'/keys.txt', mode="r",) as f:       #read data from keys
            keys = eval(f.read())
        D = {}
        for key in keys:
            val = pd.read_csv('./dictionaries/dic'+invaders+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=0)
            if len(val) == 0 :                                                   # it's an integer or a float (len(val)=0)
                with open('./dictionaries/dic'+invaders+seed_index+'/data_{}.csv'.format(str(key)),'r') as f:
                    for line in f:
                        if isFloat(line) == True:                               # it's a float
                            D[key] = np.float32(line)
                        else :                                                  # it's an integer
                            D[key] = np.int16(line)
            elif len(val.shape) == 1 :                                          # it's a series (len(val.shape)=1)
                D[key] = val.astype('float32')
            elif val.shape[1] == 0 :                                            # it's an array (len(val.shape)=2 but val.shape[1]=0)
                D[key] = np.loadtxt('./dictionaries/dic'+invaders+seed_index+'/data_{}.csv'.format(str(key)), dtype='float32')
            elif isFloatorInteger(val.iloc[0,0]) == False:                      # it's a multiindex dataframe (val.iloc[0,0]=not a float)
                #print(key, 'is a multiindex dataframe')
                val = pd.read_csv('./dictionaries/dic'+invaders+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=[0,1])
                D[key] = val.astype('float32')
            else :                                                              # it's a 1 index dataframe
                #print(key, 'is a 1 index dataframe')
                val = pd.read_csv('./dictionaries/dic'+invaders+seed_index+'/data_{}.csv'.format(str(key)), squeeze=True, index_col=0)
                D[key] = val.astype('float32')
        Data_Dictionary_Invaders = D

    print('Initial number of residents taxa = ', Data_Dictionary_Residents['n_taxa'])
    print('Initial number of invaders taxa = ', Data_Dictionary_Invaders['n_taxa'])
            
            
            
    return Data_Dictionary_Residents, Data_Dictionary_Invaders
