# This module, enzyme.py, deals with exclusively only enzyme-related parameters, holding:
# ............Enzyme():           class    
# ............Boltzman_Arrhenius: function
# Bin Wang in Janunary, 2020 

import os
import sys
import csv
import pandas as pd
import numpy as np
from utility import LHS

class Enzyme():
    """
    This class deals with all properties associated with exoenzymes.
    
    Including methods:
    1. enzyme_pool_initialization():
    2. enzyme_attributes():
    3. enzyme_Ea():
    4. enzyme_uptake_Ea():
    5. enzyme_Vmax():
    6. enzyme_uptake_Vmax():
    7. enzyme_Km():
    8. enzyme_uptake_Km():
    """
    
    def __init__(self,runtime,parameters,data_preinit,substrate_index):
        """
        The constructor of Enzyme class.

        Parameters:
            data_preinit: dictionaries;preinitialized data from the module 'preinitialization.py'
            parameters: model parameters
        """
        
        self.substrate_index    = substrate_index                           # index of substrates in their actual names
        
        # Parameters from runtime
        self.n_substrates       = int(runtime.loc['n_substrates',1])
        self.n_monomers         = self.n_substrates + 2
        
        # Parameters from parameters.csv (this is why I need to import parameters.csv in initialization.py before calling Enzyme)
        self.Enz_min            = parameters.loc['Enz_min',1]               # Initial min. enzyme present in terms of carbon
        self.Enz_max            = parameters.loc['Enz_max',1]               # Initial max. enzyme present in terms of carbon
        
        # Import Data Dictionary residents and invaders from preinitialization.py
        data_preinit_residents = data_preinit[0]
        data_preinit_invaders = data_preinit[1]
        
        # Parameters from data_preinit for Degradation
        self.n_enzymes_res = data_preinit_residents['n_enzymes']
        self.n_enzymes_inv = data_preinit_invaders['n_enzymes']
        self.n_enzymes = self.n_enzymes_res + self.n_enzymes_inv
        self.Vmax0_residents = data_preinit_residents['Vmax0'][0:self.n_enzymes_res]
        self.Vmax0_invaders = data_preinit_invaders['Vmax0'][0:self.n_enzymes_inv]
        self.Km0_residents = data_preinit_residents['Km0'][0:self.n_substrates]
        self.Km0_invaders = data_preinit_invaders['Km0'][0:self.n_substrates]
        
        # Parameters from data_preinit for Uptake
        self.n_uptake_res = data_preinit_residents['n_uptake']
        self.n_uptake_inv = data_preinit_invaders['n_uptake']
        self.n_uptake = self.n_uptake_res + self.n_uptake_inv
        self.Uptake_Vmax0_residents = data_preinit_residents['Uptake_Vmax0'][0:self.n_monomers]
        self.Uptake_Vmax0_invaders = data_preinit_invaders['Uptake_Vmax0'][0:self.n_monomers]
        self.Uptake_Ea_residents = data_preinit_residents['Uptake_Ea'][0:self.n_monomers]
        self.Uptake_Ea_invaders = data_preinit_invaders['Uptake_Ea'][0:self.n_monomers]
        self.Uptake_Km0_residents = data_preinit_residents['Uptake_Km0'][0:self.n_monomers]
        self.Uptake_Km0_invaders = data_preinit_invaders['Uptake_Km0'][0:self.n_monomers]
        
        # Parameters from data_preinit for Metabolism
        self.EnzAttrib_residents = data_preinit_residents['EnzAttrib'][0:self.n_enzymes_res]
        self.EnzAttrib_invaders = data_preinit_invaders['EnzAttrib'][0:self.n_enzymes_inv]


    def enzyme_pool_initialization(self):
        """
        Initialize the pool sizes of different enzymes.

        Parameters:
            n_enzymes:  number of resident enzymes + number of invader enzymes (calculated in __init__)
            Enz_min:    min. of initial enzyme in C
            Enz_max:    max. of initial enzyme in C
        Return:
            Enzymes_df: series; index:Enz
        """

        Enzymes_array = np.random.uniform(self.Enz_min,self.Enz_max,self.n_enzymes)
        index = ['Enz'+str(i) for i in range(1,self.n_enzymes+1)]
        Enzymes_df = pd.Series(data=Enzymes_array, index=index, name='C', dtype='float32')
        n_enzymes = self.n_enzymes

        return n_enzymes, Enzymes_df
    
    
    def enzyme_Vmax(self):
        """
        Pre-exponential constants for enzymes.
        
        Parameters:
          n_enzymes_res:     number of resident enzymes; from Dictionary_Residents
          n_enzymes_inv:     number of invaders enzymes; from Dictionary_Invaders
          Vmax0_residents:   dataframe(enzyme*substrate); from Dictionary_Residents
          Vmax0_invaders:    dataframe(enzyme*substrate); from Dictionary_Invaders
        Returns:
          Vmax0:             dataframe(substrates*enzymes) of the invaded community; for Dictionary_Common

        """
        
        # Change index of Enz in Vmax0 invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        self.Vmax0_invaders.index = enzyme_inv_index
        
        # Concatenate Vmax0 residents and invaders
        Vmax0 = pd.concat([self.Vmax0_residents,self.Vmax0_invaders], sort=False)
        
        return Vmax0
        
        
    def enzyme_Ea(self,Ea_input):
        # I kept the same code as Bin's because degradation Ea is only substrate dependent (doesn't depend on enzyme)
        """
        Enzyme specificity matrix of activation energies.
       
        Parameter:
            Ea_input: dataframe; substrate-specific activation energy range (min = max for now)
        Return:
            Ea_df.T:  dataframe; Rows:enzymes; cols: substrates
        """
        
        Ea_series = Ea_input.apply(lambda df: np.random.uniform(df['Ea_min'],df['Ea_max'],self.n_enzymes),axis=1) # series of 1D array
        columns   = ['Enz' + str(i) for i in range(1,self.n_enzymes + 1)]
        Ea_df     = pd.DataFrame(data=Ea_series.tolist(), index=self.substrate_index, columns=columns, dtype='float32') # NOTE: .tolist()
        
        return Ea_df.T
     
     
    def enzyme_Km(self):
        """
        Derive Km through implementing a Vmax-Km tradeoff.

        The tradeoff is based on a linear correlation between Km and Vmax; minimum Km constrained to Km_min
        
        Parameters:
          Km_residents:   dataframe(substrate*enzyme); from Dictionary_Invaders
          Km_invaders:    dataframe(substrate*enzyme); from Dictionary_Invaders
        Returns:
          Km:             dataframe(substrate*enzyme); for Dictionary_Common
        """
        
        # Change index of Enz in Km invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        Km0_invaders_T = self.Km0_invaders.T
        Km0_invaders_T.index = enzyme_inv_index
        self.Km0_invaders = Km0_invaders_T.T
        
        # Concatenate Km) residents and invaders
        Km = pd.concat([self.Km0_residents,self.Km0_invaders], axis=1, sort=False)
        
        return Km
        
        
    def enzyme_uptake_Vmax(self):
        """
        Pre-exponential constants for uptake.
        
        Parameters:
            n_uptake:                  number of resident transporters + number of invader transporters (calculated in __init__)
            Uptake_Vmax0_residents:    dataframe(monomers*enzymes); from Dictionary_Residents
            Uptake_Vmax0_invaders: 1   dataframe(monomers*enzymes); from Dictionary_Invaders
        Return:
            Uptake_Vmax0:              dataframe(monomers*enzymes); for Dictionary_Common
        """
        
        # Change index of Upt in Uptake_Vmax0 invaders
        uptake_inv_index = []
        for i in np.arange(self.n_uptake_inv):
            uptake_inv_index.append("Upt" + str(self.n_uptake_res+i+1))
        Uptake_Vmax0_invaders_T = self.Uptake_Vmax0_invaders.T
        Uptake_Vmax0_invaders_T.index = uptake_inv_index
        self.Uptake_Vmax0_invaders = Uptake_Vmax0_invaders_T.T
        
        # Concatenate Uptake_Vmax0 residents and invaders
        Uptake_Vmax0 = pd.concat([self.Uptake_Vmax0_residents,self.Uptake_Vmax0_invaders], axis=1, sort=False)
        
        return Uptake_Vmax0
        
        
    def enzyme_uptake_Ea(self):
        """
        Uptake activation energies constant across monomers.
        
        Parameters:
            n_uptake:               number of resident transporters + number of invader transporters (calculated in __init__)
            Uptake_Ea_residents:    dataframe(monomers*transporters); from Dictionary_Residents
            Uptake_Ea_invaders:     dataframe(monomers*transporters); from Dictionary_Invaders
        Return:
            Uptake_Ea:              dataframe(monomers*transporters); for Dictionary_Common
        """
        
        # Change index of Upt in Uptake_Ea invaders
        uptake_inv_index = []
        for i in np.arange(self.n_uptake_inv):
            uptake_inv_index.append("Upt" + str(self.n_uptake_res+i+1))
        Uptake_Ea_invaders_T = self.Uptake_Ea_invaders.T
        Uptake_Ea_invaders_T.index = uptake_inv_index
        self.Uptake_Ea_invaders = Uptake_Ea_invaders_T.T
        
        # Concatenate Uptake_Vmax0 residents and invaders
        Uptake_Ea = pd.concat([self.Uptake_Ea_residents,self.Uptake_Ea_invaders], axis=1, sort=False)
        
        return Uptake_Ea
        
        
    def enzyme_uptake_Km(self):
        """
        Derive Uptake Km by implementing Vmax-Km tradeoff.

        The tradeoff is based on a linear correlation between Km and Vmax; minimum Uptake_Km constrained to Uptake_Km_min

        Parameters:
            n_uptake_res:           number of resident transporters
            n_uptake_inv:           number of invader transporters
            Uptake_Km0_residents:   dataframe(monomers*transporters); from Dictionary_Residents
            Uptake_Km0_invaders:    dataframe(monomers*transporters); from Dictionary_Invaders
        Return:
            Uptake_Km:              dataframe(monomers*transporters); for Dictionary_Common
        """
        
        # Change index of Upt in Uptake_Ea invaders
        uptake_inv_index = []
        for i in np.arange(self.n_uptake_inv):
            uptake_inv_index.append("Upt" + str(self.n_uptake_res+i+1))
        Uptake_Km0_invaders_T = self.Uptake_Km0_invaders.T
        Uptake_Km0_invaders_T.index = uptake_inv_index
        self.Uptake_Km0_invaders = Uptake_Km0_invaders_T.T
        
        # Concatenate Uptake_Vmax0 residents and invaders
        Uptake_Km = pd.concat([self.Uptake_Km0_residents,self.Uptake_Km0_invaders], axis=1, sort=False)
        
        return Uptake_Km
        
        
    def enzyme_attributes(self):
        """
        Derive enzyme stoichiometry and maintenence cost.

        Parameters:
            Enz_C_cost:     scalar;
            Enz_N_cost:     scalar;
            Enz_P_cost:     scalar;
            Enz_Maint_cost: scalar;
        Return:
            EnzAttrib_df: dataframe; row: Enz; col:["C_cost","N_cost","P_cost","Maint_cost"]
        """
        
        # Change index of Enz in EnzAttrib invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        self.EnzAttrib_invaders.index = enzyme_inv_index
        
        # Concatenate EnzAttrib residents and invaders
        EnzAttrib_df = pd.concat([self.EnzAttrib_residents,self.EnzAttrib_invaders], axis=0, sort=False)
        
        return EnzAttrib_df


def Arrhenius(Ea,temperature):
    """
    Temperature dependence of rate constant. 

    Parameters:
       Ea:          dataframe/scalar; activation energy;
       temperature: scalar;           daily temperature; 
    Return:
       BA:          dataframe/scalar; dependent on Ea
    Reference:
        Wikipedia: https://en.wikipedia.org/wiki/Arrhenius_equation
    """

    Tref = 293
    k = np.exp((-Ea/0.008314)*(1/(temperature+273) - 1/Tref)).astype('float32')

    return k
