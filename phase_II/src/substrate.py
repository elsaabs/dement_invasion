# substrate.py module holding a class, Substrate()
# by Bin Wang on Dec. 26th, 2019ls

import sys
import pandas as pd
import numpy as np


class Substrate():
    """
    Deals with all properties related to substrates.
    
    Includs methods of:
        1) substrate_input(): input rates of varying substrates to the system
        2) substrate_produced_monomer():  substrate-produced monomers
        3) substrate_degradation_enzyme(): substrate-required enzymes
    """
    
    def __init__(self,runtime,substrates_init,data_preinit):
        """
        The constructor Substrate class.

        Parameters:
            runtime:    dataframe; user-specified parameters when running the model
            parameters: dataframe; parameters related to microbes, substrates,enzymes, and monomers
        """
        
        # Parameters from runtime
        self.n_substrates      = int(runtime.loc['n_substrates',1])
        
        # Parameters from parameters.csv
        self.Substrates_start  = substrates_init.astype('float32')          # initial substrate concentrations
        
        # Parameters from data_preinit
        data_preinit_residents = data_preinit[0]
        data_preinit_invaders = data_preinit[1]
        self.n_enzymes_res = data_preinit_residents['n_enzymes']
        self.n_enzymes_inv = data_preinit_invaders['n_enzymes']
        self.n_enzymes = self.n_enzymes_res + self.n_enzymes_inv
        self.ReqEnz_residents = data_preinit_residents['ReqEnz']
        self.ReqEnz_invaders = data_preinit_invaders['ReqEnz']
    
    
    def substrate_input(self,sub_mon_input):
        """
        Substrate inputs during simulation.

        Parameter:
            sub_mon_input: dataframe; user-provided substrates and monomers input data
        Return:
            SubInput_df: dataframe; datatype: float32
        """

        
        # Substrate input rates
        SubInputC = sub_mon_input['Sub'] # access only the Substrate column
        SubInputN = SubInputC * self.Substrates_start["N"]/self.Substrates_start["C"]
        SubInputP = SubInputC * self.Substrates_start["P"]/self.Substrates_start["C"]
        # Rename series name
        SubInputC.name = 'C' 
        SubInputN.name = 'N'
        SubInputP.name = 'P'
        SubInput_df = pd.concat([SubInputC,SubInputN,SubInputP],axis=1,sort=False)
        SubInput_df['DeadMic'] = SubInput_df['DeadEnz'] = 0  # Change NAs to 0
        SubInput_df = SubInput_df.astype('float32')

        return SubInput_df
        
  
    def substrate_produced_monomer(self):
        """
        Monomers produced by each substrate.

        Return:
            MonomersProduced_df: dataframe; row:substrate; col:monomers; all rows should sum to 1
        """
        
        MonomersProduced_array = np.concatenate((np.array([0]*self.n_substrates*2).reshape((self.n_substrates,2),order='F'),np.diagflat([1]*self.n_substrates)),axis=1)
        index   = ['Sub'+str(i) for i in range(1,self.n_substrates+1)]
        columns = ['Mon'+str(i) for i in range(-1,self.n_substrates+1)]
        MonomersProduced_df = pd.DataFrame(data=MonomersProduced_array,index=index,columns=columns,dtype='int8')
        MonomersProduced_df.rename(columns = {'Mon-1':"NH4",'Mon0':"PO4",'Mon1':"DeadMic",'Mon2':"DeadEnz"},inplace=True)

        return MonomersProduced_df
        
    
    def substrate_degradation_enzyme(self):
        # I'm not doing the cleanest way here, which would be to take [0:n_substrates] of each set, then to multiply by gridsize each set, and to recombine the 2 sets
        # This will be needed if I change the gridsize between before and after invasion
        # Here I'm only copying ReqEnz from residents and invaders and concatenating them
        """
        Derive the required enzymes of each substrate.

        Return:
            ReqEnz_df: 3-D DataFrame: sets (2) * substrates * enzymes
        .............................................................
        Structure illustration:
        .............................................................
                   Enz1 Enz2 ... Enzn
        Set 1 Sub1 ------------------
        Set 1 Sub1 ------------------
        .     .    .     .       .
        .     .    .     .       .
        Set 1 Subn ------------------
        Set 2 Sub1 ------------------
        Set 2 Sub1 ------------------
        .     .    .     .       .
        .     .    .     .       .
        Set 2 Subn ------------------
        .............................................................
        """

        # Change index of ReqEnz in ReqEnz invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        ReqEnz_invaders_T = self.ReqEnz_invaders.T
        ReqEnz_invaders_T.index = enzyme_inv_index
        self.ReqEnz_invaders = ReqEnz_invaders_T.T
        
        # Concatenate Km) residents and invaders
        ReqEnz_df = pd.concat([self.ReqEnz_residents,self.ReqEnz_invaders], axis=1, sort=False)
        

        return ReqEnz_df
