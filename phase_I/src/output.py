# output.py module dealing with outputs of DEMENTpy.
# Elsa Abs, April 2022
# Goal: Get series of overflow of C,N,P due to stoichiometric unbalance per degrader type

import numpy as np
import pandas as pd
import sys
import csv
import re
from utility import isFloat
from utility import isFloatorInteger
import os
from os import path

class Output():
    """
    This class deals with outputs.
    
    Accepts data derived from the initialization.py and grid.py modules, and have two methods:
        output():             stores all time series
        microbes_abundance(): a special method
    """
    
    def __init__(self,runtime,data_init):
        
        # Pass all runtime parameters to Runtime
        self.Runtime = runtime
        # Pass all initialized data (a dictionary) to 'Initialization'
        self.Initialization = data_init
        
        job_ID            = sys.argv[1]   # jobID
        input_folder      = sys.argv[2]   # input folder name
        output_folder     = sys.argv[3]   # output folder name
        litter            = sys.argv[4]   # litter folder and name
        seed_index        = sys.argv[5]   # seed number
        enzmax = sys.argv[6]   # max number of enzymes per taxon
        pulse             = sys.argv[7]   # number of years simulated
        seed_number       = 2102121234 + int(seed_index)
        
        # DICTIONARIES
        if not os.path.exists(r'../output/dic'+litter+seed_index):
            os.mkdir(r'../output/dic'+litter+seed_index)
            os. chdir(r'../output')
            # Export data_init
            A = data_init
            # Transform data_init into csv files
            for key, val in A.items():
                if isinstance(val, pd.DataFrame):
                    val.to_csv(r'./dic'+litter+seed_index+'/data_{}.csv'.format(str(key)))
                elif isinstance(val, pd.Series):
                    val.to_csv(r'./dic'+litter+seed_index+'/data_{}.csv'.format(str(key)), header=True)
                elif isinstance(val, np.ndarray):
                    np.savetxt(r'./dic'+litter+seed_index+'/data_{}.csv'.format(str(key)), val)
                else :
                    with open(r'./dic'+litter+seed_index+'/data_{}.csv'.format(str(key)), 'w') as f:      # this is an integer or a float
                        f.write('%s'%val)              # the title of the table inside the file is "da" because the table is small
            with open(r'./dic'+litter+seed_index+'/keys.txt', mode="w") as f:        #write keys in file
                f.write(str(list(A.keys())))
        
        # DECAY RATES
        n_substrates              = int(runtime.loc['n_substrates',1])
        ## Substrate
        Substrates_grid           = data_init['Substrates'].groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = 0
        self.SubstratesSeries     = Substrates_grid['C']
        ## Decay
        Sub_index                 = Substrates_grid.index
        DecayRates_grid           = pd.DataFrame(data=[0]*n_substrates, index=Sub_index)
        DecayRates_grid.name      = 0  # set the series name to 0
        self.DecayRatesSeries     = DecayRates_grid
        # Total over grid [C]
        Microbes_grid           = data_init['Microbes'].groupby(level=0,sort=False).sum()  # Sum biomass per taxon
        Microbes_grid['C'].name = 0                                                    # Only keep the carbon mass, and rename the column name 0
        self.MicrobesSeries     = Microbes_grid['C']                                   # Taxon-specific total biomass summed over the grid
        self.TotalBiomass_PerYear = []          # Prepare list of total biomass over taxa and grid averaged per year
        # Prepare the table to gather year data
        self.decay_results = pd.DataFrame([[litter+seed_index, litter, enzmax, seed_number]], columns=['treatment_name', 'litter', 'enz_max', 'seed_number'])
        self.decay_results_firstcolumns = pd.DataFrame([[litter+seed_index, litter, enzmax, seed_number]], columns=['treatment_name', 'litter', 'enz_max', 'seed_number'])
          
            
       
    def output(self,ecosystem,day,runtime):
        """
        Records outputs in various variables of each iteration.

        Parameters:
            ecosystem: object from the grid.py module
            day:       the day to record to outputs
        Returns:
        """
        
        end_time = int(runtime.loc['end_time',1])
        
        # Arguments
        job_ID            = sys.argv[1]   # jobID
        input_folder      = sys.argv[2]   # input folder name
        output_folder     = sys.argv[3]   # output folder name
        litter            = sys.argv[4]   # litter folder and name
        seed_index        = sys.argv[5]   # seed number
        enzmax = sys.argv[6]   # max number of enzymes per taxon
        pulse             = sys.argv[7]   # number of years simulated
        seed_number       = 2102121234 + int(seed_index)
        
        
        # GRIDS
        if day == int(pulse)*end_time - 1:
            litter     = sys.argv[4]
            seed_index = sys.argv[5]
            n_taxa     = int(runtime.loc['n_taxa',1])
            gridsize   = int(runtime.loc['gridsize',1])
            export_csv = ecosystem.Microbes.to_csv (r'../output/residents_grid_'+litter+seed_index+'.csv', header=True, index=True)
        
        
        # DECAY RATES
        # Microbe
        # Total biomass per taxon per day
        Microbes_grid  = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        Microbes_grid['C'].name = day + 1
        self.MicrobesSeries  = pd.concat([self.MicrobesSeries, Microbes_grid['C']], axis=1, sort=False)
        # Substrate
        Substrates_grid  = ecosystem.Substrates.groupby(level=0,sort=False).sum()
        Substrates_grid['C'].name = day + 1
        self.SubstratesSeries  = pd.concat([self.SubstratesSeries, Substrates_grid['C']], axis=1, sort=False)
        # Decay rate (averaged over the grid, but not over time, by substrate)
        DecayRates_grid           = ecosystem.DecayRates.groupby(level=0, sort=False).mean()
        DecayRates_grid.name      = day + 1     # index the output by day
        self.DecayRatesSeries     = pd.concat([self.DecayRatesSeries,DecayRates_grid], axis=1, sort=False)
        # Decay rate (averaged over grid, time)                                                                         (one value per substrate)
        self.DecayRates_Grid_Time   = self.DecayRatesSeries.mean(axis=1)
        # Decay rate (averaged over grid, time, substrate)                                                              (one value)
        self.DecayRates_Grid_Time_Substrate   = self.DecayRates_Grid_Time.mean()
        # PER YEAR
        for i in range(int(pulse)):
            if day == (i+1)*end_time - 1:
                # Total biomass
                Biomass_Total_Grid = self.MicrobesSeries.sum(axis=0)
                Biomass_Total_Grid_Taxa_Average_Year = Biomass_Total_Grid[i*end_time+1:(i+1)*end_time+1].mean()
                Biomass_df = pd.DataFrame([[Biomass_Total_Grid_Taxa_Average_Year]], columns=['biomass_year_' + str(i+1)])
                # Decay rate
                DecayRate_Average_Substrates = self.DecayRatesSeries.mean(axis=0)               # Average over the 12 substrates, i.e. Sum of all decay rates /12
                DecayRate_Average_Year = DecayRate_Average_Substrates[i*end_time+1:(i+1)*end_time+1].mean()
                DecayRate_df = pd.DataFrame([[DecayRate_Average_Year]], columns=['decay_year_' + str(i+1)])
                # Result dataframe
                self.decay_results = pd.concat([self.decay_results,Biomass_df,DecayRate_df], axis=1, sort=False)
        # AVERAGE OVER ALL YEARS
        if day == int(pulse)*end_time - 1:
            # Calculate average over all years
            ## Biomass
            Biomass_Total_Grid = self.MicrobesSeries.sum(axis=0)
            Biomass_Average_Time = Biomass_Total_Grid[1:].mean()
            Biomass_df = pd.DataFrame([[Biomass_Average_Time]], columns=['biomass'])
            # Decay
            DecayRate_Average_Substrates = self.DecayRatesSeries.mean(axis=0)
            Decay_Average_Time = DecayRate_Average_Substrates[1:].mean()
            DecayRate_df = pd.DataFrame([[Decay_Average_Time]], columns=['decay'])
            # Add both to dataframe
            self.decay_results = pd.concat([self.decay_results,Biomass_df,DecayRate_df], axis=1, sort=False)
            # Export to csv in output folder
            decay_results_1job = self.decay_results.copy(deep=True)
            export_csv = decay_results_1job.to_csv(r'../output/decay_results_'+litter+seed_index+'.csv', header=True, index=True)
