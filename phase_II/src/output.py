# output.py module dealing with outputs of DEMENTpy.
# Elsa Abs, Sep 2020
# Code Litter B - Seed S1
# Goal: Keep trait table + biomasse table

import numpy as np
import pandas as pd
import sys

class Output():
    
    def __init__(self,runtime,data_init):
        
        # Pass all runtime parameters to Runtime
        self.Runtime = runtime
        # Pass all initialized data (a dictionary) to 'Initialization'
        self.Initialization = data_init
        
        end_time     = int(runtime.loc['end_time',1])
        final_litter = sys.argv[4]
        residents      = sys.argv[5]
        invaders       = sys.argv[6]
        enzmax         = sys.argv[7]
        seed_index     = sys.argv[8]
        pulse        = int(sys.argv[9])
        seed_number = 2102121234 + int(seed_index)
        
        # Decay rate
        n_substrates              = int(runtime.loc['n_substrates',1])
        Substrates_grid           = data_init['Substrates'].groupby(level=0,sort=False).sum()
        Sub_index                 = Substrates_grid.index
        DecayRates_grid           = pd.DataFrame(data=[0]*n_substrates, index=Sub_index)
        DecayRates_grid.name      = 0  # set the series name to 0
        self.DecayRatesSeries     = DecayRates_grid

        # Total over grid [C]
        Microbes_grid           = data_init['Microbes'].groupby(level=0,sort=False).sum()  # Sum biomass per taxon
        Microbes_grid['C'].name = 0                                                    # Only keep the carbon mass, and rename the column name 0
        self.MicrobesSeries     = Microbes_grid['C']                                   # Taxon-specific total biomass summed over the grid
        self.TotalBiomass_PerYear = []          # Prepare list of total biomass over taxa and grid averaged per year
        
        print('invaders output',invaders)
        print('type(invaders)',type(invaders))
        print('type(residents+enzmax+seed_index)',type(residents+enzmax+seed_index))
        
        # Prepare the table to gather year data
        if invaders != '0':
            self.decay_results = pd.DataFrame([[final_litter+residents+invaders+enzmax+seed_index, final_litter, 'Res'+residents, 'Inv'+invaders, enzmax, seed_number]], columns=['treatment_name', 'final_litter', 'residents', 'invaders', 'enz_max', 'seed_number'])
        else:
            self.decay_results = pd.DataFrame([[final_litter+residents+invaders+enzmax+seed_index, final_litter, 'Res'+residents, 'NoInv', enzmax, seed_number]], columns=['treatment_name', 'final_litter', 'residents', 'invaders', 'enz_max', 'seed_number'])
        


       
    def output(self,ecosystem,day,runtime):
        
        end_time     = int(runtime.loc['end_time',1])
        final_litter = sys.argv[4]
        residents      = sys.argv[5]
        invaders       = sys.argv[6]
        enzmax         = sys.argv[7]
        seed_index     = sys.argv[8]
        pulse        = int(sys.argv[9])
        
        # Microbe
        # Total biomass per taxon per day
        Microbes_grid  = ecosystem.Microbes.groupby(level=0,sort=False).sum()
        Microbes_grid['C'].name = day + 1
        self.MicrobesSeries  = pd.concat([self.MicrobesSeries, Microbes_grid['C']], axis=1, sort=False)
        
        # Decay rate (averaged over the grid, but not over time, by substrate)
        DecayRates_grid           = ecosystem.DecayRates.groupby(level=0, sort=False).mean()
        DecayRates_grid.name      = day + 1     # index the output by day
        self.DecayRatesSeries     = pd.concat([self.DecayRatesSeries,DecayRates_grid], axis=1, sort=False)
        
        # Decay rate (averaged over grid, time)                                                                         (one value per substrate)
        self.DecayRates_Grid_Time   = self.DecayRatesSeries.mean(axis=1)
        # Decay rate (averaged over grid, time, substrate)                                                              (one value)
        self.DecayRates_Grid_Time_Substrate   = self.DecayRates_Grid_Time.mean()
        
        # One mean biomass and decay rate value per year
        for i in range(pulse):
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
        
        if day == pulse*end_time - 1:
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
            export_csv = decay_results_1job.to_csv(r'../output/decay_results_'+final_litter+residents+invaders+enzmax+seed_index+'.csv', header=True, index=True)
