# code B > src > dementpy
# Elsa Abs, Feb 2020
"""
-------------------------------------------------------------------------------
      DEMENTpy--Decomposition Model of Enzymatic Traits in Python,v1.0
                              Bin Wang, Ph.D.
              Department of Ecology and Evolutionary Biology
                       University of California Irvine
                  Emails: wbwenwu@gmail.com or bwang7@uci.edu
                          Twitter: @bio_atmosphere
-------------------------------------------------------------------------------
"""
import os
import sys
import time
import pandas as pd
import numpy  as np

from preinitialization import preinitialize_data
from initialization import initialize_data
from grid import Grid
from output import Output
from utility import export

def main():
    
    
    print("""
    ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                               Version 1.0
               Department of Ecology and Evolutionary Biology
                     University of California Irvine
    ---------------------------------------------------------------------------       
    """)
      
    #...Obtain the command line arguments
    job_ID         = sys.argv[1]   # jobID
    input_folder   = sys.argv[2]   # input folder name
    output_folder  = sys.argv[3]   # output folder name
    final_litter   = sys.argv[4]   # phase 2 litter
    residents      = sys.argv[5]   # litter the residents come from
    invaders       = sys.argv[6]   # litter the invaders come from
    enzmax         = sys.argv[7]   # enzmax the residents and invaders come from
    seed_index     = sys.argv[8]   # seed_index the residents and invaders come from (seed_index 1 = 2102121235, seed_index 2 = seed_index + 1 = 2102121236)
    pulse          = sys.argv[9]   # number of years I'm running
    
    #...Get constants in runtime
    os.chdir('../'+input_folder)
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    
    #...Go get the right initial_substrates for the invasion simulation (can be different from both the invaders' and the resident's ones)
    os.chdir('./litters/'+final_litter)
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    
    #...Go get the resident grid that's going to be the residents, and the resident grid that I'm going to use to make the invaders
    os.chdir('../../grids')
    residents_grid = pd.read_csv('residents_grid_'+residents+seed_index+'.csv', header=0, index_col=0)
    if invaders != '0':
        invaders_pre_grid = pd.read_csv('residents_grid_'+invaders+seed_index+'.csv',  header=0, index_col=0)
    else:
        invaders_pre_grid = pd.read_csv('invaders_grid_0.csv', header=0, index_col=0)

    #...grow a seed of random number generator
    seed_number = 2102121234 + int(seed_index)
    np.random.seed(seed_number)

    mic_reinit = False    # indicate reinitialization of microbial community
    
    #...Set up the working directory
    os.chdir('../')
    
    #...start timing
    start_time = time.time()
    
    #...Preinitialize data by calling the Function: Preinitialize_Data()
    data_preinitialization = preinitialize_data()
    
    #...Initialize data by calling the Function: Initialize_Data()
    data_initialization = initialize_data(runtime,substrates_init,data_preinitialization,residents_grid,invaders_pre_grid)

    #...Prepare for output by creating an instance of the Output class
    Output_init = Output(runtime,data_initialization)

    #...Create an instance of the Grid class
    Ecosystem = Grid(runtime,data_initialization)

    #...Run the model
    for p in range(int(pulse)):
        
        for i in range(p*cycle, (p+1)*cycle):
        
            # substrates degradation
            Ecosystem.degradation(p,i)
        
            # monomers uptake
            Ecosystem.uptake(p,i)
        
            # microbial metabolism
            Ecosystem.metabolism(i)
            
            # microbial death
            Ecosystem.mortality(i)
        
            # microbial reproduction and dispersal
            Ecosystem.reproduction(i)
        
            # re-initialize microbial community in each new pulse
            if i == (p+1)*cycle-1:
                Ecosystem.repopulation(i)
        
            # output data using the "output" method in the Output class
            if i == 0:
                #Output_init.output(Ecosystem,i)  # day 1
                Output_init.output(Ecosystem,i,runtime)  # day 1
            elif i%interval==interval-1:
                #Output_init.output(Ecosystem,i)  # interval
                Output_init.output(Ecosystem,i,runtime)  # interval
    
    #...export the Output_init object to the output_folder using the export() funtion in the utility module 
    os.chdir('../'+output_folder)
    export(Output_init,job_ID)
    
    #print end_time
    #...Final microbial community
    #Outputs = Output(runtime, initialize_data)
    #Outputs = Output(runtime, data_init)
    #output = Outputs.ouput(Ecosystem,day)
    #print('microbes_C_fin = ', microbes_C_fin)
    
    #microbial_community_final_C = Microbes.microbial_community_initialization()
    #print('end time = ', cycle)
    
    #print out time used
    print('   ',"Run time:", time.time() - start_time, "seconds")
    
main()
