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
    job_ID            = sys.argv[1]   # jobID
    input_folder      = sys.argv[2]   # input folder name
    output_folder     = sys.argv[3]   # output folder name
    litter            = sys.argv[4]   # litter folder and name
    seed_index        = sys.argv[5]   # seed number
    enz_per_taxon_max = sys.argv[6]   # max number of enzymes per taxon
    pulse             = sys.argv[7]   # number of years simulated
    
    #...Go get the right initial_substrates
    os.chdir('../'+input_folder)
    os.chdir('./'+litter)
    substrates_init = pd.read_csv('initial_substrates.csv', header=0,    index_col=0).astype('float32')   # initial substrates
    
    #...Set up the working directory
    os.chdir('../')

    #...grow a seed of random number generator
    #np.random.seed(int(outname[:-4]))
    seed_number = 2102121234 + int(seed_index)
    np.random.seed(seed_number)
    #np.random.seed(int(seed_number[:]))

    #...a few system constants
    runtime    = pd.read_csv('runtime.txt',header=None,index_col=0,sep='\t')
    #pulse      = int(runtime.loc['pulse',1])         # number of pulses
    cycle      = int(runtime.loc['end_time',1])      # number of time steps in each pulse
    interval   = int(runtime.loc['interval',1])      # interval of time step to record outputs
    mic_reinit = False    # indicate reinitialization of microbial community
    
    #...start timing
    start_time = time.time()
    
    #...Initialize data by calling the Function: Initialize_Data()
    data_initialization = initialize_data(substrates_init,runtime)

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
    
    # prints
    print('litter =',litter)
    print('seed_index =',seed_index)
    print('seed_number =',seed_number)
    print('enz_per_taxon_max =', enz_per_taxon_max)
    print('simulation end time (number of years) =', pulse)
    
    #print out time used
    print('   ',"Run time:", time.time() - start_time, "seconds")
    
main()
