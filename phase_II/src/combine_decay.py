# combine_outputs.py module dealing with outputs of DEMENTpy.
# Elsa Abs, May 2021
# Goal: Open all phase2_results.csv and combine them into 1

import os
import time
import pandas as pd
import glob

def main():


    print("""
    ---------------------------------------------------------------------------
         DEMENTpy (DEcomposition Model of Enzymatic Traits in Python)
                               Version 1.0
               Department of Ecology and Evolutionary Biology
                     University of California Irvine
    ---------------------------------------------------------------------------
    """)
    
    #...start timing
    start_time = time.time()
    
    path = r'../output' # use your path
    all_files = glob.glob(path + '/decay_results_*.csv')

    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)
    export_csv = frame.to_csv(r'../output/decay_results_combined.csv', header=True, index=True)
    
    #print out time used
    print('   ',"Run time:", time.time() - start_time, "seconds")

main()
