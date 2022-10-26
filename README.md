# dement_invasion

DEMENTpy_invasion code and Google Colaboratory notebook for all figures (manuscript and supplementary)

## DEPENTpy_invasion 

The invasion code is made of 2 phases (phase I and phse II). Both only differ by the initialization. In phase I, microbial communities are built randomly. In phase II, we use communities (trait combination, biomass, spatial distriution on grid) coming from end of phase I (thanks to the preinitialization.py python file), and we can add or not invaders (through the bash command, see below).

How to run it:
1. Run phase I with bash command "for i in {A,B,C,D,E}; do export i; sbatch script.sh; done"
2. Once phase I is over, move (manually using cyberduck or through terminal with command "mv") phaseI>output>dics to phaseII>input>dictionaries of phase II, and all phaseI>output>residents* to phaseII>input>grids. This is for the initialization of phase II. 
3. Run phase II with bash command "for i in {A,B,C,D,E}; do for j in {A,B,C,D,E}; do for k in 0; do export i j k; sbatch script.sh; done; done; done".

You'll find the script.sh in this repository. They contain mostly 2 arguments: the number of years of simulation run (here 3), and the number of replicates (with job-arrray, here 10). 

## Google Colaboratory Notebook

You will need to import the csv output files (e.g. decay_results, resident_grids) in a google drive, and to replace the path in the Google Colaboratory Notebook that corresponds to where your output files are.
