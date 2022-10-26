# code B > src > initialization
# Elsa Abs, Feb 2020
"""
This module, with only one function "initialized_data(), initializes data related to 
substrate, monomer, enzyme, and microbe, as well as their distribution on the spatial grid,
preceding the actual decompostion-related computations.
"""

import pandas as pd
import numpy as np

from substrate import Substrate
from monomer   import Monomer
from enzyme    import Enzyme
from microbe   import Microbe
from utility   import expand



def initialize_data(runtime_parameters,substrates_init,data_preinit,residents_grid,invaders_pre_grid):
    """
    Parameters:
        runtime_parameters: user-specified parameters setting up the system;
                            all other paras loaded by reading the parameters.csv
    Return:
        Data_Dictionary: a dictionary of all variables that feeds the grid.py module
    """
    
    # Load all input files
    parameters      = pd.read_csv('parameters.csv',         header=None, index_col=0).astype('float32')   # parameters
    Ea_input        = pd.read_csv("enzyme_ea.csv",          header=0,    index_col=0).astype('float32')   # enzyme activation
    
    # daily temperature and water potential
    daily_temp = 15.7
    daily_psi  = -0.01

    #...an instance of Substrate class
    Substrates = Substrate(runtime_parameters,substrates_init,data_preinit)
    #...substrate initial pool size
    substrates_initial_pool = Substrates.Substrates_start
    #...substrates degradation required enzymes
    substrates_req_enzyme = Substrates.substrate_degradation_enzyme()
    
    #...an instance of Monomer class
    Monomers = Monomer(runtime_parameters,parameters,data_preinit)
    #...monomers initial pool size
    monomers_initial_pool = Monomers.monomer_initialization(substrates_initial_pool)
    #...initial monomer ratios
    monomer_ratio_inital = Monomers.monomer_ratios(monomers_initial_pool)
    #...monomers uptake required enzymes
    monomers_uptake_reqenzyme = Monomers.monomer_uptake_reqenzyme()
    
    #...an instance of Enzyme class
    Enzymes = Enzyme(runtime_parameters,parameters,data_preinit,substrates_initial_pool.index)
    #...enzyme initial pool size:0
    enzymes_initial_pool = Enzymes.enzyme_pool_initialization()
    #...enzymes of substrate degradation Vmax
    enzymes_Vmax = Enzymes.enzyme_Vmax()
    #...enzymes of substrate degradation Ea
    enzymes_Ea = Enzymes.enzyme_Ea(Ea_input)
    #...enzymes of substrate degradation Km
    enzymes_Km = Enzymes.enzyme_Km()
    #...monomers uptake enzyme Vmax
    enzymes_uptake_Vmax= Enzymes.enzyme_uptake_Vmax()
    #...monomers uptake enzyme Ea
    enzymes_uptake_Ea = Enzymes.enzyme_uptake_Ea()
    #...monomers uptake enzyme Km
    enzymes_uptake_Km = Enzymes.enzyme_uptake_Km()
    #...enzyme attributes
    enzymes_attributes = Enzymes.enzyme_attributes()
    
    #...an instance of Microbe class
    Microbes = Microbe(runtime_parameters,data_preinit,residents_grid,invaders_pre_grid)
    #...Microbial enzyme genes
    microbial_enzyme_gene = Microbes.microbe_enzyme_gene()
    #...Microbial transporter and monomer per taxa
    microbial_uptake_monomer = Microbes.microbe_uptake_monomer()
    #...Microbial uptake cost
    microbial_uptake_cost = Microbes.microbe_uptake_cost()
    #...Below is "MicrobeS" (with an S) because I defined Microbe class as such just above. And I imported Microbe class from microbe.py line 15.
    microbial_community = Microbes.microbial_community_initialization()
    #...Microbial osmolyte productoin rate
    microbial_osmolyte_prod_rate = Microbes.microbe_osmoproduction_rate()
    #...Microbial enzyme production rate
    microbial_enzyme_prod_rate = Microbes.microbe_enzproduction_rate()
    #...Microbial minimum ratios
    microbial_min_ratios = Microbes.minimum_cell_quota()
    #...Microbial mortality
    microbial_mortality = Microbes.microbe_mortality()
    #...Microbial drought tolerance
    microbial_drought_tol = Microbes.microbe_drought_tol()
    
    # Parameters from runtime
    n_substrates = int(runtime_parameters.loc['n_substrates',1])
    gridsize = int(runtime_parameters.loc['gridsize',1])
    x = int(runtime_parameters.loc['x',1])
    y = int(runtime_parameters.loc['y',1])

    # Kept the same order as Bin (doesn't follow the order of call in grid.py. I start following appearance in grid starting TaxUpt)
    Data_Dictionary_Common = {"n_taxa": microbial_community[0],                               # tuple[0]: number of resident + invader taxa (regardless of biomass=0)
                       "n_substrates": n_substrates,                                          # number of substrates
                       "n_enzymes": enzymes_initial_pool[0],                                  # number of resident + invader enzymes
                       "gridsize": gridsize,                                                  # number of grid boxes
                       "grid_x": x,                                                           # length of grid
                       "grid_y": y,                                                           # width of grid
                       "Substrates": expand(substrates_initial_pool,gridsize),                # spatial distribution of substrates at the beginning of invasion
                       "ReqEnz": substrates_req_enzyme,                                       #
                       "Monomers": expand(monomers_initial_pool,gridsize),                    # spatial distribution of monomers at the beginning of invasion (zero)
                       "Monomer_ratio":expand(monomer_ratio_inital,gridsize),                 #
                       "Uptake_ReqEnz":expand(monomers_uptake_reqenzyme,gridsize),            #
                       "Enzymes": expand(enzymes_initial_pool[1],gridsize),                   # spatial distribution of enzymes at the beginning of invasion (zero)
                       "Km0": expand(enzymes_Km,gridsize),                                    # enzyme half-saturation constant
                       "Uptake_Km0": expand(enzymes_uptake_Km,gridsize),                      # transporter half-saturation constant
                       "Uptake_Ea": expand(enzymes_uptake_Ea,gridsize),                       # transporter acitivation energy
                       "Uptake_Vmax0": expand(enzymes_uptake_Vmax,gridsize),                  # transporter reaction rate
                       "Ea": expand(enzymes_Ea,gridsize),                                     # enzyme activation energy
                       "Vmax0": expand(enzymes_Vmax,gridsize),                                # enzyme reaction rate
                       "Microbes": microbial_community[1],                                    # tuple[1]: spatial microbes at the beginning of invasion
                       "EnzGenes": microbial_enzyme_gene,                                     # enzymes each taxon has
                       "TaxMon_final":microbial_uptake_monomer[1],                               # all monomers that each taxon can take up
                       "UptakeGenes_trait": expand(microbial_uptake_cost[0],gridsize),        # single gene cost of transporter
                       "UptakeGenesCost": expand(microbial_uptake_cost[1],gridsize),          # distribution of transporter gene cost across taxa
                       "TaxUpt": expand(microbial_uptake_monomer[0],gridsize),                   # transporter that each taxon has, expanded on grid
                       'alpha': parameters.loc['alpha',1],                                    # factor delineating curve concavity of microbial response to drought
                       'wp_fc': parameters.loc['wp_fc',1],                                    # threshold below which microbes start to respond to drought
                       'wp_th': parameters.loc['wp_th',1],                                    # threshold below which microbes in full swing to respond to drought
                       'Uptake_Maint_cost': parameters.loc['Uptake_Maint_cost',1],            # constant of transporter maintenence cost
                       "OsmoProdConsti": expand(microbial_osmolyte_prod_rate[0],gridsize),    # distribution of consti. osmolyte gene cost across taxa
                       "EnzProdConstit": expand(microbial_enzyme_prod_rate[0],gridsize),      # distribution of consti. enzyme gene cost across taxa
                       "EnzAttrib": enzymes_attributes,                                       # enzyme stoichiometry and energy cost
                       "AE_ref": parameters.loc["CUE_ref",1],                                 # Reference assimilation efficiency: 0.5
                       "AE_temp": parameters.loc["CUE_temp",1],                               # AE temperature sensitivity; default: -0.016
                       "OsmoProdInduci": expand(microbial_osmolyte_prod_rate[1],gridsize),    # distribution of induci. osmolyte gene cost across taxa
                       "EnzProdInduce": expand(microbial_enzyme_prod_rate[1],gridsize),       # distribution of induci. enzyme gene cost across taxa
                       "MinRatios": expand(microbial_min_ratios,gridsize),                    # microbial cell min. ratios
                       'C_min': parameters.loc['C_min',1],                                    # C threshold of cell lysis
                       'N_min': parameters.loc['N_min',1],                                    # N threshold of cell lysis
                       'P_min': parameters.loc['P_min',1],                                    # P threshold of cell lysis
                       'basal_death_prob': np.tile(microbial_mortality[0],gridsize),          # basal death probability
                       'death_rate': np.tile(microbial_mortality[1],gridsize),                # change rate of death prob. agaist mositure
                       "TaxDroughtTol": expand(microbial_drought_tol,gridsize),               # distribution of taxon-specific drought tolerance
                       "fb": np.tile(microbial_community[2],gridsize),                        # tuple[3]: fungi index
                       'max_size_b': parameters.loc['max_size_b',1],                          # C quota threshold for bacterial cell division
                       'max_size_f': parameters.loc['max_size_f',1],                          # C quota threshold for fungal cell division
                      }
    
    #print('fb = ', Data_Dictionary_Common['fb'])
    #print('fb type = ', type(Data_Dictionary_Common['fb']))
    
    #print('Km0 = ', Data_Dictionary_Common['Km0']*10000)
    
    print('Gridsize = ', Data_Dictionary_Common['gridsize'])
    print('Constant temperature =', daily_temp)
    print('Constant water potential =', daily_psi)
        
    return Data_Dictionary_Common
