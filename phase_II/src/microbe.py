# code B > src > microbe
# Elsa Abs, Feb 2020
"""
This microbe.py module has one class and two functions.
    
    Microbe():                class
    microbe_osmo_psi():       function; moisture modifier of inducible osmolyte production efficiency
    microbe_mortality_prob(): function; cell mortality probability
"""

import sys
import numpy as np
import pandas as pd
from utility import LHS

class Microbe():
    """
    This class holds all variables related to microbes.

    Methods involving composition,stoichiometry, enzyme and gene production, as well as responses to environmental factors.
    These methods include:
        1) microbial_community_initialization(): initialize microbial community on the spatial grid
        2) minimum_cell_quota():          get minimum ratios
        3) microbe_enzyme_gene():         derive taxon-specific genes for enzyme
        4) microbe_osmolyte_gene():       derive genes encoding for osmolytes
        5) microbe_uptake_gene():         derive transporter genes
        6) microbe_uptake_cost():         metabolic cost of producing transporter
        7) microbe_enzproduction_rate():  cost of producing enzymes
        8) microbe_osmoproduction_rate(): cost of producing osmolytes
        9) microbe_drought_tol():         microbial drougth tolerance
       10) microbe_mortality():           paramters pertaining to microbial mortality
    """
    
    def __init__(self,runtime,data_preinit,residents_grid,invaders_pre_grid):
        """
        The constructor of Microbe class.

        Parameters:
            runtime:     user-specified runtime parameters
            paramteters: input dataframe of all parameters
        """
        
        # Import Data Dictionary residents and invaders from preinitialization.py
        data_preinit_residents = data_preinit[0]
        data_preinit_invaders = data_preinit[1]
        
        # Parameters from data_preinit for EnzGenes
        self.n_taxa_res = data_preinit_residents['n_taxa']
        self.n_taxa_inv = data_preinit_invaders['n_taxa']
        self.n_taxa = self.n_taxa_res + self.n_taxa_inv
        self.n_enzymez_res = data_preinit_residents['n_enzymes']
        self.n_enzymez_inv = data_preinit_invaders['n_enzymes']
        self.EnzGenes_residents = data_preinit_residents['EnzGenes'][0:self.n_taxa_res]
        self.EnzGenes_invaders = data_preinit_invaders['EnzGenes'][0:self.n_taxa_inv]
        self.TaxMon_residents = data_preinit_residents['TaxMon_final'][0:self.n_taxa_res]
        self.TaxMon_invaders = data_preinit_invaders['TaxMon_final'][0:self.n_taxa_inv]
        
        # Parameters from data_preinit for Uptake
        self.UptakeGenes_trait_residents = data_preinit_residents['UptakeGenes_trait'][0:self.n_taxa_res]
        self.UptakeGenes_trait_invaders = data_preinit_invaders['UptakeGenes_trait'][0:self.n_taxa_inv]
        self.UptakeGenesCost_residents = data_preinit_residents['UptakeGenesCost'][0:self.n_taxa_res]
        self.UptakeGenesCost_invaders = data_preinit_invaders['UptakeGenesCost'][0:self.n_taxa_inv]
        self.Microbes_residents  = residents_grid.astype('float64')
        self.Microbes_pre_invaders  = invaders_pre_grid.astype('float64')
        self.gridsize = data_preinit_residents['gridsize']
        self.n_uptake_res = data_preinit_residents['n_uptake']
        self.n_uptake_inv = data_preinit_invaders['n_uptake']
        self.n_uptake = self.n_uptake_res + self.n_uptake_inv
        self.TaxUpt_residents = data_preinit_residents['TaxUpt'][0:self.n_taxa_res]
        self.TaxUpt_invaders = data_preinit_invaders['TaxUpt'][0:self.n_taxa_inv]
        
        # Parameters from data_preinit for Metabolism
        self.n_osmolytes_res = data_preinit_residents['n_osmolytes']
        self.n_osmolytes_inv = data_preinit_invaders['n_osmolytes']
        self.n_osmolytes = self.n_osmolytes_res + self.n_osmolytes_inv
        self.OsmoProdConsti_residents = data_preinit_residents['OsmoProdConsti'][0:self.n_taxa_res]
        self.OsmoProdConsti_invaders = data_preinit_invaders['OsmoProdConsti'][0:self.n_taxa_inv]
        self.n_enzymes_res = data_preinit_residents['n_enzymes']
        self.n_enzymes_inv = data_preinit_invaders['n_enzymes']
        self.n_enzymes = self.n_enzymes_res + self.n_enzymes_inv
        self.EnzProdConstit_residents = data_preinit_residents['EnzProdConstit'][0:self.n_taxa_res]
        self.EnzProdConstit_invaders = data_preinit_invaders['EnzProdConstit'][0:self.n_taxa_inv]
        self.OsmoProdInduci_residents = data_preinit_residents['OsmoProdInduci'][0:self.n_taxa_res]
        self.OsmoProdInduci_invaders = data_preinit_invaders['OsmoProdInduci'][0:self.n_taxa_inv]
        self.EnzProdInduce_residents = data_preinit_residents['EnzProdInduce'][0:self.n_taxa_res]
        self.EnzProdInduce_invaders = data_preinit_invaders['EnzProdInduce'][0:self.n_taxa_inv]
                
        # Parameters from data_preinit for Mortality
        self.MinRatios_residents = data_preinit_residents['MinRatios'][0:self.n_taxa_res]
        self.MinRatios_invaders = data_preinit_invaders['MinRatios'][0:self.n_taxa_inv]
        self.basal_death_prob_residents = data_preinit_residents['basal_death_prob'][0:self.n_taxa_res]
        self.basal_death_prob_invaders = data_preinit_invaders['basal_death_prob'][0:self.n_taxa_inv]
        self.death_rate_residents = data_preinit_residents['death_rate'][0:self.n_taxa_res]
        self.death_rate_invaders = data_preinit_invaders['death_rate'][0:self.n_taxa_inv]
        self.TaxDroughtTol_residents = data_preinit_residents['TaxDroughtTol'][0:self.n_taxa_res]
        self.TaxDroughtTol_invaders = data_preinit_invaders['TaxDroughtTol'][0:self.n_taxa_inv]
        
        # Parameters from data_preinit for Reproduction
        self.fb_residents = data_preinit_residents['fb'][0:self.n_taxa_res]
        self.fb_invaders = data_preinit_invaders['fb'][0:self.n_taxa_inv]
        
    
    
    def microbial_community_initialization(self):
        """
        Merge residents and invaders communities.
                
        Parameters:
            n_taxa_res:          number of resident taxa; from Dictionary_Residents
            n_taxa_inv:          number of invaders taxa; from Dictionary_Invaders
            gridsize:            gridsize taken from runtime (I use the same for codes A, B, C)
            Microbes_residents:  dataframe(taxa*3); 3 columns are C, N, P; from Dictionary_Residents
            Microbes_invaders:   dataframe(taxa*3); 3 columns are C, N, P; from Dictionary_Invaders
            fb_residents:        1D array(taxa);index of bacteria (0) vs fungi (1) over the grid; from Dictionary_Residents
            fb_invaders:         1D array(taxa);index of bacteria (0) vs fungi (1) over the grid; from Dictionary_Invaders
        Return:
            n_taxa:              number of resident + invader taxa (regardless of whether biomass=0); for Dictionary_Common
            Microbes_df:         dataframe(taxa*3); 3 columns are C, N, P; for Dictionary_Common
            fb:                  1D array(taxa);index of bacteria (0) vs fungi (1) over the grid; for Dictionary_Common
        """
        
            
        residents      = sys.argv[5]
        invaders       = sys.argv[6]
        enzmax         = sys.argv[7]
        seed_index     = sys.argv[8]
        
        ### n_taxa
        n_taxa = self.n_taxa
        
        ### Create invaders_grid.csv from residents_grid.csv
        self.Microbes_invaders = self.Microbes_pre_invaders.copy(deep=True)
        if invaders != '0':
            pi = 0.4
            choose_invader = np.random.choice([1,0], int(self.n_taxa_inv)*int(self.gridsize),replace=True, p=[pi,(1-pi)]).astype('int32')
            self.Microbes_invaders.loc[choose_invader==0] = 0
            export_csv = self.Microbes_invaders.to_csv (r'./grids/invaders_grid_'+invaders+enzmax+seed_index+'.csv', header=True, index=True)
        
        ### Microbes_df
        # Change index of Tax in Microbes invaders (more complicated than other tables because here, data are different in different boxes)
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        microbes_inv_index = []
        for i in np.arange(self.gridsize):
            microbes_inv_index += taxa_inv_index
        self.Microbes_invaders.index = microbes_inv_index
        
        # Goal: integrate invaders biomass per box
        # (1) Add index level to both grids to be then able to sort by index level
        ## (1.a) Residents
        list_Res = []
        Res = np.zeros((self.n_taxa_res,self.n_taxa_res))
        for i in np.arange(self.gridsize):
            Res = self.Microbes_residents.iloc[i*self.n_taxa_res:(i+1)*self.n_taxa_res]
            list_Res.append(Res)
        out_Res = pd.concat(list_Res, keys=list(range(0,self.gridsize)))
        ## (1.b) Invaders
        list_Inv = []
        Inv = np.zeros((self.n_taxa_inv,self.n_taxa_inv))
        for i in np.arange(self.gridsize):
            Inv = self.Microbes_invaders.iloc[i*self.n_taxa_inv:(i+1)*self.n_taxa_inv]
            list_Inv.append(Inv)
        out_Inv = pd.concat(list_Inv, keys=list(range(0,self.gridsize)),sort=False)
        # (2) Concatenate Microbes residents and invaders by index level
        InvRes = pd.concat([out_Res, out_Inv]).sort_index(level=0,sort_remaining=False).reset_index(level=0, drop=True)
        Microbes_df = InvRes.copy(deep=True)

        ### fb
        # Concatenate fb residents and invaders (no need to change index because fb is an array)
        fb = np.concatenate((self.fb_residents,self.fb_invaders),axis=None).astype(int)
        
        return n_taxa, Microbes_df, fb
        
        
    def microbe_enzyme_gene(self):
        """
        Derive taxon-specific enzyme genes, so I can know which taxa has which enzyme

        Return:
            EnzGenes: Rows are taxa; cols are genes;values: 0/1
        """
      
        # Change index of Tax in EnzGenes invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.EnzGenes_invaders.index = taxa_inv_index
        
        # Change index of Enz in EnzGenes invaders
        enzymes_inv_index = []
        for i in np.arange(self.n_enzymez_inv):
            enzymes_inv_index.append("Enz" + str(self.n_enzymez_res+i+1))
        EnzGenes_invaders_T = self.EnzGenes_invaders.T
        EnzGenes_invaders_T.index = enzymes_inv_index
        self.EnzGenes_invaders = EnzGenes_invaders_T.T
        
        # Concatenate TaxUpt residents and invaders
        EnzGenes_NA = pd.concat([self.EnzGenes_residents,self.EnzGenes_invaders], axis=0, sort=False)
        EnzGenes = EnzGenes_NA.fillna(0)
        
        return EnzGenes
        
        
    def microbe_uptake_monomer(self):
        """
        Derive the which transporters each taxon has, and which monomers each taxon has access to through enzymes + randomly added + transporters.
        
        Parameters:
            n_taxa_res:          number of resident taxa; from Dictionary_Residents
            n_taxa_inv:          number of invaders taxa; from Dictionary_Invaders
            TaxUpt_residents:    dataframe(taxa*transporters); from Dictionary_Residents
            TaxUpt_invaders:     dataframe(taxa*transporters); from Dictionary_Invaders
            TaxMon_residents
            TaxMon_invaders
        Return:
            TaxUpt:              dataframe(taxa*transporters); taxon- and transporter-specific association; for Dictionary_Common
            TaxMon:
        """
    
        # Change index of Tax in TaxUpt invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.TaxUpt_invaders.index = taxa_inv_index

        # Change index of Upt in TaxUpt invaders
        uptake_inv_index = []
        for i in np.arange(self.n_uptake_inv):
            uptake_inv_index.append("Upt" + str(self.n_uptake_res+i+1))
        TaxUpt_invaders_T = self.TaxUpt_invaders.T
        TaxUpt_invaders_T.index = uptake_inv_index
        self.TaxUpt_invaders = TaxUpt_invaders_T.T
        
        # Concatenate TaxUpt residents and invaders
        TaxUpt_NA = pd.concat([self.TaxUpt_residents,self.TaxUpt_invaders], axis=0, sort=False)
        TaxUpt = TaxUpt_NA.fillna(0)
        
        # Change index of Tax in TaxMon invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.TaxMon_invaders.index = taxa_inv_index
        
        # Concatenate EnzGenes residents and invaders
        TaxMon = pd.concat([self.TaxMon_residents,self.TaxMon_invaders], axis=0, sort=False)
        
        return TaxUpt, TaxMon
        
    
    def microbe_uptake_cost(self):
        """
        Derive the taxon-specific cost of every single gene of uptake transporter.
        
        Note this cost (in terms of Fraction of Biomass C) is same among different genes from the same taxon
        
        Parameters:
            n_taxa_res:                    number of resident taxa; from Dictionary_Residents
            n_taxa_inv:                    number of invaders taxa; from Dictionary_Invaders
            UptakeGenes_trait_residents:   series(taxa); from Dictionary_Residents
            UptakeGenes_trait_invaders:    series(taxa); from Dictionary_Invaders
            UptakeGenesCost_residents:     dataframe(taxa*monomers); from Dictionary_Residents
            UptakeGenesCost_invaders:      dataframe(taxa*monomers); from Dictionary_Invaders
        Returns:
            UptakeProd_series:             series(taxa); taxon-specific transporter cost; for Dictionary_Common
            UptakeGenes_Cost:              dataframe(taxa*monomers); taxon- and gene-specific transporter cost; for Dictionary_Common
        """                         
        
        # Change index of Tax in UptakeGenes_trait invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.UptakeGenes_trait_invaders.index = taxa_inv_index
        
        # Concatenate UptakeGenes_trait residents and invaders --> to make UptakeProd_series
        UptakeProd_series = pd.concat([self.UptakeGenes_trait_residents,self.UptakeGenes_trait_invaders], axis=0, sort=False)
        
        # Change index of Tax in UptakeGenesCost_invaders_T invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.UptakeGenesCost_invaders.index = taxa_inv_index
        
        # Concatenate UptakeGenesCost residents and invaders --> to make UptakeGenes_Cost
        UptakeGenes_Cost = pd.concat([self.UptakeGenesCost_residents,self.UptakeGenesCost_invaders], axis=0, sort=False)
        
        return UptakeProd_series, UptakeGenes_Cost
        
    
    def microbe_enzproduction_rate(self):
        """
        Derive the taxon-specific fraction of 'available C' as enzymes: note that.

        This fraction only varies with taxon, which is independent of gene within a taxon
        
        Parameters:
            n_taxa_res:                  number of residents taxa; from Dictionary_Residents
            n_taxa_inv:                  number of invaders taxa; from Dictionary_Invaders
            n_enzymes_res:               number of residents enzymes; from Dictionary_Residents
            n_enzymes_inv:               number of invaders enzymes; from Dictionary_Invaders
            EnzProdConstit_residents:    dataframe(taxa*enzymes); from Dictionary_Residents
            EnzProdConstit_invaders:     dataframe(taxa*enzymes); from Dictionary_Invaders
            EnzProdInduce_residents:     dataframe(taxa*enzymes); from Dictionary_Residents
            EnzProdInduce_invaders:      dataframe(taxa*enzymes); from Dictionary_Invaders
        Returns:
            Tax_Consti_Enzyme_C:         dataframe(taxa*enzymes); for Dictionary_Common
            Tax_Induce_Enzyme_C:         dataframe(taxa*enzymes); for Dictionary_Common
        """
        
        # Change index of Tax in EnzProdConstit invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.EnzProdConstit_invaders.index = taxa_inv_index
        
        # Change index of Enz in EnzProdConstit invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        EnzProdConstit_invaders_T = self.EnzProdConstit_invaders.T
        EnzProdConstit_invaders_T.index = enzyme_inv_index
        self.EnzProdConstit_invaders = EnzProdConstit_invaders_T.T
        
        # Concatenate EnzProdConstit residents and invaders
        Tax_Consti_Enzyme_C_NA = pd.concat([self.EnzProdConstit_residents,self.EnzProdConstit_invaders], axis=0, sort=False)
        Tax_Consti_Enzyme_C = Tax_Consti_Enzyme_C_NA.fillna(0)
        
        # Change index of Tax in EnzProdInduce invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.EnzProdInduce_invaders.index = taxa_inv_index
        
        # Change index of Enz in EnzProdInduce invaders
        enzyme_inv_index = []
        for i in np.arange(self.n_enzymes_inv):
            enzyme_inv_index.append("Enz" + str(self.n_enzymes_res+i+1))
        EnzProdInduce_invaders_T = self.EnzProdInduce_invaders.T
        EnzProdInduce_invaders_T.index = enzyme_inv_index
        self.EnzProdInduce_invaders = EnzProdInduce_invaders_T.T
        
        # Concatenate EnzProdInduce residents and invaders
        Tax_Induce_Enzyme_C_NA = pd.concat([self.EnzProdInduce_residents,self.EnzProdInduce_invaders], axis=0, sort=False)
        Tax_Induce_Enzyme_C = Tax_Induce_Enzyme_C_NA.fillna(0)
    
        return Tax_Consti_Enzyme_C, Tax_Induce_Enzyme_C
    
    
    def microbe_osmoproduction_rate(self):
        """
        Distribution of osmolyte production rate (i.e.,proportion of available C as osmolytes).

        Parameters:
            n_taxa_res:                  number of residents taxa; from Dictionary_Residents
            n_taxa_inv:                  number of invaders taxa; from Dictionary_Invaders
            n_osmolytes_res:             number of residents osmolytes; from Dictionary_Residents
            n_osmolytes_inv:             number of invaders osmolytes; from Dictionary_Invaders
            OsmoProdConsti_residents:    dataframe(taxa*osmolytes); from Dictionary_Residents
            OsmoProdConsti_invaders:     dataframe(taxa*osmolytes); from Dictionary_Invaders
        Returns:
            Tax_Consti_Osmo_C:           dataframe(taxa*osmolytes); for Dictionary_Common
            Tax_Induci_Osmo_C:           dataframe(taxa*osmolytes); for Dictionary_Common
        """
        
         # Change index of Tax in OsmoProdConsti invaders
        taxa_inv_index = []
        print('self.n_osmolytes_res',self.n_osmolytes_res)
        print('self.n_osmolytes_inv',self.n_osmolytes_inv)
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.OsmoProdConsti_invaders.index = taxa_inv_index

        # Change index of Osmo in OsmoProdConsti invaders
        osmo_inv_index = []
        for i in np.arange(self.n_osmolytes_inv):
            osmo_inv_index.append("Osmo" + str(self.n_osmolytes_res+i+1))
        OsmoProdConsti_invaders_T = self.OsmoProdConsti_invaders.T
        OsmoProdConsti_invaders_T.index = osmo_inv_index
        self.OsmoProdConsti_invaders = OsmoProdConsti_invaders_T.T
        
        # Concatenate TaxUpt residents and invaders
        Tax_Consti_Osmo_C_NA = pd.concat([self.OsmoProdConsti_residents,self.OsmoProdConsti_invaders], axis=0, sort=False)
        Tax_Consti_Osmo_C = Tax_Consti_Osmo_C_NA.fillna(0)
        
        # Change index of Tax in OsmoProdInduci invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.OsmoProdInduci_invaders.index = taxa_inv_index

        # Change index of Osmo in OsmoProdInduci invaders
        osmo_inv_index = []
        for i in np.arange(self.n_osmolytes_inv):
            osmo_inv_index.append("Osmo" + str(self.n_osmolytes_res+i+1))
        OsmoProdInduci_invaders_T = self.OsmoProdInduci_invaders.T
        OsmoProdInduci_invaders_T.index = osmo_inv_index
        self.OsmoProdInduci_invaders = OsmoProdInduci_invaders_T.T
        
        # Concatenate TaxUpt residents and invaders
        Tax_Induci_Osmo_C_NA = pd.concat([self.OsmoProdInduci_residents,self.OsmoProdInduci_invaders], axis=0, sort=False)
        Tax_Induci_Osmo_C = Tax_Induci_Osmo_C_NA.fillna(0)
        
        return Tax_Consti_Osmo_C, Tax_Induci_Osmo_C
        
        
    def minimum_cell_quota(self):
        """
        This will be used in the mortality calculation.
        
        Parameters:
            n_taxa_res:             number of residents taxa; from Dictionary_Residents
            n_taxa_inv:             number of invaders taxa; from Dictionary_Invaders
            MinRatios_residents:    dataframe(taxa*3); 3 columns are C,N,P; from Dictionary_Residents
            MinRatios_invaders:     dataframe(taxa*3); 3 columns are C,N,P; from Dictionary_Invaders
        Return:
            MinRatios:              dataframe(taxa*3); 3 columns are C,N,P; Optimal stoechiometry - Tolerance; for Dictionary_Common
        """

        # Change index of Tax in MinRatios invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.MinRatios_invaders.index = taxa_inv_index
        
        # Concatenate MinRatios residents and invaders
        MinRatios = pd.concat([self.MinRatios_residents,self.MinRatios_invaders], axis=0, sort=False)
        
        return MinRatios
        
        
    def microbe_mortality(self):
        """
        Derive taxon-specific microbial mortality parameters.
        
        Parameters:
            basal_death_prob_residents: 1D array(taxa); from Dictionary_Residents
            basal_death_prob_invaders:  1D array(taxa); from Dictionary_Invaders
        Returns:
            basal_death_prob:           1D array(taxa); basal death probability; for Dictionary_Common
            death_rate:                 spatial death_rate; 1D array
        """
        
        # Concatenate basal_death_prob residents and invaders
        # No need to change index of invaders this time because basal_death_prob is an array (not a series or a dataframe)
        basal_death_prob = np.concatenate((self.basal_death_prob_residents,self.basal_death_prob_invaders),axis=None)
        
        # Concatenate death_rate residents and invaders
        # No need to change index of invaders because death_rate is an array
        death_rate = np.concatenate((self.death_rate_residents,self.death_rate_invaders),axis=None).astype(int)

        return basal_death_prob, death_rate
        
        
    def microbe_drought_tol(self):
        """
        Derive taxon-specific drought tolerance value.
        
        Drought tolerance is postively correlated with taxon-specific inducible osmotic allocation efficiency.
        
        Parameter:
            n_taxa_res:                  number of residents taxa; from Dictionary_Residents
            n_taxa_inv:                  number of invaders taxa; from Dictionary_Invaders
            TaxDroughtTol_residents:     series(taxa); from Dictionary_Residents
            TaxDroughtTol_invaders:      series(taxa); from Dictionary_Invaders
        Return:
            Tax_tolerance:               series(taxa); for Dictionary_Common
        """

        # Change index of Tax in TaxDroughtTol invaders
        taxa_inv_index = []
        for i in np.arange(self.n_taxa_inv):
            taxa_inv_index.append("Tax" + str(self.n_taxa_res+i+1))
        self.TaxDroughtTol_invaders.index = taxa_inv_index
        
        # Concatenate MinRatios residents and invaders
        Tax_tolerance = pd.concat([self.TaxDroughtTol_residents,self.TaxDroughtTol_invaders], axis=0, sort=False)
        
        return Tax_tolerance
        


def microbe_osmo_psi(wp,alfa,wp_fc,wp_th):
    """
    Derive water potential modifier of inducible osmolyte production.

    Inducible production of osmolytes triggered when PSI declines to a **threshold** value,wp_fc,
    below which the production increases and reaches maxima at water potential of wp_th
    ---------------------------------------------------------------------------
    Parameters:
        wp:    scalar;water potential at a daiy step 
        alfa:  scalar;shape factor quantifying curve concavity; could be distinguished btw bacteria and fungi
        wp_fc: scalar;water potential at field capacity
        wp_th: scalar;water potential threshold   
    Returns:
        f_osmo:scalar; modifier of inducible production of osmoylte   
    References:
        Based on Manzoni et al. 2012 Ecology, a synthesis study. 
    """
    
    if wp >= wp_fc:
        f_osmo = 0.0
    #elif wp <= wp_th:
    #    f_osmo = 1.0
    else:
        x = np.log(wp/wp_fc)
        y = np.log(wp_th/wp_fc)
        
        f_osmo = (x/y)**alfa

    return np.float32(f_osmo)
   
   
   
def microbe_mortality_prob(wp,wp_fc,basal_death_prob,death_rate,Tax_tolerance):
    """
    Microbial mortality probability as a function of water potential and drought tolerance.

    Paramters:
        wp:               scalar; water potential
        wp_fc:            scalar; field capacity
        basal_death_prob: array; basal mortality prob. distinguished btw fungi and bacteria
        death_rate:       array; mortality change rate with moisture
        Tax_tolerance:    dataframe; taxon-specific drought tolerance
    Returns:
        mortality_rate:   taxon-specific mortality probability
    References:
        Allison and Goulden,2017,Soil Biology and Biochemistry
    """
    
    if wp >= wp_fc:
        mortality_rate = basal_death_prob
    else:
        tolerance = Tax_tolerance.to_numpy()
        # option 1
        mortality_rate = basal_death_prob * (1 - (1-tolerance)*(wp-wp_fc)*death_rate)
        # option 2
        #mortality_rate = death_rate * (1/np.exp(tolerance)) * (1 - beta*(wp-wp_fc))
    
    return mortality_rate.astype('float32')
