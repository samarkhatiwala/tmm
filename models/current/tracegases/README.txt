This model simulates trace gases in the ocean. Currently, N2O, CFC-11, CFC-12, and SF6 
are implemented following the OMIP protocols described in Orr et al. (2016): 
Geosci. Model Dev. Discuss., doi:10.5194/gmd-2016-155, 2016.

To specify the gas being simulated change the following command line options in the run script: 
-gas_id XX, where XX = 1 (N2O), 2 (CFC11), 3 (CFC12) and 4 (SF6)
-mixing_ratio_scale_factor XX, where XX = 1.e-9 (N2O), 1.e-12 (CFC11, CFC12 or SF6)
