IK, 2016-03-09

This  file documents MOPS-Version-2, for use within an optimization environment. 
It has developed from MOPS-Version-1_2, available under github.com/samarkhatiwala/tmm.

See file README.changes_MOPS2.0-MOPS1.2 to see detailed differences
between the code versions.

To use MOPS (or any other BGC model) within an optimization environment,
two tasks have to be accomplished: 

(1) Reading a vector of biogeochemical parameters to be tested (provided by the optimizer)
(2) Computing and writing a misfit function (required by the optimizer)

(1) The procedure of reading parameters is started with command line
    option 
    
    	-bgc_params_file <filename>

    This will read parameter values from a  file, whose name is given by <filename>.
    The parameter file can either be binary (default) or an an ASCII file 
    (compile option -DASCIIPARAMS). 
    
    NOTE: Currently, CMAES (the optimizer) is set up such that it wants to 
    read an ASCII file; so always compile with that option. 
     
    One also has to provide the number of parameters to be read, by 
    command line option 
    
    	-num_bgc_params <number of parameters>. 
    
    NOTE: the number and order of parameters in the parameter file MUST 
    match the number and order in routine mops_biogeochem_set_params.F! 
    Their upper and lower boundaries, as well as the number of parameters 
    to be read and optimized should also match the values given in file "nIter.txt" 

(2) The procedure of calculating the misfit function is started with
    command line option
     
    	-calc_misfit
    
    In this case, the code (tmm_misfit.c) currently reads three observation files, 
    namely "po4_obs.petsc", "no3_obs.petsc", "o2_obs.petsc", and calculates the 
    global RMSE (volume weighted) of annual means. 
    Volume weighting additionally requires a file called "volume_fraction.petsc", 
    containing the fractional volume of each box. 
    The code computes the deviation between annual means; this is 
    invoked by option 
    
    	-average_cost
    
    Currently, there is no other option available (subject to future changes).
    One has to specify over which time interval the model tracers should be 
    integrated to be comparable to and when to start computation. This is done by 
    
	-cost_start_time_step <starttime> -cost_time_steps <averaging interval>
    
    The output of misfit (cost) will be to a file specified by command line option 
    
    	-misfit_file <filename>. 
        
    NOTE: The observation files, as well as "volume_fraction.petsc" MUST 
    match the model geometry (i.e., be of the same format as e.g. input or 
    output vectors of the state variables). 


Most of the tasks associated with reading observations, calculating and
writing the misfit, etc., are done in tmm_misfit.F. The connection to
the actual BGC model (currently MOPS) is via mops_biogeochem_misfit.F.
The biogeochemical model in BGC_MODEL.F writes the tracers to be used
for assessment into the appropriate arrays. Communication among the functions
is by header files mops_biogeochem.h, mops_biogeochem_misfit_data.h and 
tmm_misfit.h.


Further/minor changes to the code: 
----------------------------------

To speed up computation, the code now omits the continuous output of transient runoff.
This can be switched on again via compile option "WRITE_RUNOFF". 

The code further provides the option to output global burial (to be used as
runoff/resupply in the next call to/continuation of the simulation) to a
named file:

	-pickup_runoff_out <filename>

