This is the Transport Matrix Method (TMM) code repository. It includes both 
the core TMM time-stepping driver code (under driver/), as well as various 
biogeochemical models (under models/) adapted to the TMM framework. The driver 
code and interface to models are written using the PETSc framework 
(http://www.mcs.anl.gov/petsc/) but you don't need this code to use the TMM. 
Simply skip to step (3) below. Otherwise keep reading and if you have any 
questions feel free to email: Samar Khatiwala <samar.khatiwala@earth.ox.ac.uk>

IMPORTANT: Please do NOT post the TMM code or any data files downloaded from this or related 
websites on your own github or other website. You do NOT have permission to do so.

If you use the TMM, please cite Khatiwala et al. (2005; https://doi.org/10.1016/j.ocemod.2004.04.002) 
and Khatiwala (2007; https://doi.org/10.1029/2007GB002923). Furthermore, if you use this code 
please also cite Khatiwala (2018; https://doi.org/10.5281/zenodo.1246300). If you use any of 
my transport matrices please cite the relevant paper (see 
http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/). Thank you!

For a quick overview of the TMM and the PETSc driver also have a look at this excellent 
presentation by Iris Kriest: https://ftp.geomar.de/users/ikriest/TMM/MOPS-TMM-2016-June.pdf

Quick-start instructions:

1) Install PETSc (http://www.mcs.anl.gov/petsc/) and set the PETSC_DIR and PETSC_ARCH 
environment variables. The TMM driver code is compatible with PETSc version 3.16.x 
(as of Jan 20, 2022). For the older version of the TMM code compatible with PETSc 3.6.x, 
you can checkout branch petsc3.6.

2) Download Matlab scripts and add to your Matlab path:
http://kelvin.earth.ox.ac.uk/spk/Research/TMM/tmm_matlab_code.tar.gz

3) Download transport matrices and related data for the ocean model of your 
choice: http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
Currently, there are 3 configurations of MITgcm available online (and several 
others based on the UVic Earth System Model that I am happy to make available). 
For each, download the TMs and other associated data (e.g., MITgcm_ECCO.tar). 
Unpack. Make a note of the path to this directory (e.g., /mydisk/somewhere/MITgcm_ECCO). 
We will need it later. For some experiments you may also find it useful to download 
some miscellaneous data here (and adjust paths accordingly in the provided Matlab scripts): 
http://kelvin.earth.ox.ac.uk/spk/Research/TMM/MiscData/

4) Make a local directory and checkout the TMM driver and model codes:
cd $HOME
mkdir TMM
cd TMM/
git clone https://github.com/samarkhatiwala/tmm.git

Set the environment variable TMMROOT to point to the top level of the TMM directory.

5) For each model, e.g., $TMMROOT/models/current/mops2.0/ there is model-specific code (in src/); 
Matlab script(s) in matlab/ to generate input data (e.g., make_input_files_for_mops_model.m) and 
read model output (e.g., load_output.m); and run scripts and other runtime data such as 
namelists in runscripts/. 

To try out one of these models, make a directory, e.g., Test/, and copy the following to it:

cp -p $TMMROOT/models/current/mops2.0/src/Makefile .

(If there is a file *_TMM_OPTIONS.h in src/ copy that as well, e.g.,
cp -p $TMMROOT/models/current/mobi2.0/src/MOBI_TMM_OPTIONS.h .
Edit this file to set compile-time C/Fortran preprocessor options.)

Compile the code:
(For some models you may have to first set additional environment variables as described 
in the corresponding Makefile.)

make mops

Copy the matlab and run scripts:

cp -p -R $TMMROOT/models/current/mops2.0/matlab/* .
cp -p $TMMROOT/models/current/mops2.0/runscripts/* .

Change the variable base_path at the very top of the make_input*.m and load_output*.m scripts 
to point to the top level of any of the transport matrix configurations you downloaded in #3 
above. Execute the make_input* script (e.g., make_input_files_for_mops_model.m). It should 
generate all necessary input data. With luck! (If there is missing data email me for it.)

Execute model using the example run scripts.

Load output using the example load_output.m script.

