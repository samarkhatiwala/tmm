This is the Transport Matrix Method (TMM) code repository. It includes both 
the core TMM time-stepping driver code (under driver/), as well as various 
biogeochemical models (under models/) adapted to the TMM framework. The driver 
code and interface to models are written using the PETSc framework 
(http://www.mcs.anl.gov/petsc/) but you don't need the code to use the TMM. 
Simply skip to step (3) below. Otherwise keep reading and if you have any 
questions feel free to email: Samar Khatiwala <samark@earth.ox.ac.uk>

For a quick overview of the TMM and the PETSc driver also have a look at this 
excellent presentation by Iris Kriest:
https://ftp.geomar.de/users/ikriest/TMM/TMM-2015-April.pdf

Quick-start instructions:

1) Install PETSc (http://www.mcs.anl.gov/petsc/). The TMM driver code is compatible 
with PETSc version 3.4.x and up. Note that if you're using version 3.6.x you will 
need to modify the makefiles in each example directory by replacing the include 
lines at the top with:
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

2) Download Matlab scripts and add to your Matlab path:
http://kelvin.earth.ox.ac.uk/spk/Research/TMM/tmm_matlab_code.tar.gz

3) Download transport matrices and related data for the ocean model of your 
choice: http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
Currently, there are 3 configurations of MITgcm available. For each, download 
the TMs (e.g., MITgcm_ECCO.tar) and other associated data (MITgcm_ECCO_model_data.tar). 
Unpack both. Then move the contents of the second directory into that of the first. 
Make a note of the path to this directory (e.g., /mydisk/somewhere/MITgcm_ECCO).
We will need it later.

4) Make a local directory and checkout the TMM driver and model codes:
cd $HOME
mkdir TMM
cd TMM/
git clonehttps://github.com/samarkhatiwala/tmm.git

5) For each model, e.g., tmm/models/petsc3.4/mitgchem/ there is a Matlab script to generate 
all input data (e.g., tmm/models/petsc3.4/mitgchem/matlab/make_input_files_for_migchem_dic_biotic_model.m). 
Change the top-level path at the very top of the script (and, depending on the model, some paths to other 
data files), change any other options you want, and execute. It should generate all necessary 
input data. With luck! (If there is missing data email me for it.)

6) Compile code (e.g., for mitgchem):
cd $HOME/TMM/tmm/
mkdir Work/
cd Work
cp -p $HOME/TMM/tmm/driver/petsc3.4/* .
cp -p $HOME/TMM/tmm/models/petsc3.4/mitgchem/src/* .
make tmmmitgchem

7) Copy all input data created in step 5 above and run scripts 
(e.g., $HOME/TMM/tmm/models/petsc3.4/mitgchem/runscripts) to Work/

8) Execute model using the example run scripts.
