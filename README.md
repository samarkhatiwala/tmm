**ANNOUNCEMENT**: I'm pleased to announce the release of version 3 of the Transport Matrix Method (TMM) software. Included in this release is **tmm4py**, a full-featured wrapper exposing all of the TMM's functionality in Python. tmm4py enables one to write biogeochemical models in Python, Fortran and/or C and run them entirely from Python. To make this possible, the TMM software has been redesigned from the ground up along object oriented principles as a library callable from other languages.

You can read a description of tmm4py and the new TMM library in this *J. Adv. Model. Earth Sys.* article: http://dx.doi.org/10.1029/2025MS005028.

------------------------------------------------------------------------------------------------

This is the Transport Matrix Method (TMM) code repository. It includes the core TMM software library, a Python wrapper tmm4py, and various biogeochemical models adapted to the TMM framework. The TMM software is written using the PETSc framework (http://www.mcs.anl.gov/petsc/). 

**How to cite**: If you use the TMM, please cite Khatiwala et al. (2005; https://doi.org/10.1016/j.ocemod.2004.04.002) and Khatiwala (2007; https://doi.org/10.1029/2007GB002923). If you use tmm4py, please cite Khatiwala (2025; http://dx.doi.org/10.1029/2025MS005028). If you use any of my transport matrices or related data please cite this GitHub page as the source, as well as Khatiwala (2007) and the papers for the relevant transport matrix configurations listed here: https://sites.google.com/view/samarkhatiwala-research-tmm. Thank you!

**IMPORTANT**: Please do NOT post the TMM code or any data files downloaded from this or related websites on your own github or other website. See LICENSE.txt for licensing information.

Feel free to email if you have any questions: <samkat6@gmail.com>

------------------------------------------------------------------------------------------------

To install TMM/tmm4py and any required libraries you will need a C compiler. For tmm4py you will also need Python 3.9 and up, and the numpy and Cython modules. A Fortran compiler is required to run the Fortran examples. See installation instructions below.

To run the examples in the instructions below, download this archive which contains all the necessary input data: https://drive.google.com/file/d/1ou1kamcI-ZE6BWHpnAs2LaSmO1uZeXyd/view?usp=sharing. Unpack and place the 'examples' directory in $TOPDIR/ (see below). Otherwise, you will need to download transport matrices and other data and scripts from https://sites.google.com/view/samarkhatiwala-research-tmm and generate the input data yourself (instructions below).

In the following, it is assumed you're using bash as your shell. Change the commands appropriately if you're using a different shell (e.g., 'export' to 'setenv', etc for tcsh).

###### Notes:

(1) The dependencies for TMM are: PETSc. The dependencies for tmm4py are: PETSc, petsc4py and TMM. PETSc in turn depends on MPI and BLAS/LAPACK.

(2) On many HPC systems, optimized MPI and BLAS/LAPACK libraries will be already installed. PETSc may also be already installed and you just have to set the PETSC_DIR and PETSC_ARCH environment variables (see https://petsc.org/main/install/).

(3) You can install MPI either from source, using a package manager (conda, brew, etc) or with pip. Basic instructions for building from source are provided below. You can also install MPI via pip, assuming binaries are available for your platform (they appear to be for the most common ones). However, it seems that this doesn't install compiler wrappers for Fortran. So if you want to use models written in Fortran then you should install MPI from source or with a package manager.

(4) PETSc and petsc4py can be installed with pip but instructions for building from source are also given below (preferable if you want to optimize it for your hardware and software libraries).

##### Quick Start guide:

(1) Install MPI:

pip install mpich
OR
pip install openmpi

###### Notes: 

(i) If you already have MPI installed then skip this step. Just make sure before you go to the next step that the wrapper compilers (mpicc etc) are being correctly picked up (my experience is that when using pip to install MPI, the correct directory is automatically added to your search path).

(ii) If you want to use Fortran then follow the instructions below for installing MPI from source.

(2) Install PETSc and petsc4py:

pip install petsc4py

(3) Determine path of installed PETSc and set environment variables for subsequent steps:

export PETSC_DIR=\`python3 -c "from petsc import get_petsc_dir; print(get_petsc_dir())"\`  
export PETSC_ARCH=''

(4) Build the TMM library. In the following, the path of the directory in which this README.txt file is located is in the environment variable 'TOPDIR'.

cd $TOPDIR  
make  
make install

By default this will install the TMM library and include files in a subdirectory 'TMM', i.e., $TOPDIR/TMM. If you want to install it in a different location, set the PREFIX variable, e.g.:

make install PREFIX=/home/userid/TMM

When installation is complete, the path to the installation location will be printed (the same as PREFIX if you provide it above). Set the environment variable TMM_DIR to point to this directory, e.g., for the default location:

export TMM_DIR=$TOPDIR/TMM

Note: You can delete the object files created above and in any of the subsequent examples with 'make cleanall'

(5) Lastly, build the tmm4py extension module (replace python3 with python if necessary):

cd $TOPDIR  
make tmm4py PREFIX=$TMM_DIR

This will produce tmm4py.cpython-XXX.so, where XXX is architecture dependent, and install it in the same location as the TMM library. Also installed is tmm.py, a generic driver script that can be used to run most models, whether implemented in Python or Fortran. 

Modify PYTHONPATH so that these can be found:

export PYTHONPATH=$PYTHONPATH:$TMM_DIR/lib

And that's it. For a model in pure Python, libtmm.so and tmm4py.cpython-XXX.so are all you need.

(6) To check if everything is working, run the Age example:

cd $TOPDIR/examples/idealage  
cp -p $TOPDIR/models/idealage/code/python/Age.py \.  
cp -p $TOPDIR/models/idealage/runscript/* \.  

mpiexec -np 4 python3 run_model_age_py.py

(If using a different MPI launcher modify the above command accordingly. For example, if you installed MPI via PETSc as in the instructions below, this would be $PETSC_DIR/$PETSC_ARCH/bin/mpiexec.)

(7) You can similarly run the other (pure) Python examples in the examples directory. These require the numba Python module.

MOPS:

cd $TOPDIR/examples/ironmops  
cp -p $TOPDIR/models/ironmops/code/python/* .  
cp -p $TOPDIR/models/ironmops/runscript/* .  

and follow the instructions in 'runscript'

Pa/Th:

cd $TOPDIR/examples/protactinium_thorium  
cp -p $TOPDIR/models/protactinium_thorium/code/python/* .  
cp -p $TOPDIR/models/protactinium_thorium/input/MOBI_tracer_names.txt .  
cp -p $TOPDIR/models/protactinium_thorium/input/get_tracer_names.sh .  
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .  

and follow the instructions in 'runscript'

##### Running the Fortran examples and/or using the 'classic' TMM C driver:

The following assumes you have installed PETSc with Fortran support. Whether that is the case will depend on how you installed MPI. If you installed MPI using pip as per the instructions above then it is likely MPI was installed *without* the Fortran wrapper compilers (mpif90 etc). In that case you should install MPI with another package manager or build it from source as per the instructions below. Then install PETSc before proceeding with these instructions.

(1) Build TMM and tmm4py as per instructions in the Quick Start guide above.

(2) If you're planning to use the 'classic' TMM C driver, do:

cd $TOPDIR  
make -f Makefile_driver tmm

This compiles tmm_main.c along with stubs for tmm_external*, tmm_monitor.c and tmm_misfit.c into a standalone executable 'tmm'. You can use this for the Age example (next step) as it reads in the forcing term from file.

(3) Age: To run the 'classic' TMM C driver, copy the executable $TOPDIR/tmm you built in Step 2 above:

cd $TOPDIR/examples/idealage  
cp -p $TOPDIR/tmm .  
cp -p $TOPDIR/models/idealage/runscript/* .

and follow the instructions in 'runscript'.

(4) Fortran models: Models written in Fortran can be run with either the 'classic' TMM C driver or via tmm4py. In both cases, you first compile the model as a dynamic library and then link it to either the TMM C driver or a Python gateway module (which makes the model callable from Python).

(i) MOPS:

First build MOPS into a dynamic library:

cd $TOPDIR/models/ironmops/code  
make rmops

The above creates libmops.so. (The Makefile targets several different versions of MOPS with different compiler directives; 'rmops' is the version with runoff, carbon and alkalinity switched on.) 

To run the model from Python, compile the Python gateway interface into a Python extension module:

export ARCHFLAGS=''  
python3 setup.py build_ext --inplace

This compiles model_tmm_interface.pyx (which contains the gateway functions that call the actual external forcing and other C wrapper functions) into mops.cpython-XXX.so. Set PYTHONPATH so that this module can be found:

export PYTHONPATH=$PYTHONPATH:$TOPDIR/models/ironmops/code

To run the model with the TMM C driver, compile the C driver tmm_main.c and link to libmops.so to create a standalone executable tmmmops:

make -f Makefile_standalone tmmmops  
make -f Makefile_standalone cleanall

To run via the Python driver:

cd $TOPDIR/examples/ironmops  
cp -p $TOPDIR/models/ironmops/runscript/* .

and follow the instructions in 'runscript'.

To run via the C driver:

cd $TOPDIR/examples/ironmops  
cp -p $TOPDIR/models/ironmops/runscript/* .  
cp -p $TOPDIR/models/ironmops/code/tmmmops .

and follow the instructions in 'runscript'.

(ii) PaTh:

First download MOBI, which is embedded within the UVic ESCM model. You can get it from https://github.com/OSU-CEOAS-Schmittner/UVic2.9. Download the 'source' and 'updates' directories into a directory (e.g., UVic_ESCM). Set the environment variable UVICESCMROOT to point to this directory, e.g.:

export UVICESCMROOT=$HOME/UVic_ESCM

Next, build MOBI into a dynamic library:

cd $TOPDIR/models/mobi/code  
cp -p MOBI_TMM_OPTIONS.h_PaTh MOBI_TMM_OPTIONS.h  
make

The above creates libmobi.so.

To run the model from Python, compile the Python gateway interface into a Python extension module:

export ARCHFLAGS=''   
python3 setup.py build_ext --inplace

This compiles model_tmm_interface.pyx (which contains the gateway functions that call the actual external forcing and other C wrapper functions) into mobi.cpython-XXX.so. Set PYTHONPATH so that this module can be found:

export PYTHONPATH=$PYTHONPATH:$TOPDIR/models/mobi/code

To run the model with the TMM C driver, compile the C driver tmm_main.c and link to libmobi.so to create a standalone executable tmmmobi:

make -f Makefile_standalone tmmmobi  
make -f Makefile_standalone cleanall

To run with the Python driver:

cd $TOPDIR/examples/protactinium_thorium  
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .

and follow the instructions in 'runscript'.

To run with the C driver:

cd $TOPDIR/examples/protactinium_thorium  
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .  
cp -p $TOPDIR/models/mobi/code/tmmmobi .

and follow the instructions in 'runscript'.

##### Instructions for building MPI from source

(1) Download and unpack the source: https://www.mpich.org/

(2) Configure (with gcc and gfortran in this example):

cd mpich-4.3.0 # modify for the version you downloaded  
FC=gfortran  
./configure CC=gcc --prefix=/usr/local/mpich

(3) Build and install:

make  
sudo make install  
export PATH=$PATH:/usr/local/mpich/bin

##### Instructions for building PETSc from source (also installs MPI and other libraries)

###### Notes:

(i) See https://petsc.org/main/install/ for detailed installation instructions

(ii) On HPC systems, optimized BLAS/LAPACK (Intel MKL, Cray LibSci, etc) and MPI libraries are usually already available. If so, you should use these when configuring PETSc. In particular, use the compiler 'wrappers' mpicc and mpif90 instead of the actual compilers, and specify the location of the BLAS/LAPACK libraries with: --with-blaslapack-dir=...

(iii) If you plan to use tmm4py for 'production' calculations then configure it with the optimization options appropriate for your compiler (the default is to build PETSc in 'debugging' mode). For the GNU compiler, it would be something like: --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native'.

(1) Download a recent tarball from https://petsc.org/main/install/download/ and unpack.

(2) Set the environment variable PETSC_DIR to the path of the directory you just unpacked, e.g.:

export PETSC_DIR=$HOME/petsc-3.23.3

(3) Set the environment variable PETSC_ARCH to some convenient name for your installation (e.g., arch-darwin or arch-linux). You can have different installations under the same PETSC_DIR (e.g., one for debugging and one for production work) and switch between them by setting PETSC_ARCH accordingly, e.g.:

export PETSC_ARCH=arch-arm-darwin-opt

(4) Configure and build PETSc/petsc4py (the command below assumes your C and Fortran compilers are gcc  and gfortran, respectively, but you can change them, e.g., clang, icc, ifort etc. As mentioned above, replace these with mpicc and mpif90 if you already have MPI installed (in which case also remove the --download-mpich option). If you don't have or want to use a Fortran compiler, then set --with-fc=0.

cd $PETSC_DIR/  
./configure --with-cc=gcc --with-cxx=0 --with-fc=gfortran --download-f2cblaslapack --download-mpich --with-shared-libraries=1 --with-petsc4py=1

Then build PETSc by executing the make command printed out by configure above. It will be something like:

make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all

(5) When compilation is finished, PETSc will print out the path to where it installed the petsc4py library (typically it is $PETSC_DIR/$PETSC_ARCH/lib). Add this location to the PYTHONPATH environment variable, e.g.:

export PYTHONPATH=$PYTHONPATH:$PETSC_DIR/$PETSC_ARCH/lib

You're done!

Note: The above configure command will download and build BLAS/LAPACK and the MPICH MPI library for you. You just have to use the PETSc-installed MPI when launching your program with mpiexec or mpirun. PETSc will install them in: $PETSC_DIR/$PETSC_ARCH/bin. You can add this directory to your path or specify the full path, for example:

$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -np 4 ...

##### Instructions for generating input data

To perform simulations you need transport matrices (TMs) and other forcing data. The relevant data and scripts to generate the input files can be downloaded from https://sites.google.com/view/samarkhatiwala-research-tmm.

(1) Download the Matlab scripts and add to your Matlab path. (Python scripts to perform similar operations are currently under development. Email me if you want them.)

(2) Download transport matrices and related data for the ocean model of your choice.  Currently, there are three configurations of MITgcm available online (and several others based on the UVic Earth System Model that I am happy to make available if you email me). For each, download the TMs and other associated data (e.g., MITgcm_ECCO.tar). Unpack. Make a note of the path to this directory (e.g., /mydisk/somewhere/MITgcm_ECCO). We will need it later. For some experiments you may also find it useful to download various miscellaneous data (and adjust paths accordingly in the provided Matlab scripts). 

(3) For each model, e.g., $TOPDIR/models/mops2.0, the input/ directory contains Matlab script(s) to generate input data (e.g., make_input_files_for_mops_model.m) and read model output (e.g., load_output.m).

Copy the matlab and run scripts to the directory in which you want to generate the files.

cp -p -R $TOPDIR/models/mops2.0/input/* .

Change the variable base_path at the very top of the make_input...m script to point to the top level of any of the transport matrix configurations you downloaded in #2 above. Execute the make_input* script (e.g., make_input_files_for_mops_model.m). It should generate all necessary input data. With luck! (If there is missing data email me for it.)