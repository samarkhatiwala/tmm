**ANNOUNCEMENT**: TMM and tmm4py are now available on PyPI and can be installed in a single step that automatically builds and installs all their dependencies: `pip install tmm4py`. This was a nontrivial undertaking that would not have been possible without [Jamie Carr's](https://github.com/colourfulLanguage) help.

**ANNOUNCEMENT**: I'm pleased to announce the release of version 3 of the Transport Matrix Method (TMM) software. Included in this release is **tmm4py**, a full-featured wrapper exposing all of the TMM's functionality in Python. tmm4py enables one to write biogeochemical models in Python, Fortran and/or C and run them entirely from Python. To make this possible, the TMM software has been redesigned from the ground up along object oriented principles as a library callable from other languages. A huge thank you to [Jamie Carr](https://github.com/colourfulLanguage), without whose help and expertise the tmm4py project would not have gotten off the ground.

You can read a description of tmm4py and the new TMM library in this *J. Adv. Model. Earth Sys.* article: http://dx.doi.org/10.1029/2025MS005028.

------------------------------------------------------------------------------------------------

This is the Transport Matrix Method (TMM) code repository. It includes the core TMM software library, a Python wrapper tmm4py, and various biogeochemical models adapted to the TMM framework. The TMM software is written using the PETSc framework (http://www.mcs.anl.gov/petsc/).

**How to cite**:

This repository and all versions of the TMM and tmm4py software can be cited as: Khatiwala (2018; https://doi.org/10.5281/zenodo.1246299). 

If you use the TMM, please cite Khatiwala et al. (2005; https://doi.org/10.1016/j.ocemod.2004.04.002) and Khatiwala (2007; https://doi.org/10.1029/2007GB002923). If you use tmm4py, please cite Khatiwala (2025; http://dx.doi.org/10.1029/2025MS005028). If you use any of my transport matrices or related data please cite this GitHub page (Khatiwala, 2018;  https://doi.org/10.5281/zenodo.1246299) as the source, as well as Khatiwala (2007) and the papers for the relevant transport matrix configurations listed here: https://sites.google.com/view/samarkhatiwala-research-tmm. Thank you!

**IMPORTANT**: Please do NOT post the TMM code or any data files downloaded from this or related websites on your own github or other website. See LICENSE.txt for licensing information.

Feel free to email if you have any questions: <samkat6@gmail.com>

------------------------------------------------------------------------------------------------

To install TMM/tmm4py and any required libraries you will need at a minimum need a C compiler. For tmm4py you will also need Python 3.9 and up. A Fortran compiler is required to run the Fortran examples.

In the following, it is assumed you're using bash as your shell. Change the commands appropriately if you're using a different shell (e.g.,change  `export` to `setenv`, etc for tcsh).

###### Notes:

(1) The dependencies for TMM are: PETSc. The dependencies for tmm4py are: PETSc, petsc4py and TMM. PETSc in turn depends on MPI and BLAS/LAPACK. petsc4py and tmm4py also depend on NumPy and Cython.

(2) TMM and tmm4py are available as independent packages on PyPI for installation with pip. (TMM is called tmmlib as there was already another package called TMM.) If PETSc and petsc4py are not already available on your system, the installers will build and install them for you.

(3) On many HPC systems, optimized MPI and BLAS/LAPACK libraries will be already installed. PETSc may also be already installed and you just have to set the PETSC_DIR and PETSC_ARCH environment variables (see https://petsc.org/main/install/). While PETSc and petsc4py can be installed with pip, for production runs optimized for your hardware and software libraries it is recommended you build them from source (see instructions below). That said, as described below, even with pip you can control in detail how PETSc is configured and built by setting the PETSC_CONFIGURE_OPTIONS environment variable.

(4) If MPI is not already available, it is recommended that you install it either from source or using a package manager (conda, brew, etc). Basic instructions for this are provided below. You can also install MPI via pip, assuming binaries are available for your platform (they appear to be for the most common ones). However, as of now this doesn't install compiler wrappers for Fortran. So if you want to use models written in Fortran then you should install MPI from source or with a different package manager.

(5) **MacOS users: The PETSc build system runs checks to ensure you have a functioning MPI. On MacOS, because of Apple's annoying security policies, this will throw up a dialog box asking you to allow this. Just click OK. If you don't, the checks will likely fail and, while PETSc will still install, it will switch off its MPI functionality. Bottomline is don't walk away from your computer when performing any of the installation steps below.**

#### Installation using pip

###### Notes:

(i) In the following, if you only want to use the TMM and don't want to use Python/tmm4py, then replace `pip install tmm4py` with `pip install tmmlib` . You can also first install tmmlib with `pip install tmmlib`, set its path with:

```bash
export TMM_DIR=`python3 -c "from tmmlib import get_tmm_dir; print(get_tmm_dir())"`
```

and then follow the instructions below for installing tmm4py.

(ii) It is highly recommended you do a ``pip cache purge`` before you run the install commands below.

(iii) If you run into issues installing with pip add "-v" to the 'pip install xxx' commands below and email me the full output.

(iv) When installing PETSc/petsc4py via pip (either directly or during installation of tmmlib/tmm4py), the PETSc build system will look for any compilers in its path, starting with wrappers named `mpicc`,`mpif90`etc. If your system uses different wrappers (e.g., `mpiifort` instead of `mpif90`) or if you want to use specific compilers, you can specify them in the environment variable `PETSC_CONFIGURE_OPTIONS`before invoking `pip install ...` below. You can also specify optimization flags, paths to optimized BLAS/LAPACK libraries, and any other configuration options that PETSc accepts (see instructions for building PETSc from source below) For example:

```bash
export PETSC_CONFIGURE_OPTIONS='--with-cc=mpicc --with-fc=mpiifort --with-cxx=0 --with-blas-lapack-dir=$MKLROOT --COPTFLAGS="-O2" --FOPTFLAGS="-O2"'
pip install tmm4py
```

Or:

```bash
PETSC_CONFIGURE_OPTIONS='--with-cc=mpicc --with-fc=mpiifort --with-blas-lapack-dir=$MKLROOT --COPTFLAGS="-O2" --FOPTFLAGS="-O2"' pip install tmm4py
```

**(1) If you don't have MPI or PETSc and only want to run models written in Python or use the C driver:**

```bash
pip install mpich
pip install tmm4py
```

**(2) You already have MPI but don't have PETSc:**

Assuming the MPI compiler wrappers (mpicc, mpif90) are in your path, do:

```bash
pip install tmm4py
```

**(3) You already have PETSc installed:**

Set PETSC_DIR/PETSC_ARCH as per PETSc [instructions](https://petsc.org/main/install/). If you used pip to install PETSc then paths can be set like so:

```bash
export PETSC_DIR=`python3 -c "from petsc import get_petsc_dir; print(get_petsc_dir())"`
export PETSC_ARCH=''
```

And then:

```bash
pip install tmm4py
```

**(4) You already have PETSc *and* petsc4py installed:**

The previous instructions will install petsc4py. If you also already have petsc4py installed, either when building PETSc from source or via pip, you need to tell pip about it (see below) so that it won't build it again:

```bash
PETSC_HAS_PETSC4PY=1 pip install tmm4py
```

Notes:

(i) The reason pip needs to be told about an existing petsc4py is because it uses 'build isolation' when compiling packages. Each required dependency package is installed in an isolated environment that knows nothing about what you may already have installed in your system when you invoked 'pip install ...'. In fact, this process is recursive so that the dependencies of each dependency are similarly built in an isolated environment. To get around this and deal with all the possible permutations of how PETSc/petsc4py/TMM/tmm4py can be installed, a custom backend was created that dynamically figures out what the actual dependencies are and what needs to be installed based on your configuration. Thus, if the PETSC_DIR environment variable is detected everything will be built against that version.

(ii) petsc4py and tmm4py depend on NumPy at both build and run time. They also depend on Cython at build time. Ideally we want to use same versions of NumPy and Cython for petsc4py and tmm4py, as well as for building and at runtime. Unfortunately, this is difficult to enforce for petsc4py when using pip because its build system always downloads and  builds against the latest compatible versions regardless of NumPy and Cython, regardless of what you already have installed in your Python environment (becaue of build isolation). This can be problematic on HPC systems with preinstalled versions of NumPy and other packages built against it. (That said, NumPy libraries at least are ABI-compatible across versions and unless you have very old versions of either NumPy and Cython none of this may pose an issue.) However, should you run into problems, you can ensure that petsc4py is built against your desired NumPy and Cython by turning off build isolation. First, install some required packages and (if they're not already installed) the NumPy and Cython versions you desire, for example:

```bash
pip install wheel setuptools numpy==1.24.4 cython==3.0.12
```

Then, install petsc4py *without* build isolation:

```bash
pip install petsc4py --no-build-isolation
```

The above trick also works for tmm4py, but its build system gives you much more control by letting you specify NumPy and Cython versions via environment variables, for example:

```bash
NUMPY_VER=1.24.4 CYTHON_VER=3.0.12 pip install tmm4py
```

#### Installation from source

You can also clone or download this repository and build TMM/tmm4py from source. In the following, it is assumed you already have PETSc (and petsc4py to use tmm4py) installed and have set PETSC_DIR/PETSC_ARCH. If you installed PETSc with pip these paths can be set with:

```bash
export PETSC_DIR=`python3 -c "from petsc import get_petsc_dir; print(get_petsc_dir())"`
export PETSC_ARCH=''
```
For building tmm4py, if you installed PETSc from source *and* installed petsc4py at the same time, the latter needs to be in your path:

```bash
export PYTHONPATH=$PYTHONPATH:$PETSC_DIR/$PETSC_ARCH/lib
```

(1) Build the TMM library. In the following, the path of the directory in which this README file is located is in the environment variable TOPDIR.

```bash
cd $TOPDIR  
make  
make install
```

By default this will install the TMM library and include files in a subdirectory 'TMM', i.e., $TOPDIR/TMM. If you want to install it in a different location, set the PREFIX variable, e.g.:

```bash
make install PREFIX=/home/userid/TMM
```

When installation is complete, the path to the installation location will be printed (the same as PREFIX if you provide it above). Set the environment variable TMM_DIR to point to this directory, e.g., for the default location:

```bash
export TMM_DIR=$TOPDIR/TMM
```

Note: You can delete the object files created above and in any of the subsequent examples with `make cleanall`

(2) Build the tmm4py extension module:

This step switches off build isolation (see above) so first make sure all the necessary Python modules are already installed. In addition to `petsc4py`, this includes `wheel`, `setuptools`, `numpy`and `cython` (for the latter two, the versions need to match the ones used to build petsc4py). Then:

```bash
cd $TOPDIR
make tmm4py
```

This will build and install tmm4py in site packages. (This step invokes pip to install tmm4py. By default it calls `pip3` but you can override this by setting the variable PIP: `make tmm4py PIP=pip`.)  You don't need to do anything else. However, if you prefer you can specify a different location with the PREFIX option (which can be the same or different from where you installed TMM in the previous step). For example, using the same location:

```bash
make tmm4py PREFIX=$TMM_DIR
```

You must then add this location to your PYTHONPATH so that the tmm4py dynamic library can be found:

```bash
export PYTHONPATH=$PYTHONPATH:$TMM_DIR
```

In addition to the dynamic libraries and header files, this will also install tmm.py, a generic driver script that can be used to run most models, whether implemented in Python or Fortran.

Note: The above make call invokes `pip3` to perform the installation. On some systems you may need to use `pip` instead. You can pass this to make with:

```bash
PIP=pip make tmm4py ...
```

#### Examples

Notes: 

(i) To run the examples in the instructions below, download this archive which contains all the necessary input data: https://drive.google.com/file/d/1ou1kamcI-ZE6BWHpnAs2LaSmO1uZeXyd/view?usp=sharing. Unpack and place the 'examples' directory in a convenient directory, the path to which is hereafter in the environment variable TOPDIR. (If you installed TMM from source then this can be the same as the TOPDIR referred to in the corresponding instructions.) Otherwise, you will need to download transport matrices and other data and scripts from https://sites.google.com/view/samarkhatiwala-research-tmm and generate the input data yourself (instructions below).

(ii) When using pip to install TMM/tmm4py, Makefile_driver and the models/ directory, which contains source code and various scripts for different biogeochemical models, are not installed. To run the examples, download and place them next to the examples/ directory (i.e., in $TOPDIR).

**Running the Python examples**

(1) To check if your installation is working, run the Age example:

```bash
cd $TOPDIR/examples/idealage  
cp -p $TOPDIR/models/idealage/code/python/Age.py .  
cp -p $TOPDIR/models/idealage/runscript/* .
mpiexec -np 4 python3 run_model_age_py.py
```

(If using a different MPI launcher modify the above command accordingly. For example, if you installed MPI via PETSc as in the instructions below, this would be $PETSC_DIR/$PETSC_ARCH/bin/mpiexec.)

(2) You can similarly run the other (pure) Python examples in the examples directory. These require the numba Python module which you can install with `pip install numba`.

MOPS:

```bash
cd $TOPDIR/examples/ironmops
cp -p $TOPDIR/models/ironmops/code/python/* .
cp -p $TOPDIR/models/ironmops/runscript/* .
```

and follow the instructions in 'runscript'

Pa/Th:

```bash
cd $TOPDIR/examples/protactinium_thorium
cp -p $TOPDIR/models/protactinium_thorium/code/python/* .
cp -p $TOPDIR/models/protactinium_thorium/input/MOBI_tracer_names.txt .
cp -p $TOPDIR/models/protactinium_thorium/input/get_tracer_names.sh .
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .
```

and follow the instructions in 'runscript'

**Running the Fortran examples and/or using the 'classic' TMM C driver:**

The following assumes you have installed PETSc with Fortran support.

(1) Regardless of how you installed PETSc/TMM, you will need to specify the paths to them by setting PETSC_DIR/PETSC_ARCH and TMM_DIR. 

If PETSc was installed with pip (either directly or as part of installing TMM):

```bash
export PETSC_DIR=`python3 -c "from petsc import get_petsc_dir; print(get_petsc_dir())"`
export PETSC_ARCH=''
```

If pip was used to install TMM:

```bash
export TMM_DIR=`python3 -c "from tmmlib import get_tmm_dir; print(get_tmm_dir())"`
```

(2) To use the 'classic' TMM C driver:

```bash
cd $TOPDIR
make -f Makefile_driver
```

This compiles tmm_main.c along with stubs for tmm_external*, tmm_monitor.c and tmm_misfit.c into a standalone executable 'tmm'. You can use this for the Age example (next step) as it reads in the forcing term from file.

(3) Age: To run the 'classic' TMM C driver, copy the executable $TOPDIR/tmm you built in Step 2 above:

```bash
cd $TOPDIR/examples/idealage
cp -p $TOPDIR/tmm .
cp -p $TOPDIR/models/idealage/runscript/* .
```

and follow the instructions in 'runscript'.

(4) Fortran models: Models written in Fortran can be run with either the 'classic' TMM C driver or via tmm4py. In both cases, you first compile the model as a dynamic library and then link it to either the TMM C driver or a Python gateway module (which makes the model callable from Python).

(i) MOPS:

First build MOPS into a dynamic library:

```bash
cd $TOPDIR/models/ironmops/code
make rmops
```

This will create libmops.so. (The Makefile targets several different versions of MOPS with different compiler directives; 'rmops' is the version with runoff, carbon and alkalinity switched on.) 

To run the model from Python, first ensure you have the Cython compiler installed:

```bash
pip install Cython
```

Then, compile the Python gateway interface into a Python extension module:

```bash
python3 setup.py build_ext --inplace
```

This compiles model_tmm_interface.pyx (which contains the gateway functions that call the actual external forcing and other C wrapper functions) into mops.cpython-XXX.so. Set PYTHONPATH so that this module can be found:

```bash
export PYTHONPATH=$PYTHONPATH:$TOPDIR/models/ironmops/code
```

To run the model with the TMM C driver, compile the C driver tmm_main.c and link to libmops.so to create a standalone executable tmmmops:

```bash
make -f Makefile_standalone tmmmops
make -f Makefile_standalone cleanall
```

To run via the Python driver:

```bash
cd $TOPDIR/examples/ironmops
cp -p $TOPDIR/models/ironmops/runscript/* .
```

and follow the instructions in 'runscript'.

To run via the C driver:

```bash
cd $TOPDIR/examples/ironmops
cp -p $TOPDIR/models/ironmops/runscript/* .
cp -p $TOPDIR/models/ironmops/code/tmmmops .
```

and follow the instructions in 'runscript'.

(ii) PaTh:

First download MOBI, which is embedded within the UVic ESCM model. You can get it from https://github.com/OSU-CEOAS-Schmittner/UVic2.9. Download the 'source' and 'updates' directories into a directory (e.g., UVic_ESCM). Set the environment variable UVICESCMROOT to point to this directory, e.g.:

```bash
export UVICESCMROOT=$HOME/UVic_ESCM
```

Next, build MOBI into a dynamic library:

```bash
cd $TOPDIR/models/mobi/code
cp -p MOBI_TMM_OPTIONS.h_PaTh MOBI_TMM_OPTIONS.h
make
```

This will create libmobi.so.

To run the model from Python, compile the Python gateway interface into a Python extension module:

```bash
python3 setup.py build_ext --inplace
```

This compiles model_tmm_interface.pyx (which contains the gateway functions that call the actual external forcing and other C wrapper functions) into mobi.cpython-XXX.so. Set PYTHONPATH so that this module can be found:

```bash
export PYTHONPATH=$PYTHONPATH:$TOPDIR/models/mobi/code
```

To run the model with the TMM C driver, compile the C driver tmm_main.c and link to libmobi.so to create a standalone executable tmmmobi:

```bash
make -f Makefile_standalone tmmmobi
make -f Makefile_standalone cleanall
```

To run with the Python driver:

```bash
cd $TOPDIR/examples/protactinium_thorium
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .
```

and follow the instructions in 'runscript'.

To run with the C driver:

```bash
cd $TOPDIR/examples/protactinium_thorium
cp -p $TOPDIR/models/protactinium_thorium/runscript/* .
cp -p $TOPDIR/models/mobi/code/tmmmobi .
```

and follow the instructions in 'runscript'.

#### Instructions for building MPI from source

(1) Download and unpack the source: https://www.mpich.org/

(2) Configure (with gcc and gfortran in this example):

```bash
cd mpich-4.3.0 # modify for the version you downloaded  
FC=gfortran  
./configure CC=gcc --prefix=/usr/local/mpich
```

(3) Build and install:

```bash
make  
sudo make install  
export PATH=$PATH:/usr/local/mpich/bin
```

#### Instructions for building PETSc from source (also installs MPI and other libraries)

###### Notes:

(i) See https://petsc.org/main/install/ for detailed installation instructions

(ii) On HPC systems, optimized BLAS/LAPACK (Intel MKL, Cray LibSci, etc) and MPI libraries are usually already available. If so, you should use these when configuring PETSc. In particular, use the compiler 'wrappers' mpicc and mpif90 instead of the actual compilers (the exception here is if you use the native compilers on Cray systems, in which case use the 'actual' compilers 'cc' and 'ftn' are in fact wrappers), and specify the location of the BLAS/LAPACK libraries with: --with-blaslapack-dir=...

(iii) If you plan to use tmm4py for 'production' calculations then configure it with the optimization options appropriate for your compiler (the default is to build PETSc in 'debugging' mode). For the GNU compiler, it would be something like: --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native'.

(1) Download a recent tarball from https://petsc.org/main/install/download/ and unpack.

(2) Set the environment variable PETSC_DIR to the path of the directory you just unpacked, e.g.:

```bash
export PETSC_DIR=$HOME/petsc-3.23.3
```

(3) Set the environment variable PETSC_ARCH to some convenient name for your installation (e.g., arch-darwin or arch-linux). You can have different installations under the same PETSC_DIR (e.g., one for debugging and one for production work) and switch between them by setting PETSC_ARCH accordingly, e.g.:

```bash
export PETSC_ARCH=arch-arm-darwin-opt
```

(4) Configure and build PETSc/petsc4py (the command below assumes your C and Fortran compilers are gcc  and gfortran, respectively, but you can change them, e.g., clang, icc, ifort etc. As mentioned above, replace these with mpicc and mpif90 if you already have MPI installed (in which case also remove the --download-mpich option). If you don't have or want to use a Fortran compiler, then set --with-fc=0.

```bash
cd $PETSC_DIR/  
./configure --with-cc=gcc --with-cxx=0 --with-fc=gfortran --download-f2cblaslapack --download-mpich --with-shared-libraries=1 --with-petsc4py=1
```

Then build PETSc by executing the make command printed out by configure above. It will be something like:

```bash
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
```

(5) When compilation is finished, PETSc will print out the path to where it installed the petsc4py library (typically it is $PETSC_DIR/$PETSC_ARCH/lib). Add this location to the PYTHONPATH environment variable, e.g.:

```bash
export PYTHONPATH=$PYTHONPATH:$PETSC_DIR/$PETSC_ARCH/lib
```

You're done!

Note: The above configure command will download and build BLAS/LAPACK and the MPICH MPI library for you. You just have to use the PETSc-installed MPI when launching your program with mpiexec or mpirun. PETSc will install them in: $PETSC_DIR/$PETSC_ARCH/bin. You can add this directory to your path or specify the full path, for example:

`$PETSC_DIR/$PETSC_ARCH/bin/mpiexec -np 4 ...`

#### Instructions for generating input data

To perform simulations you need transport matrices (TMs) and other forcing data. The relevant data and scripts to generate the input files can be downloaded from https://sites.google.com/view/samarkhatiwala-research-tmm.

(1) Download the Matlab scripts and add to your Matlab path. (Python scripts to perform similar operations are currently under development. Email me if you want them.)

(2) Download transport matrices and related data for the ocean model of your choice.  Currently, there are three configurations of MITgcm available online (and several others based on the UVic Earth System Model that I am happy to make available if you email me). For each, download the TMs and other associated data (e.g., MITgcm_ECCO.tar). Unpack. Make a note of the path to this directory (e.g., /mydisk/somewhere/MITgcm_ECCO). We will need it later. For some experiments you may also find it useful to download various miscellaneous data (and adjust paths accordingly in the provided Matlab scripts). 

(3) For each model, e.g., $TOPDIR/models/mops2.0, the input/ directory contains Matlab script(s) to generate input data (e.g., make_input_files_for_mops_model.m) and read model output (e.g., load_output.m).

Copy the matlab and run scripts to the directory in which you want to generate the files.

```bash
cp -p -R $TOPDIR/models/mops2.0/input/* .
```

Change the variable base_path at the very top of the make_input...m script to point to the top level of any of the transport matrix configurations you downloaded in #2 above. Execute the make_input* script (e.g., make_input_files_for_mops_model.m). It should generate all necessary input data. With luck! (If there is missing data email me for it.)