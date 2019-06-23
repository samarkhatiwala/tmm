This is the TMM interface to version 2.0 of MOBI (Model of Ocean Biogeochemistry and Isotopes; 
http://people.oregonstate.edu/~schmita2/Models/MOBI/index.html). 

(1) To use this code, first download the compatible versions of UVic ESCM 2.9 and MOBI 2.0 
distributed from here: https://github.com/OSU-CEOAS-Schmittner/UVic2.9/wiki

(2) Next, make a directory in which you want to compile and run TMM-MOBI. We will call this 
TEST below.

(3) Copy the contents (including subdirectories) of mobi2.0/matlab/ and mobi2.0/runscripts/ 
to TEST/. 

(4) Copy mobi2.0/src/Makefile and mobi2.0/src/MOBI_TMM_OPTIONS.h to TEST/ and follow the 
instructions at the top of Makefile for setting paths.

(5) If you wish to make changes to the MOBI or interface code copy the relevant source files 
to TEST/ and make your changes.

(6) Edit namelist parameters in control.in as necessary.

(7) Edit the C preprocessor options in MOBI_TMM_OPTIONS.h as necessary. 

(7) Build with:
make cleanall
make smallf
make tmmmobiwrite
make tmmmobi

NOTE: the Makefile endeavors to detect your compiler so that the correct compiler options 
can be set. But that can sometimes fail. If you want to make absolutely certain, you can 
specify your compiler, e.g.,
make COMPILER=gfortran tmmmobi

(8) To generate the names of tracers being simulated (as per the options set in 
MOBI_TMM_OPTIONS.h) and their initial conditions (as would be set by UVic/MOBI 
run in online mode), execute:
./tmmmobiwrite

This will write out MOBI_tracer_names.txt containing the tracer names, and a corresponding 
set of *.dat text files (e.g., dic.dat) containing the vertical profile of tracer concentration. 
These names and (horizontally uniform) profiles are used by make_input_files_for_mobi_model.m 
to generate initial condition files (e.g., dicini.petsc).

Make a directory InitialConditionProfiles and copy the *.dat files into it. (By default 
make_input_files_for_mobi_model.m will look in this directory for the *.dat files.)

NOTE-1: There are already copies of MOBI_tracer_names.txt and InitialConditionProfiles/ 
in matlab/ but I recommend you still go through this step to generate tracer names consistent 
with your version of UVic/MOBI and MOBI_TMM_OPTIONS.h.

NOTE-2: The order of tracers in MOBI_tracer_names.txt is the same as that in which you 
specify initial condition, output and time average files in the TMM-MOBI run script. 
You can reduce the tedium (and chances of making a mistake) of generating the correct 
sequence to specify in the run script by executing this on the command line: 
printf "%sini.petsc," `tr '\n' ' ' < MOBI_tracer_names.txt` | sed 's/,*$//g'
(The tr replaces newlines with spaces to put all tracer names on a single line, 
the printf appends "ini.petsc," to each tracer name, and the sed strips out the 
final comma.)

NOTE-3: tmmmobiwrite uses UVic routines where initial condition profiles are set. The code 
for this is hardwired for the standard UVic grid with 19 levels. If you're using TMs from 
a different ocean model you will need to modify make_input_files_for_mobi_model.m accordingly 
to specify initial conditions corresponding to your grid. You will also need to change the 
number of levels (km) in size.h.

(9) Set the path at the top of make_input_files_for_mobi_model.m and execute in Matlab. 
(See NOTE-3 in #8.)

(10) Edit runscript as necessary (see NOTE-2 in #8) and run the model.

(11) Set the path at the top of the load*.m scripts to read in the output in Matlab.
