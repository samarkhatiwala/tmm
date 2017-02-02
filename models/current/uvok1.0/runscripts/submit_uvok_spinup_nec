#!/bin/ksh 
#PBS -T intmpi
#PBS -q clmedium                # Batchklasse   
#PBS -l cpunum_job=16           # CPUs pro Knoten: max. 16 
#PBS -b 10                       # number of nodes
#PBS -l cputim_job=1600:00:00     # CPU-Time=Coreanzahl*Walltime
#PBS -l elapstim_req=10:00:00   # Walltime (Verweilzeit)
#PBS -l memsz_job=10gb          # angeforderter Hauptspeicher (max. 128GB 256 in clbigmem)
#PBS -N uvok                    # Name des Batchjobs
#PBS -o uvok.out                # Name der Standardausgabedatei 
#PBS -j o                       # Zusammenfügen von Standardausgabe- und Fehlerdatei  

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

. /opt/modules/Modules/3.2.6/init/ksh   # init for clmedium

module load intel
module load intelmpi
module load petsc-3.6.1-intel

mpirun $NQSII_MPIOPTS -np 160 ./tmmuvok \
  -numtracers 10 \
  -i dicini.petsc,c14ini.petsc,alkini.petsc,o2ini.petsc,po4ini.petsc,phytini.petsc,zoopini.petsc,detrini.petsc,no3ini.petsc,diazini.petsc \
  -me Ae \
  -mi Ai \
  -t0 0.0 -iter0 0 \
  -deltat_clock 0.0009132420091324 \
  -max_steps 5475000 \
  -write_steps 547500 \
  -o dic.petsc,c14.petsc,alk.petsc,o2.petsc,po4.petsc,phyt.petsc,zoop.petsc,detr.petsc,no3.petsc,diaz.petsc \
  -external_forcing \
  -use_profiles \
  -biogeochem_deltat 28800.0 \
  -days_per_year 365.0 \
  -periodic_matrix \
  -matrix_cycle_period 1.0 -matrix_num_per_period 12 -matrix_periodic_times_file periodic_times_365d.bin \
  -periodic_biogeochem_forcing \
  -periodic_biogeochem_cycle_period 1.0 -periodic_biogeochem_num_per_period 12 -periodic_biogeochem_periodic_times_file periodic_times_365d.bin \
  -time_avg -avg_start_time_step 5473906 -avg_time_steps 93,84,93,90,93,90,93,93,90,93,90,93 \
  -avg_files dicmm.petsc,c14mm.petsc,alkmm.petsc,o2mm.petsc,po4mm.petsc,phytmm.petsc,zoopmm.petsc,detrmm.petsc,no3mm.petsc,diazmm.petsc \
  -calc_diagnostics -diag_start_time_step 5473906 -diag_time_steps 93,84,93,90,93,90,93,93,90,93,90,93 \
  > log

# driver options:
#  I/O:
#   -pickup_out pickup.petsc \
#   -time_file output_time.txt \
#  circulation-specific options:
#   add this option for configurations in which the forcing is rescaled because of a time-varying free surface:
#    -rescale_forcing_file Rfs \
# model-specific options:
