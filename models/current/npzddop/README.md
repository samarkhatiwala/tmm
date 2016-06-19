This model requires the following input files:

% surface forcing
wind.bin
EmP.bin
atmosp.bin
fice.bin

% full 3-d fields
Ss.petsc
Ts.petsc

% grid data
latitude.bin (or swrad.bin if #ifdef READ_SWRAD)
drF.bin
dz.petsc
surface_volume_fraction.petsc
% + if using prognostic atmospheric CO2
dA.bin

% initial conditions
po4ini.petsc
dopini.petsc
oxyini.petsc
phyini.petsc
zooini.petsc
detini.petsc
% + if using carbon
dicini.petsc
alkini.petsc
% + if using prognostic atmospheric CO2
pCO2atm_ini.bin

