This model requires the following input files:

% surface forcing
Alk.bin or Alk_**
PO4.bin or PO4_**
SiO2.bin or SiO2_**
atmosp.bin or atmosp_**
% if using EmP
EmP.bin or EmP_**

% if annual mean forcing
Vgas.bin
% if periodic forcing
fice_**
% if periodic forcing and gasExchangeType==1 (OCMIP gas exchange coefficient)
xkw_**
% if periodic forcing and gasExchangeType==2 (gas exchange computed online from winds)
uwind_**
vwind_**

% if using prescribed atmospheric CO2 history
TpCO2.bin
pCO2atm.bin
% if using C14 and NOT using prognostic atmospheric model
DC14atm.bin
% if using C14 AND using prescribed atmospheric CO2 history
TC14atm.bin
% if using emissions (must also use prognostic atmospheric CO2)
Tem.bin
fossil_fuel_emissions.bin
land_use_emissions.bin

% full 3-d fields
Ss.petsc or Ss_**
Ts.petsc or Ts_**

% grid data
% if using EmP
surface_volume_fraction.petsc
% if using prognostic atmospheric CO2
dA.bin


% initial conditions
dicini.petsc
% if using C14
dic14ini.petsc
% if using prognostic atmospheric CO2
pCO2atm_ini.bin
% if using land model (must also use prognostic atmospheric CO2)
land_ini.bin


