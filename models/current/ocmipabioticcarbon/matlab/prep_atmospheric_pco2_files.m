emissionScenario='B2'

oceanCarbonBasePath='/data2/spk/OceanCarbon';
ocmip_path=fullfile(oceanCarbonBasePath,'OCMIP2');

% Load atmospheric pCO2 history. 
% Smoothed OCMIP data (from Sara Fletcher)
pCO2_atm=load(fullfile(oceanCarbonBasePath,'OCMIP2','splco2_cis92a.dat'));
% Units: pCO2 contains the CO2 mole fraction (mixing ratio) in dry air in ppm. 
% This is numerically the same as the partial pressure in uatm of CO2 in dry air
% at total pressure of 1 atm.
Tco2=pCO2_atm(:,1);
pCO2atm=pCO2_atm(:,2);
% future pCO2 from ISAM model
[tmp,pco2f]=hdrload(fullfile(oceanCarbonBasePath,'tar-isam.txt'));
Tf=pco2f(:,1);
% 1     2    3    4    5   6   7   8    9    10   11     12     13
% Year  A1B  A1T  A1FI  A2  B1  B2  A1p  A2p  B1p  B2p  IS92a  IS92a/SAR

switch emissionScenario
  case 'A1B'
    ic=2;
  case 'A1T'
    ic=3;
  case 'A1FI'
    ic=4;
  case 'A2'
    ic=5;
  case 'B1'
    ic=6;
  case 'B2'
    ic=7;
  case 'A1p'
    ic=8;
  case 'A2p'
    ic=9;
  case 'B1p'  
    ic=10;
  case 'B2p'
    ic=11;
  case 'IS92a'
    ic=12;
  case 'IS92a/SAR'
    ic=13;
  otherwise
    error('Unknown emission scenario')
end

pco2f=pco2f(:,ic); % A1B                                                                                                                 

k1=find(Tco2==2000);
k2=find(Tf==2000);
pco2dat=[pCO2atm(1:k1-1);pco2f(k2:end)];
Tdat=[Tco2(1:k1-1);Tf(k2:end)];
Tco2=[Tco2(1:k1-1):1:Tf(end)]';
pCO2atm=interp1(Tdat,pco2dat,Tco2);

write_binary('TpCO2.bin',Tco2,'real*8')
write_binary('pCO2atm.bin',pCO2atm,'real*8')

