function [Sv,Ss,Svs,Sa]=landatmsource(Cv,Cs,Cvs,pCo2a,Ef,D,Fo,deltaTsg)

%Transform D and Ef in KgC/y in order to be consisten with E&M units (kgC/y)
% D=1e12*D;
% Ef=1e12*Ef;
% Fo=1e12*Fo;
    
%type help input_constants for units
[kd,kc,kM,deltaTo,Al,Ap,Ar,As,Q10p,Q10r,Q10s,uc,ut]=input_constants;
[mass_kg,ppmtoPgC]=atmospheric_gas_mass(12,1);%ppm to PgC
%ppmtoKgC=ppmtoPgC*1e12; %PgC to KgC

%define gf=fertility factor (adimensional)
if pCo2a<kc
 gf=0;
else
 gf=(pCo2a-kc)./(kM+pCo2a-kc);
end


% 
P=Ap.*gf.*Cvs.*(Q10p).^(deltaTsg./deltaTo);%carbon production rate of photosynthesis
Rp=Ar.*Cv.*(Q10r).^(deltaTsg./deltaTo);%Autotrophic (biota) respiration rate
L=Al.*Cv;%Littefall
Rs=As.*Cs.*(Q10s).^(deltaTsg./deltaTo);%Heterotrophic (soil) respiration rate
NPP=P-Rp;%Terrestrial vegetation net primary production

%
Fl=NPP-Rs;



%sour   ces/sinks
Sv=NPP-L-D;% carbon terrestrial vegetation stock source
Ss=L-Rs;% carbon terrestrial soil stock source
Svs=-kd.*D;

Sa=((Ef+D)-(Fl+Fo))./(ppmtoPgC);

