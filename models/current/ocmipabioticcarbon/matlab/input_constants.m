function [kd,kc,kM,deltaTo,Al,Ap,Ar,As,Q10p,Q10r,Q10s,uc,ut]=input_constants;
         %Same Values and Units as in Eliseev and Mokhov, 2007
%Cini=kgC or ppmv
%kd=adimensional
%kc=kM=ppmv
%deltaTo=K
%Al=Ap=Ar=As=1/year
%Q10p=Q10r=Q10s=adimensional
% uc=kgC/ppmv
% ut=kgC/K

% Definition of parameters 
kc=29;%ppmv
kM=150;%ppmv

deltaTo=10;% degrees K


Al=1/(11); %1/year
Ar=Al; %1/year
Ap=1/(3.4); %1/year
As=1/(30.);%1/year

Q10p=1.5; %adimensional numbers
Q10r=2.15;% // 
Q10s=2.4;%//
kd=0.27;%adimensional

uc=1.31e12;%KgC/ppmv
ut=0.33e14;%KgC/K
% uc=0;
% ut=0;