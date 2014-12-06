% SPK: the following computes the steady state solution for the 
%    land model using a fixed preindustrial pCO2a

pCO2a_pre=278;

%First integrate Cv and Cs to equilibium assuming fixed pCO2a and Cvs
E=@(t)zeros(size(t));
D=@(t)zeros(size(t));
Fo=@(t,pCO2a)zeros(size(t));

% some initial conditions
Cvini=0.55e15;%KgC
Csini=1.5e15;%KgC
Cvsini=0.55e15;%KgC
Cvini=Cvini*1e3*1e-15;%PgC
Csini=Csini*1e3*1e-15;%PgC
Cvsini=Cvsini*1e3*1e-15;%PgC
pCO2aini=pCO2a_pre;%ppmv

func=@(C)steadystatefunc2(C,Cvsini,pCO2aini,E,D,Fo);
% call nonlinear solver
sol=nsoli([0;0],func,[1e-9 1e-9]);

% this is now the steady state solution
Cvpre=sol(1);
Cspre=sol(2);
Cvspre=Cvsini;

% SPK: to do a forward integration of the combined land-atmos model, the 
% following command can be used. 
% NOTE: this will currently NOT work since you dont't have all the input data
[T,Cv,Cs,Cvs,pCO2a]=land_atm(Cvini,Csini,Cvsini,pCO2aini,tspan,E,D,Fo);
% Input arguments:
% Cvini, Csini, and Cvsini are initial conditions for the land model
% pCO2aini is the initial atmospheric pCO2
% tspan is the length of time over which you want to perform the integration
% E, D, and Fo are function handles for the fossil fuel emissions, land use emissions, 
% and ocean uptake, respectively. Here is an example of how you could define these:
% E=@(t)interp1(Tem,Ev,t);
% D=@(t)interp1(Tem,Dv,t);
% Fo=@(t,pCO2a)interp1(Toc,Fov,t);
% Output arguments:
% T is the vector of times at which the output is defined
% Cv,Cs,Cvs,pCO2a are vectors of those variables definied at times in T
