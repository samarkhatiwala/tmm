function [Cvpre,Cspre,Cvspre]=calc_steadystate_land(pCO2a_pre)

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
