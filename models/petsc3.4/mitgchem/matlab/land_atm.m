function [T,Cv,Cs,Cvs,pCo2a]=land_atm(Cvini,Csini,Cvsini,pCO2aini,tspan,E,D,Fo,usePreindustrial);


eps=1.e-3;   % tolerance
options=odeset('RelTol',eps,'AbsTol',eps); % set matlab ODE options

if nargin<9
 usePreindustrial=0;
end

deltaTsg=0;

if ~usePreindustrial
  Cini=[Cvini;Csini;Cvsini;pCO2aini];%KgC;KgC
else
  Cini=[Cvini;Csini];
end

[T,C]=ode45(@myfunc,tspan,Cini,options);  % call the ODEintegrator

Cv=C(:,1);
Cs=C(:,2);

if ~usePreindustrial
  Cvs=C(:,3);
  pCo2a=C(:,4);
else
  Cvs=Cvsini;  
  pCo2a=pCO2aini;  
end

% figure(5)
% plot(T,pCo2a)
% pause(1)

function S=myfunc(t,C)
%C=Cv,Cs,Cvs,pCo2a

if ~usePreindustrial
  [Sv,Ss,Svs,Sa]=landatmsource(C(1),C(2),C(3),C(4),E(t),D(t),Fo(t),deltaTsg);
  S=[Sv;Ss;Svs;Sa]; % sources % make sure this is a column vector!
else
  [Sv,Ss,Svs,Sa]=landatmsource(C(1),C(2),Cvsini,pCO2aini,E(t),D(t),Fo(t),deltaTsg);    
  S=[Sv;Ss];
end

end

end


