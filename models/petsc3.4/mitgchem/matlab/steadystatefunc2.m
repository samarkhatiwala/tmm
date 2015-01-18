function F=steadystatefunc2(C,Cvsini,pCO2aini,E,D,Fo);

deltaTsg=0;
t=0;

[Sv,Ss,Svs,Sa]=landatmsource(C(1),C(2),Cvsini,pCO2aini,E(t),D(t),Fo(t),deltaTsg);    
F=[Sv;Ss];

