%% Set some problem-specific variables %%%
useCoarseGrainedMatrix=1;

% Specify number of tracers in the spinup problem
numTracers=5;

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'spinup model'
iniSpinupFileName={'po4start.petsc','dopstart.petsc','phystart.petsc','zoostart.petsc','detstart.petsc'};
outSpinupFileName={'po4end.petsc','dopend.petsc','phyend.petsc','zooend.petsc','detend.petsc'};

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'Jacobian model'
iniJacobianFileName={'po4startjac.petsc','dopstartjac.petsc','phystartjac.petsc','zoostartjac.petsc','detstartjac.petsc'};
outJacobianFileName={'po4q.petsc','dopq.petsc','phyq.petsc','zooq.petsc','detq.petsc'};

% Specify Matlab data file names (best not to change these)
sparsityFile='sparsity_data';
spinupParamsFile='spinup_params';
jacobianParamsFile='jacobian_params';
pcDataFile='preconditioner_data';

%% DON'T MODIFY BELOW THIS LINE %%

if useCoarseGrainedMatrix
  load pc_cg_data A nbbcg CGgrid CG Ipcg Ibcg dt  
  nz=CGgrid.nz;  
  nb=CG.nb;
  nbb=nbbcg;
  Ib=Ibcg;
  Ip=Ipcg;
else
  load pc_data A nbb nz nb Ip Ib dt
end

save(pcDataFile,'A','Ip','nb','nz','nbb','dt','sparsityFile')

clear A

% Sparsity pattern
nnzero=nb*nz; % upper bound
nnzero=0; % better estimate
for is=1:nbb % loop over each surface point
  Ipl=Ip{is}; % indices for local profile (globally indexed)
  nnzero=nnzero+length(Ipl)^2;
end

S=spalloc(nb,nb,nnzero);
for is=1:nbb % loop over each surface point
  Ipl=Ip{is}; % indices for local profile (globally indexed)
  S(Ipl,Ipl)=1;       
end

Sstr='[';
Srow=repmat('S ',[1 numTracers]);
Srow=Srow(1:end-1);
for itr=1:numTracers-1
  Sstr=[Sstr Srow ';'];
end
Sstr=[Sstr Srow ']'];

eval(['S=' Sstr ';']);

% g=zeros(nb,numTracers);                 
% for is=1:nbb
%   Ipl=Ip{is};             
%   nzloc=length(Ipl);                   
%   ii=[1:numTracers*nzloc]';                
%   g(Ipl,:)=reshape(ii,[nzloc numTracers]);
% end
% g=reshape(g,[nb*numTracers 1]);        

g=colgroup_fast(S);

save(sparsityFile,'S','g')

% clear S g

% Spinup params
iniFileName=iniSpinupFileName;
outFileName=outSpinupFileName;

save(spinupParamsFile,'numTracers','nb','iniFileName','outFileName')

% Jacobian params
iniFileName=iniJacobianFileName;
outFileName=outJacobianFileName;

save(jacobianParamsFile,'numTracers','nb','iniFileName','outFileName')
