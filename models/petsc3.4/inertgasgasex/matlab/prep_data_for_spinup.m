%% Set some problem-specific variables %%%
useCoarseGrainedMatrix=0;

% Specify number of tracers in the spinup problem
numTracers=1;

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'spinup model'
iniSpinupFileName={'trstart.petsc'};
outSpinupFileName={'trend.petsc'};

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'Jacobian model'
iniJacobianFileName={'trstartjac.petsc'};
outJacobianFileName={'trq.petsc'};

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
% nnzero=nb*nz; % upper bound
% nnzero=0; % better estimate
% for is=1:nbb % loop over each surface point
%   Ipl=Ip{is}; % indices for local profile (globally indexed)
%   nnzero=nnzero+length(Ipl)^2;
% end

% S=spalloc(nb,nb,nnzero);
% for is=1:nbb % loop over each surface point
%   Ipl=Ip{is}; % indices for local profile (globally indexed)
%   S(Ipl,Ipl)=1;       
% end

Isp=zeros(nbb,1);
for is=1:nbb % loop over each surface point
  Ipl=Ip{is}; % indices for local profile (globally indexed)
  Isp(is)=Ipl(1); % index of surface point
end
S=sparse(Isp,Isp,ones(nbb,1),nb,nb);

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
