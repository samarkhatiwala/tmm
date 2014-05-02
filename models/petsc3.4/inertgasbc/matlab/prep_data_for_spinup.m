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
  nbi=nbicg;
  Ib=Ibcg;
  Ip=Ipcg;
else
  load pc_data A nbb nz nb nbi Ip Ib dt
end

nb=nbi;

save(pcDataFile,'A','Ip','nb','nz','nbb','dt','sparsityFile')

clear A

S=[];
g=[];

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
