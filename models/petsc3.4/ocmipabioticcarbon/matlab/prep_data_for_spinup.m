%% Set some problem-specific variables %%%
useCoarseGrainedMatrix=0;

doC14calc=1

% Specify number of tracers in the spinup problem
numTracers=1; % DON'T CHANGE THIS!

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'spinup model'
if doC14calc
  iniSpinupFileName={'dic14start.petsc'};
  outSpinupFileName={'dic14end.petsc'};
else
  iniSpinupFileName={'dicstart.petsc'};
  outSpinupFileName={'dicend.petsc'};
end

% Specify PETSc filenames for tracers included in the MFNK spinup and required by the 'Jacobian model'
if doC14calc
  iniJacobianFileName={'dic14startjac.petsc'};
  outJacobianFileName={'dic14q.petsc'};
else
  iniJacobianFileName={'dicstartjac.petsc'};
  outJacobianFileName={'dicq.petsc'};
end

% Specify Matlab data file names (best not to change these)
if doC14calc
  sparsityFile='sparsity_data_c14';
  spinupParamsFile='spinup_params_c14';
  jacobianParamsFile='jacobian_params_c14';
  pcDataFile='preconditioner_data_c14';
else
  sparsityFile='sparsity_data';
  spinupParamsFile='spinup_params';
  jacobianParamsFile='jacobian_params';
  pcDataFile='preconditioner_data';
end

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
if numTracers==1
  if doC14calc
    S=sparse([1:nb]',[1:nb]',ones(nb,1),nb,nb); % DIC14  
  else
    S=sparse(Ib,Ib,ones(size(Ib)),nb,nb); % DIC
  end
  g=repmat(1,[nb 1]);  
else
  error('Simultaneous spinup of both DIC and C14 NOT supported! Must change the above code!!')
% DIC does not depend on DIC14
  i1=[1:nb]';
  j1=[1:nb]';
  ii=[Ib;nb+Ib;nb+i1];
  jj=[Ib;Ib;nb+j1];
  S=sparse(ii,jj,ones(size(ii)),nb*numTracers,nb*numTracers); 
  g=ones(nb*numTracers,1);
  g(nb+Ib)=2;
end
% S=speye(nb,nb);
% if numTracers==2
%   S=[S S;S S];
%   S(1:nb,nb+1:end)=0; % DIC does not depend on DIC14 
% end
% g=colgroup(S);

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
