load basepath

spinupParamsFile='spinup_params';
jacobianParamsFile='jacobian_params';
pcDataFile='preconditioner_data';

useCoarseGrainedMatrix=0;
jacobianType=6; % '5' will call the 'q' function one vector at a time, '6' will call it with multiple vectors

% SET THE FOLLOWING THREE VARIABLES:
% locDir='/work_O1/smomw067/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
% remDir='/work_O1/smomw067/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
% remoteMachine='smomw067@rzcluster.rz.uni-kiel.de';
% locDir='/data2/spk/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
% remDir='/scratch1/fats/spk/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
% remoteMachine='spk@fats.ldeo.columbia.edu';
locDir='/data2/spk/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
remDir='/data3/dwang/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/KielBiogeochem_NPZDDOP_Matrix5_1';
remoteMachine='dwang@monk';

% UNCOMMENT ONE OF THE TWO SET OF FLAGS BELOW.
% (1) Options for running 'spinup model' on a remote machine:
sharedFileSystem=0;
singleMachine=0;
% (2) Options for running 'spinup model' on a local machine. In this case, Matlab is running on the 
% frontend of a cluster, so even though the model runs on a cluster, it is effectively local as 
% far as Matlab and the ClusterTools toolbox are concerned. 
% NOTE: IN THIS CASE, 'remoteMachine' ABOVE IS THE ADDRESS OF THE CLUSTER. ALSO, 'locDir' AND 
% 'remDir' MUST BE IDENTICAL.
% sharedFileSystem=1;
% singleMachine=1;

%%%
checkType=2; % 1 for waitscript, 2 for signal files
checkFile='newicsignalin';
checkToken='1';
spinupoptions=make_cluster_options('default',locDir,remDir,'exe',16,[],sharedFileSystem,singleMachine,checkType,checkFile,checkToken,'spinupcheckjob');
spinupoptions.frontend=remoteMachine;

% Options for running preconditioner on a remote machine (here assumed to be the same as that on which 
% the spinup model is run.
checkType=2; % 1 for waitscript, 2 for signal files
checkFile='newpcsignalin';
checkToken='1';
pcoptions=make_cluster_options('default',locDir,remDir,'exe',16,[],sharedFileSystem,singleMachine,checkType,checkFile,checkToken,'pccheckjob');
pcoptions.frontend=remoteMachine;

% Options for running Jacobian computation on local machine
jacobianjoboptions=[];

load(pcDataFile,'dt')

load(spinupParamsFile,'nb','numTracers')

% data for initial guess
% Specify PETSc files for variables NOT included in the MFNK spinup, but are part of the model.
% NOTE: you must manually copy these files to the machine(s) on which you plan to run the 
% spinup and Jacobian models.
OXY=300.0; % [mmol O2/m3]
if jacobianType==5
  writePetscBin('oxyini.petsc',repmat(OXY,[nb 1]))
  disp(['Set itout in the Jacobian model run script to: 1'])
elseif jacobianType==6
  load(pcDataFile,'sparsityFile')
  load(sparsityFile,'g')
  nqvecs=length(unique(g));
  writePetscBin('oxyini.petsc',repmat(OXY,[nb nqvecs+1]),1)
  disp(['Set itout in the Jacobian model run script to: ' int2str(nqvecs+1)])  
end

% variables that are part of the MFNK spinup:
PO4=2.17; % [mmol P/m3]
DOP=1e-4; % [mmol P/m3]
PHY=1e-4; % [mmol P/m3]
ZOO=1e-4; % [mmol P/m3]
DET=1e-4; % [mmol P/m3]
% initial guess
initial_guess=repmat([PO4 DOP PHY ZOO DET],[nb 1]); % the order of variables MUST match the filename set in prep_data_for_spinup.m
initial_guess=reshape(initial_guess,[prod(size(initial_guess)) 1]);
u0y=initial_guess(1:nb*numTracers);

phi=@(x)calc_solution_with_running_time_stepper(x,1,spinupoptions,45,spinupParamsFile);
func=@(x)phi(x)-x;

n=numTracers*nb;
uscale=ones(n,1); % uscale=u0y;

% traditional implicit PC
dti=86400*360;
itjac=1;
% function for biogeochemistry
qfunc=@(x)calc_q_with_running_time_stepper(x,itjac,jacobianjoboptions,[],jacobianParamsFile)/dt;
% preconditioner function
% funcp=@(x,y,z)apply_direct_inverse_with_implicit_model_pc(x,y,z,qfunc,dti,numTracers,jacobianType,pcDataFile,uscale,pcoptions,10);
funcp=@(x,y)apply_direct_inverse_with_implicit_model_pc(x,y,y,qfunc,dti,numTracers,jacobianType,pcDataFile,uscale,pcoptions,10);

funcw=@(x)check_positive(x,func); % function to monitor solution and trap NaN's

% We use KINSOL as the MFNK solver. Specify some data and options.
data.func=funcw;
data.uscale=[];
data.fscale=[];
data.funcp=funcp;
data.pcFlag=2; % 0=never recompute, 1=always recompute, 2=compute once at beginning, 3=skip computing first time

data.runTime=11*60*60; % run time in second (KINSol will pause after this time has elapsed until you hit enter)

strategy = 'LineSearch';
msbset = 1;

options = KINSetOptions('Verbose',true,...
                        'LinearSolver','GMRES',...
                        'KrylovMaxDim',40,...
                        'MaxNumRestarts',5,...
                        'MaxNumSetups', msbset,...
                        'EtaForm','Type2',...
                        'MaxNumSetups', msbset,...
                        'PrecSetupFn','kin_pset',...
                        'PrecSolveFn','kin_psol',...
                        'FuncNormTol',1e-4);

global firstTime recomputePC nonlinearIter funcEvalIter pcApplyIter pcSetupIter
global kinsol_history

error('err')

KINMalloc('kin_sys', n, options, data);  % allocate memory
[status, sol] = KINSol(u0y, strategy, ones(n,1), ones(n,1)); % call solver
s = KINGetStats; % get solver statistics
KINFree; % free memory

SOL=reshape(sol,[nb numTracers]); % final solution

% Now grid the solution
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(base_path,'Matrix5/Data','boxes');
[Ir,Ip,Irr,Ir_prev,Ip_prev,Ib,Ib_prev]=calc_profile_indices(boxFile);

% the gridded steady state solution (x,y,z)
[PO4g,xg,yg,zg]=grid_boxes3d(SOL(Irr,1),[],boxFile,gridFile);
[DOPg,xg,yg,zg]=grid_boxes3d(SOL(Irr,2),[],boxFile,gridFile);
[PHYg,xg,yg,zg]=grid_boxes3d(SOL(Irr,3),[],boxFile,gridFile);
[ZOOg,xg,yg,zg]=grid_boxes3d(SOL(Irr,4),[],boxFile,gridFile);
[DETg,xg,yg,zg]=grid_boxes3d(SOL(Irr,5),[],boxFile,gridFile);

% NOTE: The struct variable kinsol_history contains the entire convergence history of 
% the MFNK solver. Of particular interest are:
% kinsol_history.residual  <- history of residual
% kinsol_history.numfuncevals <- cumulative number of function evaluations at each nonlinear iteration
% kinsol_history.numpcapply <- cumulative number of PC applications at each nonlinear iteration
% Thus, to visualize the convergence history you can use the following matlab command:
%   semilogy(cell2mat(kinsol_history.numfuncevals),cell2mat(kinsol_history.residual),'r*-')
