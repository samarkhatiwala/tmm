load basepath

spinupParamsFile='spinup_params';
jacobianParamsFile='jacobian_params';
pcDataFile='preconditioner_data';

useCoarseGrainedMatrix=0;
jacobianType=6; % '5' will call the 'q' function one vector at a time, '6' will call it with multiple vectors

% SET THE FOLLOWING THREE VARIABLES:
% locDir='/work_O1/smomw067/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/InertGas_Matrix5_1';
% remDir='/work_O1/smomw067/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/InertGas_Matrix5_1';
remoteMachine='spk@darwin.ldeo.columbia.edu';
locDir='/home/spk/MITgcm_c57g_post/spk/MIT_Matrix_ECCO/InertGas_Matrix1_1';
remDir='/home/spk/MITgcm_c57g_post/spk/MIT_Matrix_ECCO/InertGas_Matrix1_1';
% locDir='/data2/spk/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/InertGas_Matrix5_1';
% remDir='/data2/spk/MITgcm_c57g_post/spk/MIT_Matrix_Global_2.8deg/InertGas_Matrix5_1';
% remoteMachine='localhost';

% UNCOMMENT ONE OF THE TWO SET OF FLAGS BELOW.
% (1) Options for running 'spinup model' on a remote machine:
% sharedFileSystem=0;
% singleMachine=0;
% (2) Options for running 'spinup model' on a local machine. In this case, Matlab is running on the 
% frontend of a cluster, so even though the model runs on a cluster, it is effectively local as 
% far as Matlab and the ClusterTools toolbox are concerned. 
% NOTE: IN THIS CASE, 'remoteMachine' ABOVE IS THE ADDRESS OF THE CLUSTER. ALSO, 'locDir' AND 
% 'remDir' MUST BE IDENTICAL.
sharedFileSystem=1;
singleMachine=1;

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
if jacobianType==5
  disp(['Set write_steps in the Jacobian model run script to: 1'])
elseif jacobianType==6
  load(pcDataFile,'sparsityFile')
  load(sparsityFile,'g')
  gu=unique(g);
  gu=gu(gu~=0);
  nqvecs=length(gu);
  disp(['Set write_steps in the Jacobian model run script to: ' int2str(nqvecs+1)])  
end

% variables that are part of the MFNK spinup:
TR=0;
% initial guess
initial_guess=repmat([TR],[nb 1]); % the order of variables MUST match the filename set in prep_data_for_spinup.m
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
% B0inv=@(x,y,z)apply_direct_inverse_with_implicit_model_pc(x,y,z,qfunc,dti,numTracers,jacobianType,pcDataFile,uscale,pcoptions,10);
B0inv=@(x,y)apply_direct_inverse_with_implicit_model_pc(x,y,y,qfunc,dti,numTracers,jacobianType,pcDataFile,uscale,pcoptions,10);

funcw=@(x)check_positive(x,func); % function to monitor solution and trap NaN's

gridFile=fullfile(base_path,'grid');
load(gridFile,'x','y','z')
boxFile=fullfile(base_path,'Matrix1/Data','boxes');
[Ir,Ip,Irr,Ir_prev,Ip_prev,Ib,Ib_prev]=calc_profile_indices(boxFile);

spinupTemplateFile='runpbs_spinup_darwin.signalfiles_template';
jacobianTemplateFile='runpbs_jacobian_darwin.waitscript_template';
timeAvgTemplateFile='runpbs_timeavg_darwin_template';
runFile='start_model';
jacFile='start_jacobian';

scaleFactor={[1 1 1],[0 0 1],[1 0 1],[0 1 1]};
tracers={'Ne','Ar','Kr','Xe','N2','Ar36','O2'};

numTr=length(tracers);
  
absTol=1e-12;

numSc=length(scaleFactor);

for itr=2:numTr

  gasId=itr;
  tracerID=tracers{itr};
  
  global broyden_history broyden_data 
  
  broyden_data.pcFlag=3;
  broyden_data.updatePCWithGMRES=1;
  broyden_data.B0inv=B0inv;
  broyden_data.maxsaveiter=20;
  broyden_data.nfreeze=-1;
  broyden_data.doLineSearch=1;
  broyden_data.historyfile='broyden_history_dump';
  runTime=11*60*60; % run time in second
  
  myfunc=@(x)eval_function_with_pause(x,func,runTime);

  for isc=1:numSc  
    Finj=scaleFactor{isc}(1);
    Fex=scaleFactor{isc}(2);
    Fgas=scaleFactor{isc}(3);    
  
    evalExternalCommand(['cp ' spinupTemplateFile ' ' runFile]);
    srp('GASID',sprintf('%d',gasId),runFile)
    srp('FINJ',sprintf('%1.1f',Finj),runFile)
    srp('FEX',sprintf('%1.1f',Fex),runFile)
    srp('FGAS',sprintf('%1.1f',Fgas),runFile)    
    evalExternalCommand(['chmod +x ' runFile]);

    evalExternalCommand(['cp ' jacobianTemplateFile ' ' jacFile]);
    srp('GASID',sprintf('%d',gasId),jacFile)
    srp('FINJ',sprintf('%1.1f',Finj),jacFile)
    srp('FEX',sprintf('%1.1f',Fex),jacFile)
    srp('FGAS',sprintf('%1.1f',Fgas),jacFile)    
    evalExternalCommand(['chmod +x ' jacFile]);

    input('Now START the startdirectimplicitpc script and hit enter ...','s')
    input(['Now START the ' jacFile ' script and hit enter ...'],'s')    
    input(['Now START the ' runFile ' script and hit enter ...'],'s')
    
    [sol, it_hist, ierr] = mfnk_broyden(u0y,myfunc,absTol);

    input(['Now QUIT the ' runFile ' script and hit enter ...'],'s')
    input('Now QUIT the startdirectimplicitpc script and hit enter ...','s')
    input(['Now QUIT the ' jacFile ' script and hit enter ...'],'s')    

%   calculate diagnostics    
    writePetscBin('trsteady.petsc',sol)    

    evalExternalCommand(['cp ' timeAvgTemplateFile ' ' runFile]);
    srp('GASID',sprintf('%d',gasId),runFile)
    srp('FINJ',sprintf('%1.1f',Finj),runFile)
    srp('FEX',sprintf('%1.1f',Fex),runFile)
    srp('FGAS',sprintf('%1.1f',Fgas),runFile)    
    evalExternalCommand(['chmod +x ' runFile]);
    
    input(['Now START the ' runFile ' script and hit enter WHEN FINISHED'],'s')

%   read diagnostics
    TRmm=readPetscBinVec('trmm.petsc',-1);
    TRmmg=grid_boxes3d(TRmm(Irr,:),[],boxFile,gridFile);
    TReqmm=readPetscBinVec('TReqdiagavg.petsc',-1);
    TReqmmg=grid_boxes3d(TReqmm(Irr,:),[],boxFile,gridFile);
    TRsatanommm=readPetscBinVec('TRsatanomdiagavg.petsc',-1);             
    TRsatanommmg=grid_boxes3d(TRsatanommm(Irr,:),[],boxFile,gridFile);
    
    TRavgg=mean(TRmmg,4);
    TReqavgg=mean(TReqmmg,4);
    TRsatanomavgg=mean(TRsatanommmg,4);

    fn=[tracerID '_' sprintf('%d',Finj) '_' sprintf('%d',Fex) '_' sprintf('%d',Fgas)];    
    save([fn '_mfnk_broyden_dump'],'-v7.3')
    save(fn,'x','y','z','TRmmg','TRsatanommmg','TReqmmg')    
    
    evalExternalCommand(['mkdir -p ' fn]);
    evalExternalCommand(['mv exfluxavg.bin injfluxavg.bin output_time_timeavg.txt TReqdiagavg.petsc TRsatanomdiagavg.petsc gasfluxavg.bin log_timeavg pickup_timeavg.petsc trmm.petsc ' fn '/']);
    
  end
end  
  
