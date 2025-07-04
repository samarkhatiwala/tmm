% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

forcingType=1; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numForcingFields=12*1
%
matrixType=1; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numMatrices=12*335

Ttd0=1765 % Time origin for time dependent matrices, forcing etc 

dt=43200; % time step to use
daysPerYear=360;

rearrangeProfiles=0 % DON'T CHANGE!!
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0
writePCFiles=0

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

explicitMatrixFileBase=fullfile(base_path,explicitMatrixFileBase);
implicitMatrixFileBase=fullfile(base_path,implicitMatrixFileBase);

explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);

preconditionerMatrixFile=fullfile(base_path,preconditionerMatrixFile);

gcmDataPath=fullfile(base_path,'GCM');
bgcDataPath=fullfile(base_path,'BiogeochemData');
freshWaterForcingFile=fullfile(gcmDataPath,'FreshWaterForcing_gcm');
empFixFile=fullfile(gcmDataPath,empFixFile);

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','x','y','z','deltaT','gridType')

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
Ii=find(~ismember([1:nb]',Ib));
nbb=length(Ib);
nbi=length(Ii);

if rearrangeProfiles
  error('ERROR: rearrangeProfiles must be set to 0!')
end  

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

% Use steady state T/S from GCM. Note we always load seasonal data here.
load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);

if forcingType==0
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
end

Tss=Theta(Ib,:);
Sss=Salt(Ib,:);

Ts=Theta(Ii,:);
Ss=Salt(Ii,:);

atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
if forcingType==0
  atmospb=mean(atmospb,2);
elseif forcingType==2
  atmospb=repmat(atmospb,[1 numForcingFields/size(atmospb,2)]);
end

% Initial condition: initialize with solubility equilibrium concentration if 
% functions to compute it are available. You can download these from Roberta 
% Hamme's website: http://web.uvic.ca/~rhamme/download.html. Otherwise initialize 
% with zero concentration.
% Ne, Ar, Kr, Xe, N2, O2 or Ar36
trNames={'Ne','Ar','Kr','Xe','N2','O2','Ar36'};
numTracers=length(trNames);
TRini=repmat(0,[nbi numTracers]);
if exist('Nesol')==2
  disp('Initial conditions are being set to solubility equilibrium')
  Tmean=mean(Ts,2);
  Smean=mean(Ss,2);
  rho0=1024.5; % nominal density  
  for itr=1:numTracers
    gasId=trNames{itr};
    solFunc=str2func([gasId 'sol']);
    G_eq=rho0*solFunc(Smean,Tmean)/1e6; % mol/m^3/atm
    TRini(:,itr)=G_eq*1.0; % multiply by nominal atmospheric pressure of 1 atm to get equilibrium concentration
  end
else
  disp('No functions for computing solubility equilibrium found: Initial conditions set to zero')
  disp('You can download Matlab code for computing solubility equilibrium from Roberta Hamme''s website:')
  disp('http://web.uvic.ca/~rhamme/download.html')  
end

if writeFiles
  calc_periodic_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],['periodic_times_' num2str(daysPerYear) 'd.bin']);
  if matrixType==2
    calc_time_dependent_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],numMatrices/12,Ttd0,['matrix_times_' num2str(daysPerYear) 'd_' num2str(numMatrices) 'months.bin']);
  end
  if forcingType==2
    calc_time_dependent_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],numForcingFields/12,Ttd0,['forcing_times_' num2str(daysPerYear) 'd_' num2str(numForcingFields) 'months.bin']);
  end
  
% Transport matrices
  if writeTMs
	if matrixType==0
	  numMatrices=[];
    elseif matrixType==1
	  numMatrices=12;
	end  
    write_transport_matrices(base_path,dt,rearrangeProfiles,bigMat,useCoarseGrainedMatrix,matrixType,numMatrices,Ib)
  end

% Initial conditions  
  for itr=1:numTracers
    gasId=trNames{itr};
    fn=[gasId 'ini.petsc'];
    writePetscBin(fn,TRini(:,itr))
  end  
% Forcing data
  switch forcingType
    case 0
	  writePetscBin('Tss.petsc',Tss)
	  writePetscBin('Sss.petsc',Sss)
	  writePetscBin('atmosp.petsc',atmospb)
	  writePetscBin('Ts.petsc',Ts)
	  writePetscBin('Ss.petsc',Ss)
    case 1
	  for im=1:12
		writePetscBin(['Tss_' sprintf('%02d',im-1)],Tss(:,im))
		writePetscBin(['Sss_' sprintf('%02d',im-1)],Sss(:,im))
	    writePetscBin(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im))
		writePetscBin(['Ts_' sprintf('%02d',im-1)],Ts(:,im))
		writePetscBin(['Ss_' sprintf('%02d',im-1)],Ss(:,im))
	  end    
    case 2	  
%     use the first numForcingFields time slices of each field and bracket them with the first and last fields
	  tmp=Tss(:,[1 1:numForcingFields numForcingFields]);
	  writePetscBin('Tsstd.petsc',tmp,1)
	  tmp=Sss(:,[1 1:numForcingFields numForcingFields]);
	  writePetscBin('Ssstd.petsc',tmp,1)
	  tmp=atmospb(:,[1 1:numForcingFields numForcingFields]);
	  writePetscBin('atmosptd.petsc',tmp,1)
	  tmp=Ts(:,[1 1:numForcingFields numForcingFields]);
	  writePetscBin('Tstd.petsc',tmp,1)
	  tmp=Ss(:,[1 1:numForcingFields numForcingFields]);
	  writePetscBin('Sstd.petsc',tmp,1)
	otherwise
	  error('ERROR!: Unknown forcing type')    
  end    

% Grid data

% Profile data
  if rearrangeProfiles
    if ~useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ip);
	  gEndIndices=cellfun(@(x)x(end),Ip);
    else % useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ipcg);
	  gEndIndices=cellfun(@(x)x(end),Ipcg);
    end  
    write_binary('gStartIndices.bin',[length(gStartIndices);gStartIndices],'int')
    write_binary('gEndIndices.bin',[length(gEndIndices);gEndIndices],'int')
  end
end

if useCoarseGrainedMatrix
  numProfiles=nbbcg;
else  
  numProfiles=nbb;
end
disp(['Number of Profiles in this Configuration: ' int2str(numProfiles)])

if writePCFiles
  pc=load(preconditionerMatrixFile,'Aexpms');
  if rearrangeProfiles
    A=pc.Aexpms(Ir_pre,Ir_pre);
  else
    A=pc.Aexpms;
  end
  clear pc  
  A=split_transport_matrix(A,Ib);

  if useCoarseGrainedMatrix
    save pc_cg_data A nbbcg nbicg CGgrid CG Ipcg Ibcg dt
  else
    Ip=[];
    save pc_data A nbb nz nb nbi Ip Ib dt
  end
end
