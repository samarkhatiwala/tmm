% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicForcing=1
periodicMatrix=1

dt=43200; % time step to use

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

if rearrangeProfiles || bigMat
  load(profilesFile,'Ip_pre','Ir_pre','Ip_post','Ir_post','Irr')
  Ip=Ip_pre;
  Ir=Ir_pre;
end

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

if periodicForcing
  nm=12;
else
  nm=1;
end

% Use steady state T/S from GCM. Note we always load seasonal data here.
% load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
% load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
% Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
% Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);
% Use steady state T/S propagated from surface into interior using TMM.
load(fullfile(bgcDataPath,'Theta_bc'),'Tbc')
load(fullfile(bgcDataPath,'Salt_bc'),'Sbc')
Theta=gridToMatrix(Tbc,[],boxFile,gridFile);
Salt=gridToMatrix(Sbc,[],boxFile,gridFile);

% now take annual mean if necessary
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
end

Tss=Theta(Ib,:);
Sss=Salt(Ib,:);

Ts=Theta(Ii,:);
Ss=Salt(Ii,:);

atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
if ~periodicForcing
  atmospb=mean(atmospb,2);
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

if rearrangeProfiles
  error('ERROR: rearrangeProfiles must be set to 0!')
end  

if writeFiles
  calc_periodic_times_for_tmm('monthly-365-day year','periodic_times_365d.bin');
  calc_periodic_times_for_tmm('monthly-360-day year','periodic_times_360d.bin');  
% Transport matrices
  if writeTMs
%   Explicit transport matrix
	I=speye(nb,nb);
	if ~periodicMatrix
      disp('loading annual mean explicit TM')	
      load(explicitAnnualMeanMatrixFile,'Aexpms')	
	  if rearrangeProfiles
		Aexpms=Aexpms(Ir_pre,Ir_pre); % rearrange
	  end      
	  % make discrete
	  Aexpms=dt*Aexpms;
	  Aexpms=I+Aexpms;
	  [Ae1,Be,Ii]=split_transport_matrix(Aexpms,Ib);
	  writePetscBin('Ae1.petsc',Ae1,[],1)
	  writePetscBin('Be.petsc',Be,[],1)
	else
      % load each month from separate file
      disp('loading monthly mean explicit TMs')	      
	  for im=1:12 
		fn=[explicitMatrixFileBase '_' sprintf('%02d',im)];
		load(fn,'Aexp')
		if rearrangeProfiles
		  Aexp=Aexp(Ir_pre,Ir_pre); % rearrange
		end
		% make discrete
		Aexp=dt*Aexp;
		Aexp=I+Aexp;
		[Ae1,Be,Ii]=split_transport_matrix(Aexp,Ib);
		writePetscBin(['Ae1_' sprintf('%02d',im-1)],Ae1,[],1)
		writePetscBin(['Be_' sprintf('%02d',im-1)],Be,[],1)		
		clear Aexp Ae1 Be
	  end
	end
%   Implicit transport matrix
	if ~periodicMatrix
      disp('loading annual mean implicit TM')		
      load(implicitAnnualMeanMatrixFile,'Aimpms')
      if dtMultiple~=1
		if bigMat % big matrix. do it a block at a time.
		  for is=1:nbb % change time step multiple
			Aimpms(Ip_pre{is},Ip_pre{is})=Aimpms(Ip_pre{is},Ip_pre{is})^dtMultiple;
		  end
		else
		  Aimpms=Aimpms^dtMultiple;
		end  
	  end	
	  if rearrangeProfiles
		Aimpms=Aimpms(Ir_pre,Ir_pre); % rearrange
	  end
	  [Ai1,Bi,Ii]=split_transport_matrix(Aimpms,Ib);	  
	  writePetscBin('Ai1.petsc',Ai1,[],1)
	  writePetscBin('Bi.petsc',Bi,[],1)      
	else
	  % load each month from separate file
      disp('loading monthly mean implicit TMs')	      	  
	  for im=1:12
		fn=[implicitMatrixFileBase '_' sprintf('%02d',im)];		
		load(fn,'Aimp')
		if dtMultiple~=1
		  if bigMat % big matrix. do it a block at a time.		
			for is=1:nbb % change time step multiple
			  Aimp(Ip_pre{is},Ip_pre{is})=Aimp(Ip_pre{is},Ip_pre{is})^dtMultiple;
			end
		  else
			Aimp=Aimp^dtMultiple;		
		  end
		end  
		if rearrangeProfiles
		  Aimp=Aimp(Ir_pre,Ir_pre); % rearrange
		end
		[Ai1,Bi,Ii]=split_transport_matrix(Aimp,Ib);
		writePetscBin(['Ai1_' sprintf('%02d',im-1)],Ai1,[],1)
		writePetscBin(['Bi_' sprintf('%02d',im-1)],Bi,[],1)		  
		clear Aimp Ai1 Bi
	  end
	end
  end	  	  
% Initial conditions  
  for itr=1:numTracers
    gasId=trNames{itr};
    fn=[gasId 'ini.petsc'];
    writePetscBin(fn,TRini(:,itr))
  end  
% Surface forcing data
  if ~periodicForcing
	writePetscBin('atmosp.bin',atmospb)	
  else
    for im=1:nm
	  writePetscBin(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im))	  
	end
  end
  if ~periodicForcing
	writePetscBin('Tss.petsc',Tss)
	writePetscBin('Sss.petsc',Sss)
	writePetscBin('Ts.petsc',Ts)
	writePetscBin('Ss.petsc',Ss)
  else
    for im=1:nm
	  writePetscBin(['Tss_' sprintf('%02d',im-1)],Tss(:,im))
	  writePetscBin(['Sss_' sprintf('%02d',im-1)],Sss(:,im))
	  writePetscBin(['Ts_' sprintf('%02d',im-1)],Ts(:,im))
	  writePetscBin(['Ss_' sprintf('%02d',im-1)],Ss(:,im))
    end    
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
