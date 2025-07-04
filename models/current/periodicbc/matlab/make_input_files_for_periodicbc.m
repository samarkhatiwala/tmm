% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicMatrix=1

dt=43200; % time step to use

rearrangeProfiles=0 % DON'T CHANGE!!
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

explicitMatrixFileBase=fullfile(base_path,explicitMatrixFileBase);
implicitMatrixFileBase=fullfile(base_path,implicitMatrixFileBase);

explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);

gcmDataPath=fullfile(base_path,'GCM');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','deltaT','gridType')

load(boxFile,'izBox','nb')

if rearrangeProfiles
  error('ERROR: rearrangeProfiles must be set to 0!')
end  

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

Ib=find(izBox==1);
Ii=find(izBox~=1);
nbb=length(Ib);
nbi=length(Ii);

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

% Use T/S from GCM as boundary conditions.
numTracers=2;

load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Tsteady=gridToMatrix(Tgcm,[],boxFile,gridFile);
Ssteady=gridToMatrix(Sgcm,[],boxFile,gridFile);

clear Tgcm Sgcm % make some space

% Boundary conditions
Cbc{1}=Tsteady(Ib,:);
Cbc{2}=Ssteady(Ib,:);

% Initial condition
Cini{1}=Tsteady(Ii,1);
Cini{2}=Ssteady(Ii,1);

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
    writePetscBin(['Cini_' sprintf('%02d',itr) '.petsc'],Cini{itr})
  end
% Boundary conditions
  for itr=1:numTracers
	for im=1:12
	  writePetscBin(['Cbc_' sprintf('%02d',itr) '_' sprintf('%02d',im-1)],Cbc{itr}(:,im))
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
