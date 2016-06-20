% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/disks/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/disks/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';

regNums=[1:7];
regionFile='global_plus_regions';

load(regionFile,'numRegions','regionb')

% simulate a subset of regions
regionb=regionb(:,regNums);
numRegions=size(regionb,2);

periodicMatrix=0

dt=43200; % time step to use

bigMat=0
rearrangeProfiles=0 % DO NOT CHANGE!
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

load(boxFile,'nb','ixBox','iyBox','izBox','volb','Xboxnom','Yboxnom','Zboxnom')

Ib=find(izBox==1);
Ii=find(izBox~=1);

nbb=length(Ib);
nbi=length(Ii);

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(explicitAnnualMeanMatrixFile,'Aexpms')
load(implicitAnnualMeanMatrixFile,'Aimpms')
if exist('dtMultiple','var')
  disp(['WARNING: dtMultiple is set to ' num2str(dtMultiple)])
  disp('Modifying Aimpms')    
  if bigMat % big matrix. do it a block at a time.
    load(profilesFile,'Ip_pre')  
	for is=1:nbb
	  for is=1:nbb % change time step multiple
		Aimpms(Ip_pre{is},Ip_pre{is})=Aimpms(Ip_pre{is},Ip_pre{is})^dtMultiple;
	  end
	end
  else
    Aimpms=Aimpms^dtMultiple;
  end
end
% make discrete
I=speye(nb,nb);
Aexpms=dt*Aexpms;
Aexpms=I+Aexpms;

% Boundary conditions
Gbc=zeros([nbb numRegions]);
for ir=1:numRegions 
  kr=find(regionb(:,ir));
  Gbc(kr,ir)=1;
end  
% Initial condition
Gini=zeros([nbi 1]);

if writeFiles
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
  for ir=1:numRegions
    writePetscBin(['Gini_' sprintf('%04d',ir) '.petsc'],Gini)
  end
% Boundary conditions
  for ir=1:numRegions
	writePetscBin(['Gbc_' sprintf('%04d',ir) '.petsc'],Gbc(:,ir))
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
