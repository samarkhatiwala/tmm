% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicMatrix=0

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

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

if strcmp(gridType,'llc_v4')
  load(boxFile,'YboxnomGlob','izBoxGlob','nb')
  Yboxnom=YboxnomGlob;  
  nb=sum(nb);
  izBox=izBoxGlob;
else
  load(boxFile,'Yboxnom','izBox','nb')
end

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

% Define surface patches
% First load basin mask file
load(fullfile(gcmDataPath,'basin_mask'),'basin_mask','basin_names')
basin_mask=gridToMatrix(basin_mask,Ib,boxFile,gridFile);

Ybb=Yboxnom(Ib);

numPatches=3;
patchb=zeros([nbb numPatches]);
% (1) North Atlantic north of 50 N
basinIndex=find(strcmp(basin_names,'Atlantic'));
ipatch=find(basin_mask==basinIndex & Ybb>50);
patchb(ipatch,1)=1;
% (2) Southern Ocean south of 55 S
ipatch=find(Ybb<=-55);
patchb(ipatch,2)=2;
% (3) Global surface
patchb(:,3)=3;

% Boundary conditions
Gbc=zeros([nbb numPatches]);
for ip=1:numPatches
  Gbc(find(patchb(:,ip)),ip)=1;
end

% Initial conditions
Gini=repmat(0,[nbi numPatches]);

if rearrangeProfiles
  error('ERROR: rearrangeProfiles must be set to 0!')
end  

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
  for ip=1:numPatches
    writePetscBinOrig(['Gini_' sprintf('%02d',ip) '.petsc'],Gini(:,ip))
  end
% Boundary conditions
  for ip=1:numPatches  
    writePetscBinOrig(['Gbc_' sprintf('%02d',ip) '.petsc'],Gbc(:,ip))
  end  
  
% Grid data

% Profile data
  if rearrangeProfiles
    if ~useCoarseGrainedMatrix
      gStartIndices=repmat(0,[nbb 1]);
      gEndIndices=repmat(0,[nbb 1]);
      for is=1:nbb % loop over each surface point
        Ipl=Ip{is}; % indices for local profile (globally indexed)  
        gStartIndices(is)=Ipl(1);
        gEndIndices(is)=Ipl(end);
      end
    else % useCoarseGrainedMatrix
      gStartIndices=repmat(0,[nbbcg 1]);
      gEndIndices=repmat(0,[nbbcg 1]);
      for is=1:nbbcg % loop over each surface point
        Ipl=Ipcg{is}; % indices for local profile (globally indexed)  
        gStartIndices(is)=Ipl(1);
        gEndIndices(is)=Ipl(end);
      end  
    end  
    write_binary('gStartIndices.bin',[length(gStartIndices);gStartIndices],'int')
    write_binary('gEndIndices.bin',[length(gEndIndices);gEndIndices],'int')
  end
end
