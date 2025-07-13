% Set toplevel path to GCMs configuration
base_path='/Volumes/data2/spk/TransportMatrixConfigs/UVicOSUpicdefault';
% base_path='/disks/backup/spk/Backup/Yellowstone/glade/p/cosc0001/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault_transient_ssp2-4.5';

bcType=1; % 0 (constant), 1 (periodic), 2 (time dependent)
numBCFields=12*4
%
matrixType=1; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numMatrices=12*4

Ttd0=1765 % Time origin for time dependent matrices, forcing etc 

dt=28800; % time step to use
daysPerYear=365;

rearrangeProfiles=0 % DO NOT CHANGE!
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

% Boundary conditions
% Use T/S from GCM as boundary conditions
numTracers=2;
load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);

clear Tgcm Sgcm % make some space

switch bcType
  case 0
    Cbc{1}=mean(Theta(Ib,1:12),2);
    Cbc{2}=mean(Salt(Ib,1:12),2);    
  case 1
    Cbc{1}=Theta(Ib,1:12);
    Cbc{2}=Salt(Ib,1:12);    
  case 2
    Cbc{1}=zeros([nbb 12 numBCFields/12]);
    Cbc{1}(:,:,1)=Theta(Ib,1:12);
    Cbc{2}=zeros([nbb 12 numBCFields/12]);
    Cbc{2}(:,:,1)=Salt(Ib,1:12);    
    for it=2:(numBCFields/12)
      Cbc{1}(:,:,it)=Cbc{1}(:,:,1)*(1+0.05*it);
      Cbc{2}(:,:,it)=Cbc{2}(:,:,1)*(1-0.05*it);
    end
    Cbc{1}=reshape(Cbc{1},[nbb numBCFields]);
    Cbc{2}=reshape(Cbc{2},[nbb numBCFields]);    
  otherwise
    error('ERROR!: Unknown bcType')
end

% Initial condition
Cini{1}=zeros(nbi,1);
Cini{2}=zeros(nbi,1);

if writeFiles
  calc_periodic_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],['periodic_times_' num2str(daysPerYear) 'd.bin']);
  if matrixType==2
    calc_time_dependent_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],numMatrices/12,Ttd0,['matrix_times_' num2str(daysPerYear) 'd_' num2str(numMatrices) 'months.bin']);
  end
  if bcType==2
    calc_time_dependent_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],numBCFields/12,Ttd0,['bc_times_' num2str(daysPerYear) 'd_' num2str(numBCFields) 'months.bin']);
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
    writePetscBinVec(['Cini_' sprintf('%02d',itr) '.petsc'],Cini{itr})
  end

% Boundary conditions
  switch bcType
    case 0
	     for itr=1:numTracers
        writePetscBinVec(['Cbc_' sprintf('%02d',itr) '.petsc'],Cbc{itr})
      end
    case 1
	     for itr=1:numTracers    
		      for im=1:12
		        writePetscBinVec(['Cbc_' sprintf('%02d',itr) '_' sprintf('%02d',im-1)],Cbc{itr}(:,im))
		      end
	     end	
    case 2	  
      % use the first numBCFields time slices of each field and bracket them with the first and last fields
   	  for itr=1:numTracers
		      tmp=Cbc{itr}(:,[1 1:numBCFields numBCFields]);
      		writePetscBinVec(['Cbctd_' sprintf('%02d',itr) '.petsc'],tmp)
   	  end	
    otherwise
      error('ERROR!: Unknown bcType')
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
