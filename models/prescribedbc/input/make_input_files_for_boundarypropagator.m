% Set toplevel path to GCMs configuration
base_path='/Volumes/data2/spk/TransportMatrixConfigs/UVicOSUpicdefault';
% base_path='/disks/backup/spk/Backup/Yellowstone/glade/p/cosc0001/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault_transient_ssp2-4.5';

% First run make_region_maps.m
regNums=[1:7];
regionFile='global_plus_regions';

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

load(regionFile,'numRegions','regionb')

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
% Simulate a subset of regions
regionb=regionb(:,regNums);
numRegions=size(regionb,2);

Gbc=zeros([nbb numRegions]);
for ir=1:numRegions 
  kr=find(regionb(:,ir));
  Gbc(kr,ir)=1;
end  

% Initial condition
Gini=zeros([nbi 1]);

if writeFiles
  calc_periodic_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],['periodic_times_' num2str(daysPerYear) 'd.bin']);
  if matrixType==2
    calc_time_dependent_times_for_tmm(['monthly-' num2str(daysPerYear) '-day year'],numMatrices/12,Ttd0,['matrix_times_' num2str(daysPerYear) 'd_' num2str(numMatrices) 'months.bin']);
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
  for ir=1:numRegions
    writePetscBinVec(['Gini_' sprintf('%04d',ir) '.petsc'],Gini)
  end

% Boundary conditions
  for ir=1:numRegions
	   writePetscBinVec(['Gbc_' sprintf('%04d',ir) '.petsc'],Gbc(:,ir))
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
