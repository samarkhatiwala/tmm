% Set toplevel path to GCMs configuration
% base_path='/Volumes/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
base_path='/Volumes/data2/spk/TransportMatrixConfigs/UVicOSUpicdefault';

forcingType=2; % 0 (constant), 1 (periodic), 2 (time dependent)
numForcingFields=12*4
%
matrixType=1; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numMatrices=12*4

Ttd0=1765 % Time origin for time dependent matrices, forcing etc 

dt=28800; % time step to use
daysPerYear=365;

rearrangeProfiles=0
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
freshWaterForcingFile=fullfile(gcmDataPath,'FreshWaterForcing_gcm');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','dz','deltaT','gridType')

load(boxFile,'izBox','nb')

if rearrangeProfiles
  load(profilesFile,'Ip_pre','Ir_pre','Ip_post','Ir_post','Irr')
  Ip=Ip_pre;
  Ir=Ir_pre;
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

% Calculate forcing
dzb_surf=gridToMatrix(dz,Ib,boxFile,gridFile);
% make flux field here or load it from file
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
load(freshWaterForcingFile,'EmPgcm')
Sgcm=Sgcm(:,:,1,1:12);
EmPgcm=EmPgcm(:,:,1:12);
flux=1e-3*Sgcm.*EmPgcm;
fluxb=gridToMatrix(flux,Ib,boxFile,gridFile,1);
% enforce zero net flux
% dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);
% fluxb = fluxb - mean(dab_surf'*fluxb)/sum(dab_surf);
% convert flux to tendency applied in surface layer
qb=zeros([nb 12]);
for im=1:12
  qb(Ib,im)=fluxb(:,im)./dzb_surf;
end  
qb=qb*dt; % discretize in time

clear Sgcm EmPgcm flux fluxb % make some space

switch forcingType
  case 0
    q=mean(qb,2);
  case 1
    q=qb(:,1:12);
  case 2
    q=zeros([nb 12 numForcingFields/12]);
    q(:,:,1)=qb(:,1:12);
    for it=2:(numForcingFields/12)
      q(:,:,it)=q(:,:,1)*(1+0.05*it);
    end
    q=reshape(q,[nb numForcingFields]);
  otherwise
    error('ERROR!: Unknown forcingType')
end

% Initial condition
Cini=repmat(0,[nb 1]);

if rearrangeProfiles
  q=q(Ir,:);
  Cini=Cini(Ir,:);
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
    write_transport_matrices(base_path,dt,rearrangeProfiles,bigMat,useCoarseGrainedMatrix,matrixType,numMatrices)
  end

% Initial conditions
  writePetscBinVec('Cini.petsc',Cini)

% Forcing data
  switch forcingType
    case 0
      writePetscBinVec('q.petsc',q)
    case 1
      for im=1:12
        writePetscBinVec(['q_' sprintf('%02d',im-1)],q(:,im))
      end
    case 2	  
      % use the first numForcingFields time slices of each field and bracket them with the first and last fields
      tmp=q(:,[1 1:numForcingFields numForcingFields]);
      writePetscBinVec('qtd.petsc',tmp)
    otherwise
      error('ERROR!: Unknown forcingType')
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

