% Set toplevel path to GCMs configuration
base_path='/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault';

forcingType=0; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numForcingFields=12*335
%
matrixType=1; % 0 (annual mean), 1 (periodic), 2 (time dependent)
numMatrices=12*335

Ttd0=1765 % Time origin for time dependent matrices, forcing etc 

useEmP=0

dt=28800; % time step to use
daysPerYear=365;

rearrangeProfiles=1
bigMat=0
writeFiles=1
writeTMs=0
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

% model specific data
mobiInputDataFile=fullfile(bgcDataPath,'MOBI_input_data');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dz','dznom','da','x','y','z','deltaT','gridType')

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
nbb=length(Ib);

if rearrangeProfiles || bigMat
  load(profilesFile,'Ip_pre','Ir_pre','Ip_post','Ir_post','Irr')
  Ip=Ip_pre;
  Ir=Ir_pre;
end

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);

volFracSurf=zeros(nb,1);
volFracSurf(Ib)=volb(Ib)/sum(volb(Ib));  

% Surface forcing
if useEmP
  load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
  Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);
% Compute total E-P for virtual flux:
% Fv = TRg*(E-P), E-P in m/s
% d[TR]/dt = ... + Fv/dz
  load(freshWaterForcingFile,'EmPgcm','Srelaxgcm','saltRelaxTimegcm')
  zeroNetEmP=1;
  EmP=get_surface_emp_for_virtual_flux(gridFile,boxFile,EmPgcm,Srelaxgcm,saltRelaxTimegcm,Salt,zeroNetEmP);
%   gV=EmP./dzb;

  if fixEmP
	load(empFixFile,'empFixX','empFixY')
	nEmPFix=length(empFixX);
%   Points to fix: each SET is individually fixed.
	for k=1:nEmPFix
	  if length(empFixX{k})>1 % polygon  
		Ifix{k}=find(inpolygon(Xboxnom(Ib),Yboxnom(Ib),empFixX{k},empFixY{k})); % indexed to Ib
	  else % single point
		Ifix{k}=find(Xboxnom(Ib)==empFixX{k} & Yboxnom(Ib)==empFixY{k}); % referenced to Ib
	  end
	end
%   Rest of ocean (not fixed)
    Infix=find(~ismember([1:nbb]',cat(1,Ifix{:}))); % Everything else indexed to Ib	
%   Now fix E-P so that annual and spatial integral is 0.
	for k=1:length(Ifix) % loop over each SET of problematic points
	  if useAreaWeighting
		areabfix=dab_surf(Ifix{k});
		EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(areabfix'*EmP(Ifix{k},:))/sum(areabfix);
	  else
		EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(mean(EmP(Ifix{k},:)));
	  end
	end
	if useAreaWeighting
	  areabnfix=dab_surf(Infix);
	  EmP(Infix,:) = EmP(Infix,:) - mean(areabnfix'*EmP(Infix,:))/sum(areabnfix);
	else
	  EmP(Infix,:) = EmP(Infix,:) - mean(mean(EmP(Infix,:)));
	end
  end  
end

% Grid variables
dzb=gridToMatrix(dz,[],boxFile,gridFile);

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if forcingType==0
  if useEmP
    EmP=mean(EmP,2);
  end  
end

if rescaleForcing
  if forcingType==0
	Rfs=mean(Rfs,2);    
  end  
end

latb=Yboxnom(Ib);

% Initial condition
trNames=readtable('MOBI_tracer_names.txt','ReadVariableNames',0);
trNames=table2cell(trNames);

numTracers=length(trNames);
TRini=zeros(nb,numTracers);
for itr=1:numTracers
  fn=fullfile('InitialConditionProfiles',[trNames{itr} '.dat']);
  [hdr,dat]=hdrload(fn);
  TRini(:,itr)=interp1(dat(:,1),dat(:,2),Zboxnom);
  kk=find(Zboxnom<dat(1,1));
  TRini(kk,itr)=dat(1,2);
  kk=find(Zboxnom>dat(end,1));
  TRini(kk,itr)=dat(end,2);  
  if any(isnan(TRini(:,itr)))
    error('ERROR interpolating initial conditions!')
  end
end

if rearrangeProfiles
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  TRini=TRini(Ir,:);
  volFracSurf=volFracSurf(Ir);
  dzb=dzb(Ir);
  if rescaleForcing
    Rfs=Rfs(Ir,:);
  end    
  Ib=find(izBox==1);
%
  Ip=Ip_post;
  Ir=Ir_post;
end  

if useCoarseGrainedMatrix
% Coarse grain initial conditions and forcing data
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
  for itr=1:numTracers
    fn=[trNames{itr} 'ini.petsc'];
	writePetscBinVec(fn,TRini(:,itr))  
  end

% Forcing data
%   Salt=(Salt-35)/1000; % change to MOBI units
  if useEmP
    EmP=EmP*100; % change to MOBI units (cm/s)
	switch forcingType
	  case 0
		write_binary('EmP.bin',EmP,'real*8')
	  case 1
		for im=1:12
		  write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
		end    
	  case 2	  
%       use the first numForcingFields time slices of each field and bracket them with the first and last fields
		tmp=EmP(:,[1 1:numForcingFields numForcingFields]);
		write_binary('EmPtd.bin',tmp,'real*8')
	  otherwise
		error('ERROR!: Unknown forcing type')    
	end    
  end
  
  if rescaleForcing
    switch matrixType
	  case 0
		writePetscBinVec('Rfs.petsc',Rfs)
	  case 1
		for im=1:numMatrices
		  writePetscBinVec(['Rfs_' sprintf('%02d',im-1)],Rfs(:,im))
		end    
	  case 2
%       use the first numForcingFields time slices of each field and bracket them with the first and last fields
		tmp=Rfs(:,[1 1:numMatrices numMatrices]);
		writePetscBinVec('Rfs.petsc',tmp)
	  otherwise	
	    error('ERROR!: Unknown forcing type')    
	end    
  end  
  
% Grid data
  z=z*100; % convert to cm
  dznom=dznom*100; % convert to cm
  dzb=dzb*100; % convert to cm
  dab_surf=dab_surf*1e4; % convert to cm^2
  write_binary('zt.bin',nz,'int')  
  write_binary('zt.bin',z,'real*8',1)
  write_binary('drF.bin',nz,'int')  
  write_binary('drF.bin',dznom,'real*8',1)
  writePetscBinVec('dz.petsc',dzb)
  write_binary('latitude.bin',latb,'real*8')
  write_binary('dA.bin',dab_surf,'real*8')    
  writePetscBinVec('surface_volume_fraction.petsc',volFracSurf)  
  
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
  if useCoarseGrainedMatrix
    A=Beta*A*M;
    save pc_cg_data A nbbcg CGgrid CG Ipcg Ibcg dt
  else
    save pc_data A nbb nz nb Ip Ib dt
  end
end
