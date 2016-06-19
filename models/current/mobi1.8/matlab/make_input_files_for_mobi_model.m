% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicForcing=1
periodicMatrix=1

dt=43200; % time step to use

rearrangeProfiles=1
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0
writePCFiles=0

READ_SWRAD=1
useEmP=1;

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

if periodicForcing
  nm=12;
else
  nm=1;
end

% Use steady state T/S from GCM. Note we always load seasonal data here.
load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);

% Surface forcing
if useEmP
% Compute total E-P for virtual flux:
% Fv = TRg*(E-P), E-P in m/s
% d[TR]/dt = ... + Fv/dz
  load(freshWaterForcingFile,'EmPgcm','Srelaxgcm','saltRelaxTimegcm')
  zeroNetEmP=1;
  EmP=get_surface_emp_for_virtual_flux(gridFile,boxFile,EmPgcm,Srelaxgcm,saltRelaxTimegcm,Salt,zeroNetEmP);
%   gV=EmP./dzb;
  volFracSurf=zeros(nb,1);
  volFracSurf(Ib)=volb(Ib)/sum(volb(Ib));  

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
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
  if useEmP
    EmP=mean(EmP,2);
  end  
end

if ~periodicMatrix
  if rescaleForcing
	Rfs=mean(Rfs,2);    
  end  
end

% Surface forcing data
load(mobiInputDataFile,'aice')
% acie should be a fraction
aiceb=gridToMatrix(aice,Ib,boxFile,gridFile,1);
if ~periodicForcing
  aiceb=mean(aiceb,2);
end

load(mobiInputDataFile,'hice')
hiceb=gridToMatrix(hice,Ib,boxFile,gridFile,1);
if ~periodicForcing
  hiceb=mean(hiceb,2);
end

load(mobiInputDataFile,'hsno')
hsnob=gridToMatrix(hsno,Ib,boxFile,gridFile,1);
if ~periodicForcing
  hsnob=mean(hsnob,2);
end

load(mobiInputDataFile,'wind')
windb=gridToMatrix(wind,Ib,boxFile,gridFile,1);
if ~periodicForcing
  windb=mean(windb,2);
end

atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
if ~periodicForcing
  atmospb=mean(atmospb,2);
end

load(mobiInputDataFile,'Fe')
Feb=gridToMatrix(Fe,[],boxFile,gridFile);
if ~periodicForcing
  Feb=mean(Feb,2);
end

load(mobiInputDataFile,'sgbathy')
sgbathyb=gridToMatrix(sgbathy,[],boxFile,gridFile);
Ibot=cellfun(@max,Ip); % indices of bottom cell
sgbathyb(Ibot)=1; % required for PO4 conservation

if READ_SWRAD
% Read SW radiation from file
  load(mobiInputDataFile,'swrad')
  swradb=gridToMatrix(swrad,Ib,boxFile,gridFile,1);
  if ~periodicForcing
	swradb=mean(swradb,2);
  end  
end

latb=Yboxnom(Ib);

% Initial condition
trNames={'dic','dic13','c14','o2','alk','po4','dop','phyt','zoop','detr','no3','don',...
'diaz','din15','don15','phytn15','zoopn15','detrn15','diazn15','doc13','phytc13',...  
'zoopc13','detrc13','diazc13'};  

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
  Theta=Theta(Ir,:);
  Salt=Salt(Ir,:);
  Feb=Feb(Ir,:);
  sgbathyb=sgbathyb(Ir);
  TRini=TRini(Ir,:);
  if useEmP
    volFracSurf=volFracSurf(Ir);
  end  
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
	  if useCoarseGrainedMatrix
		Aexpms=Beta*Aexpms*M; % coarse-grained explicit transport matrix
	  end        
	  writePetscBin('Ae.petsc',Aexpms,[],1)
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
		if useCoarseGrainedMatrix
%         Not sure if this is really kosher!		  
		  Aexp=Beta*Aexp*M; % coarse-grained explicit transport matrix
		end
		writePetscBin(['Ae_' sprintf('%02d',im-1)],Aexp,[],1)
		clear Aexp
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
	  if useCoarseGrainedMatrix
		Aimpms=Beta*Aimpms*M; % coarse-grained implicit transport matrix
	  end
	  writePetscBin('Ai.petsc',Aimpms,[],1)
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
		if useCoarseGrainedMatrix
		  Aimp=Beta*Aimp*M; % coarse-grained implicit transport matrix		
		end
		writePetscBin(['Ai_' sprintf('%02d',im-1)],Aimp,[],1)
		clear Aimp
	  end
	end
  end	  	  
% Initial conditions  
  for itr=1:numTracers
    fn=[trNames{itr} 'ini.petsc'];
	writePetscBin(fn,TRini(:,itr))  
  end
  
  if ~periodicForcing
	write_binary('aice.bin',aiceb,'real*8')
	write_binary('hice.bin',hiceb,'real*8')
	write_binary('hsno.bin',hsnob,'real*8')	
	write_binary('wind.bin',windb,'real*8')	
	write_binary('atmosp.bin',atmospb,'real*8')	
  else
    for im=1:nm
	  write_binary(['aice_' sprintf('%02d',im-1)],aiceb(:,im),'real*8')
	  write_binary(['hice_' sprintf('%02d',im-1)],hiceb(:,im),'real*8')
	  write_binary(['hsno_' sprintf('%02d',im-1)],hsnob(:,im),'real*8')
	  write_binary(['wind_' sprintf('%02d',im-1)],windb(:,im),'real*8')	  
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')	  
	end
  end
  if READ_SWRAD
	if ~periodicForcing  
	  write_binary('swrad.bin',swradb,'real*8')
	else
	  for im=1:nm	
		write_binary(['swrad_' sprintf('%02d',im-1)],swradb(:,im),'real*8')
	  end  
    end	
  end
  Salt=(Salt-35)/1000; % change to MOBI units
  if ~periodicForcing
	writePetscBin('Ts.petsc',Theta)
	writePetscBin('Ss.petsc',Salt)
  else
    for im=1:nm
	  writePetscBin(['Ts_' sprintf('%02d',im-1)],Theta(:,im))
	  writePetscBin(['Ss_' sprintf('%02d',im-1)],Salt(:,im))
    end    
  end    
  Feb=Feb*1e9; % change to MOBI units
  if ~periodicForcing
	writePetscBin('Fe.petsc',Feb)
  else
    for im=1:nm
	  writePetscBin(['Fe_' sprintf('%02d',im-1)],Feb(:,im))
    end    
  end      

  writePetscBin('sgbathy.petsc',sgbathyb)
  
  if useEmP
    EmP=EmP*100; % change to MOBI units (cm/s)
    writePetscBin('surface_volume_fraction.petsc',volFracSurf)
    if ~periodicForcing
  	  write_binary('EmP.bin',EmP,'real*8')
    else
	  for im=1:nm
  	    write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
	  end
	end    
  end
  
  if rescaleForcing
	if ~periodicMatrix
	  writePetscBin('Rfs.petsc',Rfs)
	else
	  for im=1:nm
		writePetscBin(['Rfs_' sprintf('%02d',im-1)],Rfs(:,im))
	  end    
	end    
  end  
  
% Grid data
  z=z*100; % convert to cm
  dznom=dznom*100; % convert to cm
  dzb=dzb*100; % convert to cm  
  write_binary('zt.bin',nz,'int')  
  write_binary('zt.bin',z,'real*8',1)
  write_binary('drF.bin',nz,'int')  
  write_binary('drF.bin',dznom,'real*8',1)
  writePetscBin('dz.petsc',dzb)
  write_binary('latitude.bin',latb,'real*8')    
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
