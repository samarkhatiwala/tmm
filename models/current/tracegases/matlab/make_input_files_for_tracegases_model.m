% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

% Specify the trace gas here. Valid options are: 'N2O', 'CFC11', 'CFC12' and 'SF6'.
gasID='CFC11';

periodicForcing=1
timeDependentForcing=0
numForcingMats=4020
Ttd0=1765
periodicMatrix=1
timeDependentMatrix=0

dt=43200; % time step to use

rearrangeProfiles=1
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0
writePCFiles=0

useAtmModel=0
useEmissions=0
%
useTimeVaryingPrescribedMixingRatio=1
useVirtualFlux=0
empScaleFactor=0.0 % OCMIP-2 protocol does not require E-P contribution

oceanCarbonBasePath='/data2/spk/OceanCarbon';
%-----------------------------------
% Options for wind/gas exchange
% (1) Gas exchange velocity using OCMIP protocol (if periodicForcing==0) or 
%     exchange coefficient using OCMIP protocol (if periodicForcing==1)
% (2) Gas exchange computed online from winds; default is to read CORE-2 winds
gasExchangeType=2;
corePath=fullfile(oceanCarbonBasePath,'CORE');
%-----------------------------------

% Atmospheric mixing ratio is either prescribed or computed prognostically 
% by setting useAtmModel=1 above. In the latter case it is also possible to 
% force the atmosphere with prescribed emissions by setting useEmissions=1

% data for atmospheric mixing ratio
atmosDataPath=fullfile(oceanCarbonBasePath,'MiscData');

% data for atm model
xTRatm_ini=0
% 
emissionsFile='';

if useEmissions & ~useAtmModel
  error('Cannot prescribe emissions without turning on atmospheric model!')
end  

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

iceFile=fullfile(bgcDataPath,'ice_fraction');
uwindFile=fullfile(bgcDataPath,'uwind');
vwindFile=fullfile(bgcDataPath,'vwind');
windFile=fullfile(bgcDataPath,'wind_speed');
atmospFile=fullfile(bgcDataPath,'atmospress');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','dz','da','x','y','z','deltaT','gridType')

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
nbb=length(Ib);

if ~useAtmModel
  if useTimeVaryingPrescribedMixingRatio
%   Load atmospheric mixing ratio history. 
    tab=readtable(fullfile(atmosDataPath,'CFC_atmospheric_histories_revised_2015_Table1.csv'),'HeaderLines',6);
    xTR_atm=table2array(tab);
	Tatm=xTR_atm(:,1);	
    switch gasID
      case 'N2O'
        ig=[12:13];
        useSpatiallyVariablePrescribedMixingRatio=1;
        latS=-10;
        latN=10;
      case {'CFC11','CFC-11'}
        ig=[2:3];
        useSpatiallyVariablePrescribedMixingRatio=1;
        latS=-10;
        latN=10;
      case {'CFC12','CFC-12'}
        ig=[4:5];
        useSpatiallyVariablePrescribedMixingRatio=1;
        latS=-10;
        latN=10;
      case {'SF6'}
        ig=[10:11];
        useSpatiallyVariablePrescribedMixingRatio=1;
        latS=-10;
        latN=10;
      otherwise
        error('Unknown gas!')  
     end
     if useSpatiallyVariablePrescribedMixingRatio
       xTR_atm=xTR_atm(:,ig);     
	   Ybb=Yboxnom(Ib);
	   jeq=find(Ybb>=latS & Ybb<=latN);
	   jnh=find(Ybb>latN);
	   jsh=find(Ybb<latS);
	   xTRatmb=repmat(0,[nbb length(Tatm)]);
	   xTRN=xTR_atm(:,1)'; % 1 x nt
	   xTRS=xTR_atm(:,2)'; % 1 x nt
%      interpolate between latS and latN	   
	   xTRatmb(jeq,:)=interp1([latS;latN],[xTRS;xTRN],Ybb(jeq),'linear');
%      NH
       xTRatmb(jnh,:)=repmat(xTRN,[length(jnh) 1]);
%      SH
       xTRatmb(jsh,:)=repmat(xTRS,[length(jsh) 1]);
     else
	   xTRatmb=xTR_atm(:,ig);
	 end
  end  
else 
  if useEmissions
    [hdr,emis]=hdrload(fullfile(atmosDataPath,emissionsFile));
    Tem=emis(:,1); % time
    TRem=emis(:,2); % annual emissions [mass/year]
  end    
end

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

clear Tgcm Sgcm % make some space

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);
% surface layer thickness
dzb_surf=gridToMatrix(dz,Ib,boxFile,gridFile);

% Surface forcing
% Compute total E-P for virtual flux:
% Fv = TRg*(E-P), E-P in m/s
% d[TR]/dt = ... + Fv/dz
load(freshWaterForcingFile,'EmPgcm','Srelaxgcm','saltRelaxTimegcm')
zeroNetEmP=1;
EmP=get_surface_emp_for_virtual_flux(gridFile,boxFile,EmPgcm,Srelaxgcm,saltRelaxTimegcm,Salt,zeroNetEmP);
if fixEmP
  load(empFixFile,'empFixX','empFixY')
  nEmPFix=length(empFixX);
% Points to fix: each SET is individually fixed.
  for k=1:nEmPFix
	if length(empFixX{k})>1 % polygon  
	  Ifix{k}=find(inpolygon(Xboxnom(Ib),Yboxnom(Ib),empFixX{k},empFixY{k})); % indexed to Ib
	else % single point
	  Ifix{k}=find(Xboxnom(Ib)==empFixX{k} & Yboxnom(Ib)==empFixY{k}); % referenced to Ib
	end
  end
% Rest of ocean (not fixed)
  Infix=find(~ismember([1:nbb]',cat(1,Ifix{:}))); % Everything else indexed to Ib	
% Now fix E-P so that annual and spatial integral is 0.
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
EmP=EmP*empScaleFactor;

if useVirtualFlux
  volFracSurf=zeros(nb,1);
  volFracSurf(Ib)=volb(Ib)/sum(volb(Ib));  
end

% surface T and S (note: we do this here as the full 3-d salinity is needed above)
Thetas=Theta(Ib,:);
Salts=Salt(Ib,:);

clear Theta Salt % make some space

switch gasExchangeType
  case (1)
%   Compute gas exchange velocity/exchange coefficient using OCMIP data/protocol
	if periodicForcing
	  xkwb=load_ocmip_variable([],'XKW',Xboxnom(Ib),Yboxnom(Ib));
	  ficeb=load_ocmip_variable([],'FICE',Xboxnom(Ib),Yboxnom(Ib));
	  ficeb(ficeb<0.2)=0; % OCMIP-2 howto 	
	else 
	  Vgas660=calc_ocmip_piston_velocity([],Xboxnom(Ib),Yboxnom(Ib),Thetas);
	end
	atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
  case (2)
%   Compute gas exchange velocity online from CORE winds and OCMIP fields
	[u10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'u_10.15JUNE2009.nc'),'U_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
	[v10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'v_10.15JUNE2009.nc'),'V_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
	ficeb=load_ocmip_variable([],'FICE',Xboxnom(Ib),Yboxnom(Ib));
	ficeb(ficeb<0.2)=0; % OCMIP-2 howto  
	atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));    
  case (3)
%   Compute gas exchange velocity online from model winds and fields
	load(iceFile,'Fice')
	load(atmospFile,'atmospress')  
	if ~isempty(uwindFile)
	  load(uwindFile,'uwind')
	  load(vwindFile,'vwind')
	else  
	  load(windFile,'windspeed')
%     fake it
	  uwind=windspeed;
	  vwind=zeros(size(windspeed));
	end  
	u10b=gridToMatrix(uwind,Ib,boxFile,gridFile,1);
	v10b=gridToMatrix(vwind,Ib,boxFile,gridFile,1);  
	ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
	atmospb=gridToMatrix(atmospress,Ib,boxFile,gridFile,1);
  case (4)
%   Compute gas exchange velocity online from CORE winds and model fields
	load(iceFile,'Fice')
	[u10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'u_10.15JUNE2009.nc'),'U_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
	[v10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'v_10.15JUNE2009.nc'),'V_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
	ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
	ficeb(ficeb<0.2)=0; % OCMIP-2 howto 
	atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));    
  otherwise
	error('ERROR: Unknown option for gas exchange')
end

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if ~periodicForcing
  Thetas=mean(Thetas,2);
  Salts=mean(Salts,2);
  EmP=mean(EmP,2);
  if gasExchangeType==1
	Vgas660=mean(Vgas660,2);
  end
end

if ~periodicMatrix
  if rescaleForcing
	Rfs=mean(Rfs,2);    
  end  
end

if ~periodicForcing
  atmospb=mean(atmospb,2);
end

% Initial condition
TR=repmat(0,[nb 1]);

if rearrangeProfiles
  TR=TR(Ir); % initial condition
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  volb=volb(Ir);
  if useVirtualFlux
    volFracSurf=volFracSurf(Ir);
  end  
  if rescaleForcing
    Rfs=Rfs(Ir,:);
  end  
  Ib=find(izBox==1);
%
  Ip=Ip_post;
  Ir=Ir_post;
end  

if useCoarseGrainedMatrix
% Coarse grain initial conditions
end

if writeFiles
  calc_periodic_times_for_tmm('monthly-365-day year','periodic_times_365d.bin');
  calc_periodic_times_for_tmm('monthly-360-day year','periodic_times_360d.bin');  
  calc_periodic_times_for_tmm('6-hourly','periodic_times_6hourly.bin');
  if timeDependentForcing==1 || timeDependentMatrix==1
    bpp=[31 28 31 30 31 30 31 31 30 31 30 31]; % base year
    N0=repmat(bpp,[1 (numForcingMats/nm)]);
    T0=calc_periodic_times_for_tmm(N0)*(numForcingMats/nm);
    T=zeros([length(T0)+2 1]); % extend times on either side
    T(2:end-1)=T0;
    T(1)=-bpp(1)/2/sum(bpp);
    T(end)=T0(end)+bpp(end)/365;
    T=T+Ttd0;
    write_binary(['times_365d_' num2str(numForcingMats) 'months.bin'],T,'real*8')
  end
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
      if timeDependentMatrix==1
        ntmax=numForcingMats;
      else
        ntmax=nm;
      end    
	  for im=1:ntmax 
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
		if timeDependentMatrix==1
		  it=im;
		else
		  it=im-1;
		end
		writePetscBin(['Ae_' sprintf('%02d',it)],Aexp,[],1)
		if timeDependentMatrix==1
		  if im==1
		    it=im-1;
		  elseif im==ntmax
		    it=im+1;
		  end
		  writePetscBin(['Ae_' sprintf('%02d',it)],Aexp,[],1)		    
		end
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
      if timeDependentMatrix==1
        ntmax=numForcingMats;
      else
        ntmax=nm;
      end    
	  for im=1:ntmax 
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
		if timeDependentMatrix==1
		  it=im;
		else
		  it=im-1;
		end		
		writePetscBin(['Ai_' sprintf('%02d',it)],Aimp,[],1)
		if timeDependentMatrix==1
		  if im==1
		    it=im-1;
		  elseif im==ntmax
		    it=im+1;
		  end
		  writePetscBin(['Ai_' sprintf('%02d',it)],Aimp,[],1)		    
		end		
		clear Aimp
	  end
	end
  end	  	  
% Initial conditions  
  writePetscBin('trini.petsc',TR)
% Biogeochemical parameters
  write_binary('dzsurf.bin',dzb_surf,'real*8')    
  if ~useAtmModel 
    if useTimeVaryingPrescribedMixingRatio
	  write_binary('TxTR.bin',length(Tatm),'int')  
	  write_binary('TxTR.bin',Tatm,'real*8',1)    
      write_binary('xTRatm.bin',xTRatmb,'real*8')
    end
  else  
    write_binary('xTRatm_ini.bin',xTRatm_ini,'real*8')
    write_binary('dA.bin',dab_surf,'real*8')    
    if useEmissions
	  write_binary('Tem.bin',length(Tem),'int')  
	  write_binary('Tem.bin',Tem,'real*8',1)
      write_binary('emissions.bin',TRem,'real*8')          
    end
  end
  
% Surface forcing data
  if ~periodicForcing
	write_binary('Ts.bin',Thetas,'real*8')
	write_binary('Ss.bin',Salts,'real*8')  
    if gasExchangeType==1
      write_binary('Vgas660.bin',Vgas660,'real*8')
    end
	write_binary('atmosp.bin',atmospb,'real*8')
	write_binary('EmP.bin',EmP,'real*8')
	if rescaleForcing
	  writePetscBin('Rfs.petsc',Rfs)
    end	
  else
    if timeDependentForcing==1
%     use the first numForcingMats time slices of each field
      tmp=Thetas(:,[1 1:numForcingMats numForcingMats]);
	  write_binary('Ts.bin',tmp,'real*8')
      tmp=Salts(:,[1 1:numForcingMats numForcingMats]);
	  write_binary('Ss.bin',tmp,'real*8')
      if gasExchangeType==2 % OCMIP-2 fice; make copies here to fake time dependence
        disp('Warming: OCMIP-2 fice is being used; making copies to fake time dependence!')
        ficeb=repmat(ficeb,[1 numForcingMats/nm]);
      end  
      tmp=ficeb(:,[1 1:numForcingMats numForcingMats]);
	  write_binary('fice.bin',tmp,'real*8')
      tmp=EmP(:,[1 1:numForcingMats numForcingMats]);
	  write_binary('EmP.bin',tmp,'real*8')
	  if rescaleForcing
		tmp=Rfs(:,[1 1:numForcingMats numForcingMats]);
		writePetscBin('Rfs.petsc',tmp,1)
	  end
    else
	  for im=1:nm
		write_binary(['Ts_' sprintf('%02d',im-1)],Thetas(:,im),'real*8')
		write_binary(['Ss_' sprintf('%02d',im-1)],Salts(:,im),'real*8')	  
		write_binary(['fice_' sprintf('%02d',im-1)],ficeb(:,im),'real*8')
		write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
		if rescaleForcing
		  writePetscBin(['Rfs_' sprintf('%02d',im-1)],Rfs(:,im))		
        end				
      end    
    end    

    for im=1:nm
	  if gasExchangeType==1
        write_binary(['xkw_' sprintf('%02d',im-1)],xkwb(:,im),'real*8')
      end
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')
	end
	if gasExchangeType==2 || gasExchangeType==3 || gasExchangeType==4
	  for im=1:size(u10b,2)
        write_binary(['uwind_' sprintf('%02d',im-1)],u10b(:,im),'real*8')
        write_binary(['vwind_' sprintf('%02d',im-1)],v10b(:,im),'real*8')
      end
    end	  
  end
  if useVirtualFlux
    writePetscBin('surface_volume_fraction.petsc',volFracSurf)
  end  
  
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
