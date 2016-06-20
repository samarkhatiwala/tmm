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

useAtmModel=0
useLandModel=0
useEmissions=0
%
useTimeVaryingPrescribedCO2=1
doC14calc=1
% iniConditionFile='preindustrial_steady_state_dic'
useVirtualFlux=0
empScaleFactor=1.0

oceanCarbonBasePath='/data2/spk/OceanCarbon';
%-----------------------------------
% Options for wind/gas exchange
% (1) Gas exchange velocity using OCMIP protocol (if periodicForcing==0) or 
%     exchange coefficient using OCMIP protocol (if periodicForcing==1)
% (2) Gas exchange computed online from winds; default is to read CORE-2 winds
gasExchangeType=2;
corePath=fullfile(oceanCarbonBasePath,'CORE');
%-----------------------------------

% Atmospheric CO2 and, optionally, C14 are either prescribed or computed prognostically 
% by setting useAtmModel=1 above. In the latter case it is also possible to: 
% (1) use a land vegetation model by also setting useLandModel=1 above
% (2) force the atmosphere with prescribed emissions by setting useEmissions=1

% data for atmospheric CO2 and C14
% Available options: 'historical', 'RCP3PD', 'RCP45', 'RCP6' and 'RCP85'
co2Scenario='historical';
atmosDataPath=fullfile(oceanCarbonBasePath,'AtmosphericCarbonData');

% data for land-atm model
pCO2atm_ini=277.62
% Available emission options: 'historical', 'RCP3PD', 'RCP45', 'RCP6' and 'RCP85'
emissionScenario='RCP3PD';
atmosDataPath=fullfile(oceanCarbonBasePath,'AtmosphericCarbonData');

if useLandModel & ~useAtmModel
  error('Cannot use land model without turning on atmospheric model!')
end  
if useLandModel
  useLandAtmModel=1;
end  
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
  if useTimeVaryingPrescribedCO2
%   Load atmospheric pCO2 history. 
    if strcmp(co2Scenario,'historical')
      co2File='co2_HG11_extrap.dat'; % historical data from Heather Graven
    else
      co2File=[co2Scenario '_CO2_conc.txt']; % other scenarios
    end      
	[hdr,pCO2_atm]=hdrload(fullfile(atmosDataPath,co2File));
	Tco2=pCO2_atm(:,1);
	pCO2atm=pCO2_atm(:,2);
  end

  if doC14calc
    if useTimeVaryingPrescribedCO2
	  if strcmp(co2Scenario,'historical')    
		Ybb=Yboxnom(Ib);
		jeq=find(Ybb>=-20 & Ybb<=20);
		jnh=find(Ybb>20);
		jsh=find(Ybb<-20);
  %     Historical data from Heather Graven
		[hdr,DC14atm]=hdrload(fullfile(atmosDataPath,'c14equ_HG09.dat'));
		TC14atm=DC14atm(:,1);
		DC14atm=DC14atm(:,2);  
		DC14atmb=repmat(0,[nbb length(TC14atm)]);
  %     EQ
		DC14atmb(jeq,:)=repmat(DC14atm',[length(jeq) 1]);
  %     NH
		[hdr,DC14atm]=hdrload(fullfile(atmosDataPath,'c14nth_HG09.dat'));
		DC14atm=DC14atm(:,2);
		DC14atmb(jnh,:)=repmat(DC14atm',[length(jnh) 1]);
  %     SH
		[hdr,DC14atm]=hdrload(fullfile(atmosDataPath,'c14sth_HG09.dat'));
		DC14atm=DC14atm(:,2);
		DC14atmb(jsh,:)=repmat(DC14atm',[length(jsh) 1]);    
	  else
		error('ERROR: only historical scenario is currently available')
	  end      
    else
      DC14atmb=repmat(0,[nbb 1]); % for equilibrium run
    end
  end
  
else 
  if useEmissions
    emFile=[emissionScenario '_emissions.txt'];
    [hdr,emis]=hdrload(fullfile(atmosDataPath,emFile));
    Tem=emis(:,1); % time
    FFem=emis(:,2); % fossil fuel emissions
    LUem=emis(:,3); % land use emissions
  end  
  
  if useLandModel
%   calculate steady state land model   
    [Cvini,Csini,Cvsini]=calc_steadystate_land(pCO2atm_ini);
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

if gasExchangeType==1
% Compute gas exchange velocity/exchange coefficient using OCMIP data/protocol
  if periodicForcing
	xkwb=load_ocmip_variable([],'XKW',Xboxnom(Ib),Yboxnom(Ib));
 	ficeb=load_ocmip_variable([],'FICE',Xboxnom(Ib),Yboxnom(Ib));
    ficeb(ficeb<0.2)=0; % OCMIP-2 howto 	
  else 
    Vgas=calc_ocmip_piston_velocity([],Xboxnom(Ib),Yboxnom(Ib),Theta(Ib,:),'CO2');
  end
  atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
elseif gasExchangeType==2
% Compute gas exchange velocity online from winds
  [u10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'u_10.15JUNE2009.nc'),'U_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
  [v10b,Tcore,lon,lat]=load_core_variable(fullfile(corePath,'v_10.15JUNE2009.nc'),'V_10_MOD',Xboxnom(Ib),Yboxnom(Ib));
  ficeb=load_ocmip_variable([],'FICE',Xboxnom(Ib),Yboxnom(Ib));
  ficeb(ficeb<0.2)=0; % OCMIP-2 howto  
  atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));    
elseif gasExchangeType==3
% Compute gas exchange velocity online from model winds
  load(iceFile,'Fice')
  load(atmospFile,'atmospress')  
  if ~isempty(uwindFile)
	load(uwindFile,'uwind')
	load(vwindFile,'vwind')
  else  
	load(windFile,'windspeed')
%   fake it
    uwind=windspeed;
    vwind=zeros(size(windspeed));
  end  
  u10b=gridToMatrix(uwind,Ib,boxFile,gridFile,1);
  v10b=gridToMatrix(vwind,Ib,boxFile,gridFile,1);  
  ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
  atmospb=gridToMatrix(atmospress,Ib,boxFile,gridFile,1);
else
  error('ERROR: Unknown option for gas exchange')
end

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
  EmP=mean(EmP,2);
  if gasExchangeType==1
	Vgas=mean(Vgas,2);
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

% Biogeochem parameters
% Some constants for carbonate chemistry
rho0=1024.5; % kg/m^3 nominal density from OCMIP-2
SiO2avg=7.7e-3; % mol/m^3 
PO4avg=5.1e-4; % mol/m^3
DICavg=2230.0*rho0*1.0e-6; % mol/m^3
meanSurfaceSalinity=mean((volb(Ib)/sum(volb(Ib)))'*Salt(Ib,:)); % volume weighted, annual mean surface salinity
Sbar=meanSurfaceSalinity;
Alkbar=2310*1e-6*rho0; % ueq/kg -> eq/m^3  OCMIP-2
alkFunc = @(S) Alkbar*S/Sbar;

SiO2b=repmat(SiO2avg,[nbb nm]);
PO4b=repmat(PO4avg,[nbb nm]);
Alkb=alkFunc(Salt(Ib,:));

% Initial condition
if exist('iniConditionFile','var')
  disp(['LOADING initial condition from ' iniConditionFile])
  load(iniConditionFile,'C');
  DIC=C;
else
  DIC=repmat(DICavg,[nb 1]);
  if doC14calc
	DIC14=repmat(DICavg,[nb 1]);
  end  
end

if rearrangeProfiles
  DIC=DIC(Ir); % initial condition
  if doC14calc
	DIC14=DIC14(Ir); % initial condition
  end      
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  Theta=Theta(Ir,:);
  Salt=Salt(Ir,:);
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
  writePetscBin('dicini.petsc',DIC)
  if doC14calc
    writePetscBin('dic14ini.petsc',DIC14)
  end  
% Biogeochemical parameters
  write_binary('dzsurf.bin',dzb_surf,'real*8')    
  if ~useAtmModel 
    if useTimeVaryingPrescribedCO2
	  write_binary('TpCO2.bin',length(Tco2),'int')  
	  write_binary('TpCO2.bin',Tco2,'real*8',1)    
      write_binary('pCO2atm.bin',pCO2atm,'real*8')
    end
	if doC14calc
	  write_binary('DC14atm.bin',DC14atmb,'real*8')  
	  if useTimeVaryingPrescribedCO2
		write_binary('TC14atm.bin',length(TC14atm),'int')  
		write_binary('TC14atm.bin',TC14atm,'real*8',1)    
	  end
	end    
  else  
    write_binary('pCO2atm_ini.bin',pCO2atm_ini,'real*8')
    write_binary('dA.bin',dab_surf,'real*8')    
    if useEmissions
	  write_binary('Tem.bin',length(Tem),'int')  
	  write_binary('Tem.bin',Tem,'real*8',1)
      write_binary('fossil_fuel_emissions.bin',FFem,'real*8')          
      write_binary('land_use_emissions.bin',LUem,'real*8')          
    end
    if useLandModel
      write_binary('land_ini.bin',[Cvini Csini Cvsini]','real*8')    
    end
  end
  
% Surface forcing data
  if ~periodicForcing
    if gasExchangeType==1
      write_binary('Vgas.bin',Vgas,'real*8')
    end
	write_binary('atmosp.bin',atmospb,'real*8')
	write_binary('PO4.bin',PO4b,'real*8')
	write_binary('SiO2.bin',SiO2b,'real*8')
	write_binary('Alk.bin',Alkb,'real*8')
  else
    for im=1:nm
	  write_binary(['fice_' sprintf('%02d',im-1)],ficeb(:,im),'real*8')
	  if gasExchangeType==1
        write_binary(['xkw_' sprintf('%02d',im-1)],xkwb(:,im),'real*8')
      end
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')
	  write_binary(['PO4_' sprintf('%02d',im-1)],PO4b(:,im),'real*8')
	  write_binary(['SiO2_' sprintf('%02d',im-1)],SiO2b(:,im),'real*8')
	  write_binary(['Alk_' sprintf('%02d',im-1)],Alkb(:,im),'real*8')	  
	end
	if gasExchangeType==2 || gasExchangeType==3
	  for im=1:size(u10b,2)
        write_binary(['uwind_' sprintf('%02d',im-1)],u10b(:,im),'real*8')
        write_binary(['vwind_' sprintf('%02d',im-1)],v10b(:,im),'real*8')
      end
    end	  
  end
  if ~periodicForcing
	writePetscBin('Ts.petsc',Theta)
	writePetscBin('Ss.petsc',Salt)
  else
    for im=1:nm
	  writePetscBin(['Ts_' sprintf('%02d',im-1)],Theta(:,im))
	  writePetscBin(['Ss_' sprintf('%02d',im-1)],Salt(:,im))
    end    
  end    
  if ~periodicForcing
	write_binary('EmP.bin',EmP,'real*8')
  else
	for im=1:nm
	  write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
	end
  end      
  if useVirtualFlux
    writePetscBin('surface_volume_fraction.petsc',volFracSurf)
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
