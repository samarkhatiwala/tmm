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
doAnthroCalc=0
doC14calc=0
empScaleFactor=1.0
ALLOW_FE=1
LIGHT_CHL=0
READ_PAR=0

% Some biology params. Make sure these match values in dic_biotic_params.F
% timescale for biological activity
alpha_val=1.5e-3/(24*60*60*360);
% inorganic/organic carbon rain ratio
rain_ratio_val=0.07;

oceanCarbonBasePath='/data2/spk/OceanCarbon';
corePath=fullfile(oceanCarbonBasePath,'CORE');

% data for land-atm model
pCO2atm_ini=277.62
emissionScenario='A2'
emissionModel='ASF'
emissionPath=fullfile(oceanCarbonBasePath,'FrancescaEmissionData');

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
windFile=fullfile(bgcDataPath,'wind_speed');
silicaFile=fullfile(bgcDataPath,'silica');
if ALLOW_FE  
  ironFile=fullfile(bgcDataPath,'iron_flux');  
end  
if LIGHT_CHL
  chlFile=fullfile(bgcDataPath,'chlorophyll');
end

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
  if doAnthroCalc
% Load atmospheric pCO2 history. 
% Data from Francesca Terenzi
%   pCO2_atm=load(fullfile(oceanCarbonBasePath,'pCO2_atm_from_fra.dat'));
%   Tco2=pCO2_atm(:,1);
%   pCO2atm=pCO2_atm(:,2);
% Smoothed OCMIP data (from Sara Fletcher)
%   pCO2_atm=load(fullfile(oceanCarbonBasePath,'OCMIP2','splco2_cis92a.dat'));
%   Tco2=pCO2_atm(:,1);
%   pCO2atm=pCO2_atm(:,2);
% Construct my own history from OCMIP2 and more recent data
%     load splco2_mod.dat  
%     load co2_annmean_gl.txt
%     Tco2_dat=[splco2_mod(:,1);co2_annmean_gl(23:end,1)];   
%     co2_dat=[splco2_mod(:,2);co2_annmean_gl(23:end,2)];
%     Tf=[1765:.1:2008]';       
%     co2f=interp1(Tco2_dat,co2_dat,Tf);
%     co2s=co2f;
%     for ism=1:5
%       co2s=lapsmooth(co2s,0);
%     end
%     Tco2=[1765:.5:2008.5]';
%     pCO2atm=interp1(Tf,co2s,Tco2,'spline','extrap');
% Data from Heather Graven
  [hdr,pCO2_atm]=hdrload(fullfile(oceanCarbonBasePath,'HeatherGraven','co2_HG11_extrap.dat'));
  Tco2=pCO2_atm(:,1);
  pCO2atm=pCO2_atm(:,2);
% Units: pCO2 contains the CO2 mole fraction (mixing ratio) in dry air in ppm. 
% This is numerically the same as the partial pressure in uatm of CO2 in dry air
% at total pressure of 1 atm.
  end

  if doC14calc
    if doAnthroCalc
      Ybb=Yboxnom(Ib);
      jeq=find(Ybb>=-20 & Ybb<=20);
      jnh=find(Ybb>20);
      jsh=find(Ybb<-20);
%     Data from Heather Graven
      [hdr,DC14atm]=hdrload(fullfile(oceanCarbonBasePath,'HeatherGraven','c14equ_HG09.dat'));
      TC14atm=DC14atm(:,1);
      DC14atm=DC14atm(:,2);  
      DC14atmb=repmat(0,[nbb length(TC14atm)]);
%     EQ
      DC14atmb(jeq,:)=repmat(DC14atm',[length(jeq) 1]);
%     NH
      [hdr,DC14atm]=hdrload(fullfile(oceanCarbonBasePath,'HeatherGraven','c14nth_HG09.dat'));
      DC14atm=DC14atm(:,2);
      DC14atmb(jnh,:)=repmat(DC14atm',[length(jnh) 1]);
%     SH
      [hdr,DC14atm]=hdrload(fullfile(oceanCarbonBasePath,'HeatherGraven','c14sth_HG09.dat'));
      DC14atm=DC14atm(:,2);
      DC14atmb(jsh,:)=repmat(DC14atm',[length(jsh) 1]);    
    else
      DC14atmb=repmat(0,[nbb 1]); % for equilibrium run
    end
  end
  
else 
  if useEmissions
    emissionName=[emissionScenario '_' emissionModel];
    [T_E,Ev,T_D,Dv]=get_emissions_pCO2atm(emissionName,emissionPath);
  end  
  if useLandModel
%   calculate steady state land model here using   
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

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);

% Surface forcing
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
EmP=EmP*empScaleFactor;

% Grid variables
dzb=gridToMatrix(dz,[],boxFile,gridFile);
if strcmp(gridType,'llc_v4')
  nF=length(nx);
  tmp=cell(nF,1);
  for iF=1:nF
    tmp{iF}=repmat(0,[nx{iF} ny{iF} nz]);
	for k=1:nz
	  tmp{iF}(:,:,k)=dznom(k);
	end
  end
else  
  tmp=repmat(0,[nx ny nz]);
  for k=1:nz
	tmp(:,:,k)=dznom(k);
  end
end
dznomb=gridToMatrix(tmp,[],boxFile,gridFile);
hFacCb=dzb./dznomb;
recip_hFacCb=1./hFacCb;

% Biogeochem data
% timescale for biological activity
alpha=repmat(alpha_val,[nbb 1]);
% inorganic/organic carbon rain ratio
rain_ratio=repmat(rain_ratio_val,[nbb 1]);

load(iceFile,'Fice')
load(windFile,'windspeed')
load(silicaFile,'silica')
if ALLOW_FE
  load(ironFile,'ironflux')
end
if LIGHT_CHL
  load(chlFile,'chl')
end

Ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
windb=gridToMatrix(windspeed,Ib,boxFile,gridFile,1);
silicab=gridToMatrix(silica,Ib,boxFile,gridFile,1);
atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
if ALLOW_FE
  inputfeb=gridToMatrix(ironflux,Ib,boxFile,gridFile,1);
end
if LIGHT_CHL
  chlb=gridToMatrix(chl,Ib,boxFile,gridFile,1);
end

if READ_PAR
% specify PAR
  parb=repmat(100,[nbb 12]);
else
% specify latitude
  latb=Yboxnom(Ib);
end

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
  EmP=mean(EmP,2);
  Ficeb=mean(Ficeb,2);
  windb=mean(windb,2);  
  silicab=mean(silicab,2);    
  atmospb=mean(atmospb,2);
  if READ_PAR
    parb=mean(parb,2);
  end  
  if ALLOW_FE
    inputfeb=mean(inputfeb,2);
  end
  if LIGHT_CHL
    chlb=mean(chlb,2);
  end
end

if ~periodicMatrix
  if rescaleForcing
	Rfs=mean(Rfs,2);    
  end  
end

% Biogeochem parameters
% Some constants for carbonate chemistry
rho0=1024.5; % kg/m^3 nominal density from OCMIP-2
PO4avg=2.17*rho0*1.0e-6; % mol/m^3
DOPavg=0.02*rho0*1.0e-6; % mol/m^3
O2avg=170.0*rho0*1.0e-6; % mol/m^3
DICavg=2230.0*rho0*1.0e-6; % mol/m^3
Alkavg=2370.0*rho0*1.0e-6; % mol/m^3
if ALLOW_FE
  Feavg=5.5e-7; % mean concentration from one of Steph's simulations
end

if rearrangeProfiles
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  Theta=Theta(Ir,:);
  Salt=Salt(Ir,:);
  hFacCb=hFacCb(Ir);
  recip_hFacCb=recip_hFacCb(Ir);
  if rescaleForcing
    Rfs=Rfs(Ir,:);
  end  
  Ib=find(izBox==1);
%
  Ip=Ip_post;
  Ir=Ir_post;
end  

% Initial condition
DIC=repmat(DICavg,[nb 1]);
Alk=repmat(Alkavg,[nb 1]);
PO4=repmat(PO4avg,[nb 1]);
O2=repmat(O2avg,[nb 1]);
DOP=repmat(DOPavg,[nb 1]);
Fe=repmat(Feavg,[nb 1]);

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
  writePetscBin('po4ini.petsc',PO4)
  writePetscBin('dopini.petsc',DOP)
  writePetscBin('o2ini.petsc',O2)
  writePetscBin('dicini.petsc',DIC)
  writePetscBin('alkini.petsc',Alk)
  if ALLOW_FE
    writePetscBin('feini.petsc',Fe)
  end
% Biogeochemical parameters
  if ~useAtmModel 
    if doAnthroCalc
      write_binary('TpCO2.bin',Tco2,'real*8')
      write_binary('pCO2atm.bin',pCO2atm,'real*8')
    end
  else  
    write_binary('pCO2atm_ini.bin',pCO2atm_ini,'real*8')
    write_binary('dA.bin',dab_surf,'real*8')    
    if useEmissions
      write_binary('Tem.bin',T_E,'real*8')
      write_binary('fossil_fuel_emissions.bin',Ev,'real*8')          
      write_binary('land_use_emissions.bin',Dv,'real*8')          
    end
    if useLandModel
      write_binary('land_ini.bin',[Cvini Csini Cvsini]','real*8')    
    end
  end
  write_binary('alpha.bin',alpha,'real*8')
  write_binary('rain_ratio.bin',rain_ratio,'real*8')
  if ~periodicForcing
	write_binary('wind.bin',windb,'real*8')
	write_binary('atmosp.bin',atmospb,'real*8')
	write_binary('silica.bin',silicab,'real*8')
	write_binary('fice.bin',Ficeb,'real*8')
	if ALLOW_FE
	  write_binary('inputfe.bin',inputfeb,'real*8')
	end
	if LIGHT_CHL
	  write_binary('chl.bin',chlb,'real*8')
	end	
  else
    for im=1:nm
	  write_binary(['wind_' sprintf('%02d',im-1)],windb(:,im),'real*8')
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')
	  write_binary(['silica_' sprintf('%02d',im-1)],silicab(:,im),'real*8')
	  write_binary(['fice_' sprintf('%02d',im-1)],Ficeb(:,im),'real*8')
	  if ALLOW_FE
		write_binary(['inputfe_' sprintf('%02d',im-1)],inputfeb(:,im),'real*8')
	  end
	  if LIGHT_CHL
		write_binary(['chl_' sprintf('%02d',im-1)],chlb(:,im),'real*8')
	  end	  
	end
  end
  if READ_PAR
	if ~periodicForcing
	  write_binary('par.bin',parb,'real*8')
	else
	  for im=1:nm
		write_binary(['par_' sprintf('%02d',im-1)],parb(:,im),'real*8')
	  end
	end        
  else   
    write_binary('latitude.bin',latb,'real*8')
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
  write_binary('drF.bin',nz,'int')  
  write_binary('drF.bin',dznom,'real*8',1)
  writePetscBin('hFacC.petsc',hFacCb)
  writePetscBin('recip_hFacC.petsc',recip_hFacCb)
  
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
