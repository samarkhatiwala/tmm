riverFile='perry1996-runoff_noname.txt';
% riverFile='perry1996-runoff-noarctic_noname.txt'; % for testing

%%%%
% mindist = sqrt(2)*2*2.8;
% riverFile='perry1996-runoff-noarctic_noname.txt';
% mindist = sqrt(2)*2*.9;

% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','dx','dy','gridType')

load(boxFile,'Xboxnom','Yboxnom','izBox','nb','volb')

Ib=find(izBox==1);
nbb=length(Ib);

load(profilesFile,'Ip_pre','Ir_pre')
Ip=Ip_pre;
Ir=Ir_pre;

% figure out grid spacing
a=6371e3;
dxb=gridToMatrix(dx,Ib,boxFile,gridFile);
% dyb=gridToMatrix(dy,Ib,boxFile,gridFile);
th=Yboxnom(Ib)*pi/180; % latitude in radians
dphib=(180/pi)*dxb./(a*cos(th)); % dx in degrees
% dth=(180/pi)*dyb/a; % dy in degrees
mindist = sqrt(2)*2*median(dphib);

% ! regrid perry runoff onto MIT2.8 degree grid
% ! modifications/methods:
% ! (1) restrict Perry's data set to rivers south of 60N
% ! (2) define a large plume for the five largest rivers
% ! (3) define a small plume/point source for the other rivers

% ! Units:
% ! vol  : km3/a

[hdr,rivdata]=hdrload(riverFile);
ptnum=rivdata(:,1);
ptlat=rivdata(:,2);
ptlon=rivdata(:,3);
prunoff=rivdata(:,4);
nriv=length(prunoff);

% ! Note: Perry's longitude starts with 0 at Greenwich, and then continues
% ! Note: Perry's runoff is in m3/sec
% ! -->   Scale to km3/a

% ! DO NOT ACCOUNT FOR ARCTIC RIVERS
plon=ptlon;
plat=ptlat;

perryvol=prunoff*1e-9*86400*360;
perrytotvol=sum(perryvol);

modelsvol=volb(Ib)*1e-9;
modelsvol=repmat(modelsvol,[1 nriv]);

modelx=repmat(Xboxnom(Ib),[1 nriv]);
modely=repmat(Yboxnom(Ib),[1 nriv]);

% ! DISTANCES BETWEEN MODEL LOCATION AND OBSERVED LOCATIOS
plonext=repmat(plon',[nbb 1]);
platext=repmat(plat',[nbb 1]);

delx = abs(plonext - modelx);
dely = abs(platext - modely);

% ! EUCLIDEAN DISTANCE IN TERMS OF LAT AND LON
deleuc = sqrt(delx.^2+dely.^2);

% ! MIN. EUCLIDEAN DISTANCE AS SCALAR FOR EACH L 
mindeleuc = min(deleuc,[],1);

% ! GET THE NEAREST LOCATION OF EACH OCEAN WET BOX, 
% ! BUT ONLY IF IT IS NOT FURTHER AWAY THEN 10 DEGREES
lpnearcount=zeros(size(deleuc));
lpnearvol=zeros(size(deleuc));

for ir=1:nriv
  kr=find(deleuc(:,ir)==mindeleuc(ir) & deleuc(:,ir)<mindist);
  lpnearcount(kr,ir)=1;
  lpnearvol(kr,ir)=perryvol(ir);
end

% ! SCALE THE DISTANCE BY THE AMOUNT OF RUNOFF
scaledist = lpnearvol./modelsvol*2.0;

% ! SCALAR FOR EACH L 
scalar_lpnearcount = sum(lpnearcount,1);
scalar_lpnearvol = sum(lpnearvol,1);
scalar_scaledist = max(1,sum(scaledist,1));

% ! IF THERE IS AT LEAST ONE GOOD VALUE FOR THIS RIVER, THEN MARK ALL FIELDS WITHIN A GIVEN CIRCLE AROUND LOCATION
lpcount=zeros(size(deleuc));
for ir=1:nriv
  if scalar_lpnearcount(ir)>0
    kr=find(deleuc(:,ir) <= mindeleuc(ir)*scalar_scaledist(ir));
    lpcount(kr,ir)=1;
  end
end  

lpweight=zeros(size(lpcount));
scalar_lpcount = sum(lpcount,1);
for ir=1:nriv
  if scalar_lpcount(ir)>0
    lpweight(:,ir) = lpcount(:,ir)/scalar_lpcount(ir);
  end
end

countplon=plon;
countplat=plat;
discardplon=plon;
discardplat=plat;
countplon(scalar_lpcount==0)=NaN;
countplat(scalar_lpcount==0)=NaN;
discardplon(scalar_lpcount>0)=NaN;
discardplat(scalar_lpcount>0)=NaN;

% ! IF THERE IS AT LEAST ONE GOOD VALUE FOR THIS RIVER, THEN FILL ALL FIELDS WITHIN A GIVEN CIRCLE AROUND LOCATION
lpvol=zeros(size(lpweight));
for ir=1:nriv
  if scalar_lpcount(ir)>0
    lpvol(:,ir) = lpweight(:,ir)*scalar_lpnearvol(ir);
  end
end

lptotvol = sum(lpvol(:));

fracvolint = sum(lpvol,2)/lptotvol;

%
pfracvol=fracvolint;

vvol=zeros(nbb,1);
for is=1:nbb % loop over each surface point
  Ipl=Ip{is}; % indices for local profile (globally indexed)
  vvol(is)=sum(volb(Ipl));
end
vvol=vvol*1e-12;

% ! THE FOLLOWING GIVES A 3D FIELD OF GAIN PER VOLUME AND YEAR AT
% ! LOCATIONS AFFECTED BY RIVER INPUT.

tmp=pfracvol./vvol;
gaine12=zeros(nb,1);
for is=1:nbb % loop over each surface point
  Ipl=Ip{is}; % indices for local profile (globally indexed)  
  gaine12(Ipl)=tmp(is);
end
gaine12g=matrixToGrid(gaine12,[],boxFile,gridFile);

gaine12skipl=mean(gaine12,2);

gaine12skipl=gaine12skipl(Ir);

writePetscBin('runoff_volume_annual.petsc',gaine12skipl)
