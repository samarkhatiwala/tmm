% Script to define regions
writeRegionFiles=1
regionFile='global_plus_regions';

% Set toplevel path to GCMs configuration
base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/disks/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/disks/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
basinFile=fullfile(base_path,'GCM','basin_mask');

load(gridFile,'nx','ny','nz','x','y')

load(boxFile,'nb','izBox','Xboxnom','Yboxnom')

load(basinFile,'basin_mask');

Ib=find(izBox==1);
nbb=length(Ib);
Xbb=Xboxnom(Ib);
Ybb=Yboxnom(Ib);

basinsb=gridToMatrix(basin_mask,Ib,boxFile,gridFile);

% region map
numRegions=7;
Irb=zeros([nbb numRegions]);
for ir=1:numRegions
  switch ir
    case 1 % global
      Ireg=[1:nbb]';
    case 2 % NA
      Ireg=find(basinsb==1 & Ybb>=40);
    case 3 % ANT
      Ireg=find(Ybb<-50);
    case 4 % SUBANT
      Ireg=find(Ybb>=-50 & Ybb<=-40);
    case 5 % STROP
      Ireg=find((Ybb>-40 & Ybb<=-30) | (Ybb>=30 & Ybb<40));
    case 6 % TROP
      Ireg=find(Ybb>-30 & Ybb<30);
    case 7 % NPAC
      Ireg=find(basinsb==3 & Ybb>=40);
    otherwise
      error('ERROR: unknown region!')
  end  
  Irb(Ireg,ir)=ir;
end

regionNumber=[1:numRegions]'; 
regionb=Irb;

regionmap=matrixToGrid(Irb,Ib,boxFile,gridFile);
      
if writeRegionFiles
  save(regionFile,'numRegions','regionmap','regionb','regionNumber','x','y');
end
