% base_path='/Volumes/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
base_path='/Volumes/data2/spk/TransportMatrixConfigs/UVicOSUpicdefault';

rearrangeProfiles=0

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')
load(boxFile,'nb','izBox')

if rearrangeProfiles
  load(profilesFile,'Irr')
end

trNames={'C'};

numTr=length(trNames);

if strcmp(gridType,'llc_v4')
  load(fullfile(base_path,'llc_v4_grid'))
  gcmfaces_global
end

for itr=1:numTr
  varName=upper(trNames{itr})
  fn=[trNames{itr} '.petsc'];
  tr=readPetscBinVec(fn,1,-1);
  if rearrangeProfiles
    tr=tr(Irr,end);
  end  
  TR=matrixToGrid(tr,[],boxFile,gridFile);
  
  if strcmp(gridType,'llc_v4')
    varName=[varName '_plot'];
	tmp=gcmfaces(TR);
	[x,y,TRplot]=convert2pcol(mygrid.XC,mygrid.YC,tmp);
	[n1,n2]=size(TRplot);
	eval([trNames{itr} '=zeros([n1 n2 nz]);']);
	for iz=1:nz
	  eval(['[x,y,' varName '(:,:,iz)]=convert2pcol(mygrid.XC,mygrid.YC,tmp(:,:,iz));']);
	end
  else
	eval([varName '=TR;']);
  end

  save(varName,varName,'x','y','z')
%   write2netcdf([varName '.nc'],TR,x,y,z,[],upper(varName))
end
