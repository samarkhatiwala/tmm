base_path='/Volumes/data2/spk/TransportMatrixConfigs/UVicOSUpicdefault';

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')
load(boxFile,'nb','izBox')

Ib=find(izBox==1);
Ii=find(izBox~=1);

nbb=length(Ib);
nbi=length(Ii);

trNames={'G_0001','G_0002','G_0003','G_0004','G_0005','G_0006','G_0007'};
bcNames={'Gbc_out_0001','Gbc_out_0002','Gbc_out_0003','Gbc_out_0004','Gbc_out_0005','Gbc_out_0006','Gbc_out_0007'};

numTr=length(trNames);

if strcmp(gridType,'llc_v4')
  load(fullfile(base_path,'llc_v4_grid'))
  gcmfaces_global
end

for itr=1:numTr
  varName=upper(trNames{itr})
  fn=[bcNames{itr} '.petsc'];
  bc=readPetscBinVec(fn,-1);  
  fn=[trNames{itr} '.petsc'];
  tmptr=readPetscBinVec(fn,-1);
  nt=size(tmptr,2);
  tr=zeros([nb nt]);
  tr(Ib,:)=bc;
  tr(Ii,:)=tmptr;
  TR=matrixToGrid(tr,[],boxFile,gridFile);

  if strcmp(gridType,'llc_v4')
    varName=[varName '_plot'];
	tmp=gcmfaces(TR);
	[x,y,TRplot]=convert2pcol(mygrid.XC,mygrid.YC,tmp);
	[n1,n2]=size(TRplot);
	eval([varName '=zeros([n1 n2 nz nt]);']);
	for it=1:nt
	  for iz=1:nz
		eval(['[x,y,' varName '(:,:,iz,it)]=convert2pcol(mygrid.XC,mygrid.YC,tmp(:,:,iz,it));']);
	  end
	end  
  else
	eval([varName '=TR;']);
  end
  save(varName,varName,'x','y','z')
%  write2netcdf([varName '.nc'],TR,x,y,z,[],upper(varName))
end
