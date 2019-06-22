base_path='/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault';

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')

load(boxFile,'izBox','nb','volb')

load(profilesFile,'Irr')

Ib=find(izBox==1);
nbb=length(Ib);

diagLogFile='available_diagnostics.txt';

diagLog=readtable(diagLogFile,'Delimiter',',','ReadVariableNames',0);
diagLog=table2cell(diagLog);

load('diagnostic_output_time.txt')
T=diagnostic_output_time(:,2);

numDiags=size(diagLog,1);

for id=1:numDiags
  diagName=diagLog{id,1};
  diagLongName=diagLog{id,2};
  diagUnit=diagLog{id,3};
  diagFile=diagLog{id,4};
  if ~isempty(strfind(diagFile,'bin'))
    numDims=2;
  elseif ~isempty(strfind(diagFile,'petsc'))
    numDims=3;
  else
    error('Unknown file type!')
  end
  if numDims==2
    diagData=read_binary(diagFile,[],'float64');
    nt=length(diagData)/nbb;
    diagData=reshape(diagData,[nbb nt]);
    Diag=matrixToGrid(diagData,Ib,boxFile,gridFile);
  elseif numDims==3
    diagData=readPetscBinVec(diagFile,-1);
    Diag=matrixToGrid(diagData(Irr,:),[],boxFile,gridFile);    
  end
  eval([diagName '=Diag;']);
  if numDims==2
    save(diagName,diagName,'x','y','T','diagUnit')
    write2netcdf([diagName '.nc'],Diag,x,y,[],T,diagName,diagUnit)
  elseif numDims==3
    save(diagName,diagName,'x','y','z','T','diagUnit')
    write2netcdf([diagName '.nc'],Diag,x,y,z,T,diagName,diagUnit)
  end          
end
