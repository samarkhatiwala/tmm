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

diagDir=pwd;

diagLogFile=fullfile(diagDir,'available_diagnostics.txt');

diagLog=readtable(diagLogFile,'Delimiter',',','ReadVariableNames',0);
diagLog=table2cell(diagLog);
numDiags=size(diagLog,1);

for id=1:numDiags
  diagName=diagLog{id,1};
  [Diag,x,y,z,T,diagUnit,diagLongName]=get_mobi_diagnostic(base_path,diagDir,diagName,1);  
%   eval([diagName '=Diag;']);
%   if numDims==2
%     save(diagName,diagName,'x','y','T','diagUnit')
%     write2netcdf([diagName '.nc'],Diag,x,y,[],T,diagName,diagUnit)
%   elseif numDims==3
%     save(diagName,diagName,'x','y','z','T','diagUnit')
%     write2netcdf([diagName '.nc'],Diag,x,y,z,T,diagName,diagUnit)
%   end          
end
