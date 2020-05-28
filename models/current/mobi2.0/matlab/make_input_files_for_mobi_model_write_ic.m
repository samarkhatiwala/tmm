% Set toplevel path to GCMs configuration
base_path='/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

%
gridFile=fullfile(base_path,'grid');

load(gridFile,'nz','dznom','z')
  
% Grid data
z=z*100; % convert to cm
dznom=dznom*100; % convert to cm
write_binary('zt.bin',nz,'int')  
write_binary('zt.bin',z,'real*8',1)
write_binary('drF.bin',nz,'int')  
write_binary('drF.bin',dznom,'real*8',1)
