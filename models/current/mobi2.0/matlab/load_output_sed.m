base_path='/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault';
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')
load(boxFile,'izBox','nb')

Ib=find(izBox==1);
nbb=length(Ib);

sedmaskb=read_binary('sedmask.bin',[],'float64');
sedmaskb(sedmaskb==0)=NaN;
sedmask=matrixToGrid(sedmaskb,Ib,boxFile,gridFile);

% sedMixedBlockSize = nzmaxSed*numSedMixedTracers;
% lSedMixedSize = lNumProfiles*sedMixedBlockSize;
% ierr = VecCreate(PETSC_COMM_WORLD,&sedMixedTracers);CHKERRQ(ierr);
% ierr = VecSetSizes(sedMixedTracers,lSedMixedSize,PETSC_DECIDE);CHKERRQ(ierr);

nzmax=8;
numMixedTracers=20;

trsed=readPetscBinVec('sedmixed.petsc',1,-1);
trs=reshape(trsed(:,end),[nzmax*numMixedTracers nbb]); % nzmax*numMixedTracers x nbb
trs=reshape(trs,[nzmax numMixedTracers nbb]); % nzmax x numMixedTracers x nbb
trs=permute(trs,[3 1 2]); % nbb x nzmax x numMixedTracers
TRsedmixed=zeros([nx ny nzmax numMixedTracers]);
for itr=1:numMixedTracers
  TRsedmixed(:,:,:,itr)=matrixToGrid(trs(:,:,itr),Ib,boxFile,gridFile).*sedmask;
end

% sedBuriedBlockSize = ibmaxSed*numSedBuriedTracers;
% lSedBuriedSize = lNumProfiles*sedBuriedBlockSize;
% ierr = VecCreate(PETSC_COMM_WORLD,&sedBuriedTracers);CHKERRQ(ierr);
% ierr = VecSetSizes(sedBuriedTracers,lSedBuriedSize,PETSC_DECIDE);CHKERRQ(ierr);

ibmax=20;
numBuriedTracers=2;

trsed=readPetscBinVec('sedburied.petsc',1,-1);
trs=reshape(trsed(:,end),[ibmax*numBuriedTracers nbb]); % ibmax*numBuriedTracers x nbb
trs=reshape(trs,[ibmax numBuriedTracers nbb]); % ibmax x numBuriedTracers x nbb
trs=permute(trs,[3 1 2]); % nbb x ibmax x numBuriedTracers
for itr=1:numBuriedTracers
  TRsedburied(:,:,:,itr)=matrixToGrid(trs(:,:,itr),Ib,boxFile,gridFile).*sedmask;
end

