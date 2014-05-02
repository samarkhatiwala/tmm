TRmm=readPetscBinVec('trmm.petsc',-1);
TRmmg=grid_boxes3d(TRmm(Irr,:),[],boxFile,gridFile);
TReqmm=readPetscBinVec('TReqavg.petsc',-1);
TReqmmg=grid_boxes3d(TReqmm(Irr,:),[],boxFile,gridFile);
TRsatanommm=readPetscBinVec('TRsatanomavg.petsc',-1);             
TRsatanommmg=grid_boxes3d(TRsatanommm(Irr,:),[],boxFile,gridFile);

TRavgg=mean(TRmmg,4);
TReqavgg=mean(TReqmmg,4);
TRsatanomavgg=mean(TRsatanommmg,4);

