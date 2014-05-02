typedef struct {
  PetscScalar    *aa;
  PetscInt       *jj;
  PetscInt       Istart,Iend,*numColsPerLocRow;
  PetscInt       *allIstart,*allnumNNZ;
  off_t          off1,off2;
  PetscTruth     firstTime;
  PetscTruth     lastTime;
} MatLayout;

extern PetscErrorCode VecAXPBYmy(PetscScalar a,PetscScalar b,Vec x,Vec y,Vec *z);
extern PetscErrorCode MatAXPBYmy(PetscScalar a,PetscScalar b,Mat X,Mat Y,Mat *Z);
extern PetscErrorCode VecLoadIntoVectorRandomAccess(PetscViewer viewer,Vec vec, PetscInt length, PetscInt iRec);
extern PetscErrorCode VecLoadVecIntoArray(Vec v, const char filename[], PetscScalar *arr);
extern PetscErrorCode MatLoadIntoMatrix(PetscViewer viewer, Mat A);
extern PetscErrorCode MatLoadIntoMatrix2(const char filename[], Mat A);
extern PetscErrorCode MatLoadIntoMatrix3(const char filename[], Mat A);
extern PetscErrorCode MatLoadIntoMatrix4(const char filename[], Mat A, MatLayout *user);
extern PetscErrorCode MatGetSizeFromFile(const char filename[], PetscInt *M, PetscInt *N, PetscInt *nnz);
