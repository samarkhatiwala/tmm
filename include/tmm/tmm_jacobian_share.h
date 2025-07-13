#ifdef DEFINE_JACOBIAN_VARIABLES
#define EXTERNJACSHARE
#else
#define EXTERNJACSHARE extern
#endif

typedef struct {
  PetscInt numTracers;
  PetscScalar tc, tf;
  PetscInt Iterc, iLoop;
  PetscScalar *XTovScaleFac;
  PetscScalar *vToXScaleFac;  
  TMMState state;
} TMMJAC;

EXTERNJACSHARE PetscInt lXSize;
EXTERNJACSHARE PetscInt *gXIndices, gXLow, gXHigh;
