from petsc4py.PETSc cimport Vec, PetscVec, Viewer, PetscViewer, PetscErrorCode, Comm
from tmm4py.tmm4py_core cimport PetscStepTimer

cdef extern from * nogil:
    ctypedef long   PetscInt
    ctypedef double PetscReal
    ctypedef double PetscScalar

# include "PETSc/petscdef.pxi"
cdef extern from * nogil:
  ctypedef enum PetscBool:
      PETSC_FALSE
      PETSC_TRUE

# Add members that you want to expose below
cdef extern from "mops_forcing.h": 
    struct _p_ExternalForcingCtx:
      PetscInt efctxId;
      PetscInt stateId;
      PetscScalar *localph;
      PetscInt numpCO2atm_hist;
      PetscBool fixedAtmosCO2;
      PetscBool useAtmModel;
      PetscScalar pCO2atm;
      PetscScalar atmModelDeltaT;
      PetscScalar Focean;
      PetscScalar Foceanint;
      PetscStepTimer atmWriteTimer;
      PetscBool atmAppendOutput;
      PetscViewer atmfd;
      PetscInt atmfp;
      PetscScalar GRunoff;
      PetscScalar localFburial;
      PetscScalar Fburial;
      PetscInt burialSumSteps;
      PetscBool calcDiagnostics;
      PetscStepTimer diagTimer;
      PetscBool appendDiagnostics;
    
    ctypedef _p_ExternalForcingCtx* ExternalForcingContext "ExternalForcingContext"

