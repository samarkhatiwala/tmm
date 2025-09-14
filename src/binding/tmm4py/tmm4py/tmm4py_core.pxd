from petsc4py.PETSc cimport Vec, PetscVec, Viewer, PetscViewer, PetscErrorCode, PETSC_SUCCESS, PETSC_ERR_PYTHON, Comm

# Note: Cython will ignore the definitions in the block below. The correct definitions 
# should be picked up by the C compiler from $PETSC_DIR/$PETSC_ARCH/include/petscconf.h 
# and $PETSC_DIR/include/petscsystypes.h when it reads those header files
cdef extern from * nogil:
    ctypedef long   PetscInt
    ctypedef double PetscReal
    ctypedef double PetscScalar

# include "PETSc/petscdef.pxi"
cdef extern from * nogil:
  ctypedef enum PetscBool:
      PETSC_FALSE
      PETSC_TRUE

cdef extern from "tmm_timer.h":
    struct _p_StepTimer:
      PetscBool isInitialized;
      PetscBool fixedStep;
      PetscBool haveResetStartTimeStep;
      PetscInt count;
      PetscInt startTimeStep;
      PetscInt startTimeStepResetFreq;
      PetscInt numTimeSteps;
      PetscInt maxNumIntervals;
      PetscInt currInterval;
      PetscInt *timeIntervals;

    ctypedef _p_StepTimer* PetscStepTimer "StepTimer"

    struct _p_PeriodicTimer:
      PetscBool isInitialized;
      PetscScalar *tdp;    
      PetscScalar cyclePeriod;
      PetscInt numPerPeriod;
    
    ctypedef _p_PeriodicTimer* PetscPeriodicTimer "PeriodicTimer"

    struct _p_TimeDependentTimer:
      PetscBool isInitialized;
      PetscScalar *tdt;
      PetscInt numTimes;

    ctypedef _p_TimeDependentTimer* PetscTimeDependentTimer "TimeDependentTimer"

cdef extern from "tmm_forcing_utils.h":
    struct _p_PeriodicVec:
      PetscVec *up;
      PetscBool firstTime;
      PetscInt numPerPeriod;
    
    ctypedef _p_PeriodicVec* PetscPeriodicVec "PeriodicVec"

    struct _p_TimeDependentVec:
      PetscVec utd[2];
      PetscInt vecLength;
      PetscBool firstTime;
      PetscInt itcurr;
    
    ctypedef _p_TimeDependentVec* PetscTimeDependentVec "TimeDependentVec"

    struct _p_PeriodicArray:
      PetscScalar **up;
      PetscInt arrayLength;
      PetscInt numValsPerProfile;  
      PetscBool firstTime;
      PetscInt numPerPeriod;
    
    ctypedef _p_PeriodicArray* PetscPeriodicArray "PeriodicArray"

    struct _p_TimeDependentArray:
      PetscScalar *utd[2];
      PetscInt arrayLength;
      PetscInt numValsPerProfile;
      PetscBool firstTime;
      PetscInt itcurr;
    
    ctypedef _p_TimeDependentArray* PetscTimeDependentArray "TimeDependentArray"

# Add members that you want to expose below
cdef extern from "tmm.h":    
    struct _p_TMMState:
      PetscInt stateId;
      PetscInt numTracers;
      PetscVec *c;
      PetscVec *qf;
      PetscVec *qef;
      PetscVec *qrel;
      PetscVec *cbc;
      PetscVec *cbf;
      PetscBool isInitialized;
      PetscBool useExternalForcing;
      PetscBool useForcingFromFile;
      PetscBool usePrescribedBC;
      PetscBool applyExternalForcing;
      PetscBool applyForcingFromFile;
      PetscBool applyBC;
      PetscBool periodicForcing;
      PetscBool timeDependentForcing;
      PetscBool constantForcing;
      PetscBool periodicBC;
      PetscBool timeDependentBC;
      PetscBool constantBC;
      PetscBool doCalcBC;
      PetscBool useMonitor;
      PetscBool doMisfit;
      PetscBool relaxTracer;
      PetscBool isInitializedExternalForcing;
      PetscBool isInitializedCalcBC;
      PetscBool doOutput;
      PetscBool appendOutput;
      PetscBool writePickup;
      PetscBool doWriteBC;
      PetscBool doWriteQF;
      PetscBool doWriteQEF;
      PetscBool pickupFromFile;
      PetscBool doTimeAverage;
      PetscBool avgAppendOutput;
      PetscBool doExtraWrite;

    ctypedef _p_TMMState* PetscTMMState "TMMState"

cdef class TMMState:
    cdef PetscTMMState state
    cdef public:
      object finiex
      object fcalcex
      object fwriex
      object ffinex
      object freiex
      object finicbc
      object fcalccbc
      object fwricbc
      object ffincbc
      object freicbc
      object finimon
      object fcalcmon
      object fwrimon
      object ffinmon
      object finimis
      object fcalcmis
      object fwrimis
      object ffinmis
      object stateId
      object numTracers
      object config
      object tracerNames
      object c
      object qf
      object qef
      object qrel
      object cbc
      object cbf
      
cdef class StepTimer:
    cdef PetscStepTimer timer
    cdef public:
      object count
      object startTimeStep
      object numTimeSteps

cdef class PeriodicTimer:
    cdef PetscPeriodicTimer timer

cdef class TimeDependentTimer:
    cdef PetscTimeDependentTimer timer

cdef class PeriodicVec:
    cdef PetscPeriodicVec pvec
    cdef PetscVec c
    cdef PetscInt lDim
    cdef PetscInt frombuf
    cdef public:
      object arr

cdef class TimeDependentVec:
    cdef PetscTimeDependentVec pvec
    cdef PetscVec c
    cdef PetscInt lDim
    cdef PetscInt frombuf
    cdef public:
      object arr

cdef class PeriodicArray:
    cdef PetscPeriodicArray parr
    cdef PetscScalar *data_pointer
    cdef PetscInt lDim
    cdef PetscInt frombuf
    cdef public:
      object arr

cdef class TimeDependentArray:
    cdef PetscTimeDependentArray parr
    cdef PetscScalar *data_pointer
    cdef PetscInt lDim
    cdef PetscInt frombuf
    cdef public:
      object arr

# This is taken from PETSc.pyx; there seems to be no way to import them directly
# Also, the asReal/asInt and other functions not really necessary as the C and python data types are 
# generally the same and the inline functions don't do any type conversion and simply return the same 
# value. Nor is it strictly correct to use the *Real functions since what we want is *Scalar. But I 
# couldn't get those to work.
# Public definitions (for cimport into other extension modules)
cdef inline object toBool(PetscBool value):
    return True if value else False

cdef inline PetscBool asBool(object value) except? <PetscBool>0:
    return PETSC_TRUE if value else PETSC_FALSE

cdef inline object toInt(PetscInt value):
    return value

cdef inline PetscInt asInt(object value) except? -1:
    return value

cdef inline object toReal(PetscReal value):
    return value

cdef inline PetscReal asReal(object value) except? -1:
    return value

# Public functions (for cimport into other extension modules)
cdef getNumpyRealArrayFromPointer(PetscScalar *data_pointer, PetscInt n)
cdef getNumpyIntArrayFromPointer(PetscInt *data_pointer, PetscInt n)
cdef copyNumpyRealArrayFromPointer(PetscScalar *data_pointer, PetscInt n)
cdef copyNumpyIntArrayFromPointer(PetscInt *data_pointer, PetscInt n)