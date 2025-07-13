# from __future__ import absolute_import
cimport cython

import petsc4py
import sys
import warnings
import atexit

try:
  from mpi4py import MPI
except:
  MPI=None

# This doesn't give an error on compiling but it fails as runtime
# from petsc4py import asBool, asInt, asReal, asScalar, toBool, toInt, toReal

from petsc4py import PETSc

cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
from cpython cimport PyObject, Py_INCREF
from cpython.mem cimport PyMem_Malloc, PyMem_Free

np.import_array()

from cpython.bytes cimport PyBytes_AsString

# def mywarnings(message, category, filename, lineno, file=None, line=None):
#   print(message)
# 
# warnings.showwarning=mywarnings

debug=False

# we use this to keep track of all created states so that we can call their destroy method
_statescreated=[]

cdef extern from "tmm_share.h":
    PetscScalar deltaTClock
    PetscScalar time0
    PetscInt maxSteps
    PetscInt Iter0
    PetscInt nb
    PetscInt *lProfileLength
    PetscInt *lStartIndices
    PetscInt *lEndIndices
    PetscInt lBCSize, gBCSize
    PetscBool prescribedBCInUse, calcBCInUse;
    PetscInt lNumProfiles, lSize, numPrevProfiles, totalNumProfiles
    PetscBool useProfiles

configVars = ["useExternalForcing", 
              "useForcingFromFile", 
              "usePrescribedBC", 
              "applyExternalForcing", 
              "applyForcingFromFile", 
              "applyBC", 
              "periodicForcing", 
              "timeDependentForcing", 
              "constantForcing", 
              "periodicBC", 
              "timeDependentBC", 
              "constantBC", 
              "doCalcBC", 
              "useMonitor", 
              "doMisfit", 
              "relaxTracer", 
              "isInitializedExternalForcing",
              "isInitializedCalcBC",
              "doOutput", 
              "appendOutput", 
              "writePickup", 
              "doWriteBC", 
              "doWriteQF", 
              "doWriteQEF", 
              "pickupFromFile", 
              "doTimeAverage", 
              "avgAppendOutput", 
              "doExtraWrite",
              "numTracers"] 

profileVars = ["useProfiles", 
               "lSize",
               "lNumProfiles",
               "totalNumProfiles",
               "lProfileLength",
               "lStartIndices",
               "lEndIndices"]               

cdef extern from * nogil:
    ctypedef PetscErrorCode (*TMMExtIniFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMExtCalcFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMExtWriFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMExtFinFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMExtReiFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx)

    ctypedef PetscErrorCode (*TMMCalcBCIniFunctionFn_type)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMCalcBCCalcFunctionFn_type)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMCalcBCWriFunctionFn_type)(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMCalcBCFinFunctionFn_type)(PetscScalar tc, PetscInt Iterc, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMCalcBCReiFunctionFn_type)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscTMMState state, void *ctx)

    ctypedef PetscErrorCode (*TMMMonIniFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMonCalcFunctionFn_type)(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMonWriFunctionFn_type)(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMonFinFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)

    ctypedef PetscErrorCode (*TMMMisIniFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMisCalcFunctionFn_type)(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMisWriFunctionFn_type)(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx)
    ctypedef PetscErrorCode (*TMMMisFinFunctionFn_type)(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx)

    PetscErrorCode TMMSetIniExtForcFunction(PetscTMMState state, TMMExtIniFunctionFn_type, void*)
    PetscErrorCode TMMSetCalcExtForcFunction(PetscTMMState state, TMMExtCalcFunctionFn_type, void*)
    PetscErrorCode TMMSetWriExtForcFunction(PetscTMMState state, TMMExtWriFunctionFn_type, void*)
    PetscErrorCode TMMSetFinExtForcFunction(PetscTMMState state, TMMExtFinFunctionFn_type, void*)
    PetscErrorCode TMMSetReiExtForcFunction(PetscTMMState state, TMMExtReiFunctionFn_type, void*)

    PetscErrorCode TMMSetIniCalcBCFunction(PetscTMMState state, TMMCalcBCIniFunctionFn_type, void*)
    PetscErrorCode TMMSetCalcCalcBCFunction(PetscTMMState state, TMMCalcBCCalcFunctionFn_type, void*)
    PetscErrorCode TMMSetWriCalcBCFunction(PetscTMMState state, TMMCalcBCWriFunctionFn_type, void*)
    PetscErrorCode TMMSetFinCalcBCFunction(PetscTMMState state, TMMCalcBCFinFunctionFn_type, void*)
    PetscErrorCode TMMSetReiCalcBCFunction(PetscTMMState state, TMMCalcBCReiFunctionFn_type, void*)

    PetscErrorCode TMMSetIniMonitorFunction(PetscTMMState state, TMMMonIniFunctionFn_type, void*)
    PetscErrorCode TMMSetCalcMonitorFunction(PetscTMMState state, TMMMonCalcFunctionFn_type, void*)
    PetscErrorCode TMMSetWriMonitorFunction(PetscTMMState state, TMMMonWriFunctionFn_type, void*)
    PetscErrorCode TMMSetFinMonitorFunction(PetscTMMState state, TMMMonFinFunctionFn_type, void*)

    PetscErrorCode TMMSetIniMisfitFunction(PetscTMMState state, TMMMisIniFunctionFn_type, void*)
    PetscErrorCode TMMSetCalcMisfitFunction(PetscTMMState state, TMMMisCalcFunctionFn_type, void*)
    PetscErrorCode TMMSetWriMisfitFunction(PetscTMMState state, TMMMisWriFunctionFn_type, void*)
    PetscErrorCode TMMSetFinMisfitFunction(PetscTMMState state, TMMMisFinFunctionFn_type, void*)

cdef extern from "tmm.h":
    PetscErrorCode TMMInitialize(PetscInt *Iter0ret, PetscInt *maxStepsret, PetscScalar *time0ret, PetscScalar *deltaTClockret)
    PetscErrorCode TMMCreate(PetscTMMState *state)
    PetscErrorCode TMMSetFromOptions(PetscTMMState state, const char pre[], PetscBool doOutput)
    PetscErrorCode TMMUpdateTMs(PetscScalar tc)
    PetscErrorCode TMMForcingUpdate(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state)
    PetscErrorCode TMMTimeStep(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state)
    PetscErrorCode TMMTimeStepPost(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state)
    PetscErrorCode TMMOutput(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state)
    PetscErrorCode TMMDestroy(PetscScalar tc, PetscTMMState *state)
    PetscErrorCode TMMFinalize(PetscScalar tc)

cdef extern from "tmm_timer.h":
    PetscErrorCode StepTimerCreate(PetscStepTimer *timer);
    PetscErrorCode PeriodicTimerCreate(PetscPeriodicTimer *timer);
    PetscErrorCode TimeDependentTimerCreate(PetscTimeDependentTimer *timer);
    PetscErrorCode StepTimerIni(const char pre1[], const char pre2pre[], PetscInt Iter0, PetscStepTimer thetimer);
    PetscErrorCode StepTimerUpdate(PetscInt Iter, PetscStepTimer thetimer);
    PetscErrorCode PeriodicTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdp, PetscPeriodicTimer thetimer);
    PetscErrorCode TimeDependentTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdt, PetscTimeDependentTimer thetimer);

cdef extern from "tmm_petsc_matvec_utils.h":
    PetscErrorCode VecLoadIntoArray(PetscInt lDim, const char filename[], PetscScalar *arr);
    PetscErrorCode VecLoadIntoVectorRandomAccessFromFile(const char filename[], PetscVec vec, PetscInt length, PetscInt iRec);
    PetscErrorCode dotProdArrays(PetscScalar *xarr, PetscScalar *yarr, PetscInt n, PetscScalar *z);
    PetscErrorCode sumArray(PetscScalar *xarr, PetscInt n, PetscScalar *z);
    PetscErrorCode sumScalar(PetscScalar x, PetscScalar *tot);
    PetscErrorCode Barrier();
#     PetscErrorCode VecWriteLocalArrayToVec(PetscInt lDim, PetscScalar *arr, const char filename[], PetscFileMode mode);

cdef extern from "tmm_forcing_utils.h":
    PetscErrorCode PeriodicVecCreate(PetscPeriodicVec *c);
    PetscErrorCode PeriodicVecInterp(PetscScalar tc, PetscVec *u, PetscScalar cyclePeriod, PetscInt numPerPeriod, PetscScalar *tdp, PetscPeriodicVec c, const char *fileName);
    PetscErrorCode PeriodicVecDestroy(PetscPeriodicVec *c);
    PetscErrorCode VecCreateFromLocalSize(PetscInt lDim, PetscVec *c);
    PetscErrorCode VecCreateWithArrayFromLocalSize(PetscInt lDim, PetscScalar *arr, PetscVec *c);
    PetscErrorCode PeriodicArrayCreate(PetscPeriodicArray *arr, PetscInt arrayLength);
    PetscErrorCode PeriodicArrayDestroy(PetscPeriodicArray *arr);
    PetscErrorCode interpPeriodicProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PetscPeriodicArray user, const char *fileName);
    PetscErrorCode TimeDependentVecCreate(PetscTimeDependentVec *c);
    PetscErrorCode TimeDependentVecInterp(PetscScalar tc, PetscVec *u, PetscInt N, 
                                                PetscScalar *tdp, PetscTimeDependentVec c, const char *fileName);
    PetscErrorCode TimeDependentVecDestroy(PetscTimeDependentVec *c);
    PetscErrorCode TimeDependentArrayCreate(PetscTimeDependentArray *arr, PetscInt arrayLength);
    PetscErrorCode TimeDependentArrayDestroy(PetscTimeDependentArray *arr);
    PetscErrorCode writeBinaryScalarData(const char *fileName, PetscScalar *arr, PetscInt N, PetscBool appendToFile);

cdef extern from "tmm_profile_utils.h":
    PetscErrorCode readProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile);
    PetscErrorCode writeProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscBool appendToFile);
    PetscErrorCode dotProdProfileSurfaceScalarData(PetscScalar *xarr, PetscScalar *yarr, PetscScalar *z);
    PetscErrorCode interpTimeDependentProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscInt numTimes, PetscScalar *tdt,
                                    PetscTimeDependentArray user, const char *fileName);
    
cdef extern from * nogil:
    PetscErrorCode VecRestoreArray(PetscVec,PetscScalar*[])
    PetscErrorCode VecDuplicate(PetscVec,PetscVec*)
    PetscErrorCode VecCopy(PetscVec,PetscVec)
    PetscErrorCode VecAssemblyBegin(PetscVec)
    PetscErrorCode VecAssemblyEnd(PetscVec)
    PetscErrorCode VecGetArray(PetscVec,PetscScalar*[])
    PetscErrorCode VecGetLocalSize(PetscVec,PetscInt*)
    PetscErrorCode VecSet(PetscVec,PetscScalar)
    PetscErrorCode VecRestoreArray(PetscVec,PetscScalar*[])
    PetscErrorCode VecDestroy(PetscVec *c)

# Callback functions: these are the functions called by TMMComputeExtForcFunction and 
# TMMComputeCalcBCFunction; the actual PYTHON callback functions, arguments to it, and 
# the PYTHON State object are stashed in ctx and unpacked below. The python functions 
# can - if they are wrapped in Cython - in turn call C functions.
# Note: the *args and **kwargs unpack a tuple and dict, respectively, into arguments that 
# are passed to the python callback function.
cdef PetscErrorCode iniExternalForcingCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode calcExternalForcingCB(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode wriExternalForcingCB(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode finExternalForcingCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode reiExternalForcingCB(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode iniCalcBCCB(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iterc), toReal(tf), toInt(Iterf), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode calcCalcBCCB(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iterc), toReal(tf), toInt(Iterf), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode wriCalcBCCB(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iterc), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode finCalcBCCB(PetscScalar tc, PetscInt Iterc, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iterc), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode reiCalcBCCB(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iterc), toReal(tf), toInt(Iterf), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode iniMonitorCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode calcMonitorCB(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode wriMonitorCB(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode finMonitorCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode iniMisfitCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode calcMisfitCB(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode wriMisfitCB(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(iLoop), State, *args, **kwargs)
    return PETSC_SUCCESS

cdef PetscErrorCode finMisfitCB(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx) noexcept with gil:
    (function, args, kwargs, State) = <object>ctx
    (<object>function)(toReal(tc), toInt(Iter), State, *args, **kwargs)
    return PETSC_SUCCESS

# This is to prevent someone from accidentally modifying TMMState.c, etc. 
# You can still modify the numpy arrays in that list but not replace them 
# with something else, or append or insert something in the list.
class _ImmutableList(list):
  def __setitem__(self, index, item):
    raise RuntimeError("List is not mutable")
  def insert(self, index, item):
    raise RuntimeError("List is not insertable")
  def append(self, item):
   raise RuntimeError("List is not appendable")

class Stub:

    def __init__(self, Id, prefix):
      self.Id = Id
    def iniExternalForcingFn(self, tc, Iter, state):
      pass
    def calcExternalForcingFn(self, tc, Iter, iLoop, state):
      pass
    def writeExternalForcingFn(self, tc, Iter, iLoop, state):
      pass
    def finalizeExternalForcingFn(self, tc, Iter, state):
      pass
    def reInitializeExternalForcingFn(self, tc, Iter, iLoop, state):
      pass
    def iniCalcBCFn(self, tc, Iterc, tf, Iterf, state):
      pass
    def calcCalcBCFn(self, tc, Iterc, tf, Iterf, iLoop, state):
      pass
    def writeCalcBCFn(self, tc, Iter, iLoop, state):
      pass
    def finalizeCalcBCFn(self, tc, Iter, state):
       pass
    def reInitializeCalcBCFn(self, tc, Iter, iLoop, state):
      pass
    def iniMonitorFn(self, tc, Iter, state):
      pass
    def calcMonitorFn(self, tc, iLoop, state):
      pass
    def writeMonitorFn(self, tc, iLoop, state):
      pass
    def finalizeMonitorFn(self, tc, Iter, state):
      pass
    def iniMisfitFn(self, tc, Iter, state):
      pass
    def calcMisfitFn(self, tc, iLoop, state):
      pass
    def writeMisfitFn(self, tc, iLoop, state):
      pass
    def finalizeMisfitFn(self, tc, Iter, state):
      pass

cdef class TMMState():

#     cdef PetscTMMState state
    def __cinit__(self):
        self.state = NULL
      
    def __init__(self):
        self.finiex = None
        self.fcalcex = None
        self.fwriex = None
        self.ffinex = None
        self.freiex = None
        self.finicbc = None
        self.fcalccbc = None
        self.fwricbc = None
        self.ffincbc = None
        self.freicbc = None
        self.finimon = None
        self.fcalcmon = None
        self.fwrimon = None
        self.ffinmon = None
        self.finimis = None
        self.fcalcmis = None
        self.fwrimis = None
        self.ffinmis = None
        self.stateId = -1
        self.numTracers = 0
        self.config = None
        self.tracerNames = None
        self.c = None
        self.qf = None
        self.qef = None
        self.qrel = None
        self.cbc = None
        self.cbf = None
        _statescreated.append(self)
        
    def create(self):
        
        cdef PetscErrorCode ierr 
        
        cdef PetscTMMState newstate
        ierr = TMMCreate(&newstate)
        if ierr != 0: raise PETSc.Error(ierr)
        self.state=newstate

        return self

#   In the following setIni* functions we pack the name of the actual python 
#   callback function, its arguments and the python state object in 'context'. 
#   This is passed as a 'user context' argument to the TMMSetIni* functions, 
#   which stash a pointer to it in the PETSc state object. However, when the 
#   setIni* functions finish, context will go out of scope and be destroyed 
#   (along with its pointer). To prevent this, we also store a pointer to context 
#   in the python state object. 
#   The TMMCompute* functions will retrieve the pointer to context from the 
#   PETSc state object and pass it to the *ForcingFn callback functions, which 
#   in turn will unpack context into the actual python callback function, its 
#   arguments and the python state object.

    def setIniExternalForcingFunction(self, func, *args, **kwargs):
        if not self.config['useExternalForcing']:
            raise RuntimeError("setIniExternalForcingFunction failed: -external_forcing not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.finiex = context # this saves a pointer to context
            if debug: print(f"Ini exf forcing func is {context}")
            ierr = TMMSetIniExtForcFunction(self.state, iniExternalForcingCB, <void*>context)

        else:
            ierr = TMMSetIniExtForcFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetIniExtForcFunction failed with error code {ierr}")

    def setCalcExternalForcingFunction(self, func, *args, **kwargs):
        if not self.config['useExternalForcing']:
            raise RuntimeError("setCalcExternalForcingFunction failed: -external_forcing not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fcalcex = context # this saves a pointer to context
            ierr = TMMSetCalcExtForcFunction(self.state, calcExternalForcingCB, <void*>context)
        else:
            ierr = TMMSetCalcExtForcFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetCalcExtForcFunction failed with error code {ierr}")

    def setWriExternalForcingFunction(self, func, *args, **kwargs):
        if not self.config['useExternalForcing']:
            raise RuntimeError("setWriExternalForcingFunction failed: -external_forcing not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fwriex = context # this saves a pointer to context
            ierr = TMMSetWriExtForcFunction(self.state, wriExternalForcingCB, <void*>context)
        else:
            ierr = TMMSetWriExtForcFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetWriExtForcFunction failed with error code {ierr}")

    def setFinExternalForcingFunction(self, func, *args, **kwargs):
        if not self.config['useExternalForcing']:
            raise RuntimeError("setFinExternalForcingFunction failed: -external_forcing not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.ffinex = context # this saves a pointer to context
            ierr = TMMSetFinExtForcFunction(self.state, finExternalForcingCB, <void*>context)
        else:
            ierr = TMMSetFinExtForcFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetFinExtForcFunction failed with error code {ierr}")

    def setReiExternalForcingFunction(self, func, *args, **kwargs):
        if not self.config['useExternalForcing']:
            raise RuntimeError("setReiExternalForcingFunction failed: -external_forcing not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.freiex = context # this saves a pointer to context
            ierr = TMMSetReiExtForcFunction(self.state, reiExternalForcingCB, <void*>context)
        else:
            ierr = TMMSetReiExtForcFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetReiExtForcFunction failed with error code {ierr}")

    def setIniCalcBCFunction(self, func, *args, **kwargs):
        if not self.config['doCalcBC']:
            raise RuntimeError("setIniCalcBCFunction failed: -calc_bc not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.finicbc = context # this saves a pointer to context
            ierr = TMMSetIniCalcBCFunction(self.state, iniCalcBCCB, <void*>context)
        else:
            ierr = TMMSetIniCalcBCFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetIniCalcBCFunction failed with error code {ierr}")

    def setCalcCalcBCFunction(self, func, *args, **kwargs):
        if not self.config['doCalcBC']:
            raise RuntimeError("setCalcCalcBCFunction failed: -calc_bc not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fcalccbc = context # this saves a pointer to context
            ierr = TMMSetCalcCalcBCFunction(self.state, calcCalcBCCB, <void*>context)
        else:
            ierr = TMMSetCalcCalcBCFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetCalcCalcBCFunction failed with error code {ierr}")

    def setWriCalcBCFunction(self, func, *args, **kwargs):
        if not self.config['doCalcBC']:
            raise RuntimeError("setWriCalcBCFunction failed: -calc_bc not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fwricbc = context # this saves a pointer to context
            ierr = TMMSetWriCalcBCFunction(self.state, wriCalcBCCB, <void*>context)
        else:
            ierr = TMMSetWriCalcBCFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetWriCalcBCFunction failed with error code {ierr}")

    def setFinCalcBCFunction(self, func, *args, **kwargs):
        if not self.config['doCalcBC']:
            raise RuntimeError("setFinCalcBCFunction failed: -calc_bc not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.ffincbc = context # this saves a pointer to context
            ierr = TMMSetFinCalcBCFunction(self.state, finCalcBCCB, <void*>context)
        else:
            ierr = TMMSetFinCalcBCFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetFinCalcBCFunction failed with error code {ierr}")

    def setReiCalcBCFunction(self, func, *args, **kwargs):
        if not self.config['doCalcBC']:
            raise RuntimeError("setReiCalcBCFunction failed: -calc_bc not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.freicbc = context # this saves a pointer to context
            ierr = TMMSetReiCalcBCFunction(self.state, reiCalcBCCB, <void*>context)
        else:
            ierr = TMMSetReiCalcBCFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetReiCalcBCFunction failed with error code {ierr}")

    def setIniMonitorFunction(self, func, *args, **kwargs):
        if not self.config['useMonitor']:
            raise RuntimeError("setIniMonitorFunction failed: -use_monitor not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.finimon = context # this saves a pointer to context
            if debug: print(f"Ini monitor func is {context}")
            ierr = TMMSetIniMonitorFunction(self.state, iniMonitorCB, <void*>context)
        else:
            ierr = TMMSetIniMonitorFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetIniMonitorFunction failed with error code {ierr}")

    def setCalcMonitorFunction(self, func, *args, **kwargs):
        if not self.config['useMonitor']:
            raise RuntimeError("setCalcMonitorFunction failed: -use_monitor not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fcalcmon = context # this saves a pointer to context
            ierr = TMMSetCalcMonitorFunction(self.state, calcMonitorCB, <void*>context)
        else:
            ierr = TMMSetCalcMonitorFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetCalcMonitorFunction failed with error code {ierr}")

    def setWriMonitorFunction(self, func, *args, **kwargs):
        if not self.config['useMonitor']:
            raise RuntimeError("setWriMonitorFunction failed: -use_monitor not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fwrimon = context # this saves a pointer to context
            ierr = TMMSetWriMonitorFunction(self.state, wriMonitorCB, <void*>context)
        else:
            ierr = TMMSetWriMonitorFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetWriMonitorFunction failed with error code {ierr}")

    def setFinMonitorFunction(self, func, *args, **kwargs):
        if not self.config['useMonitor']:
            raise RuntimeError("setFinMonitorFunction failed: -use_monitor not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.ffinmon = context # this saves a pointer to context
            ierr = TMMSetFinMonitorFunction(self.state, finMonitorCB, <void*>context)
        else:
            ierr = TMMSetFinMonitorFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetFinMonitorFunction failed with error code {ierr}")

    def setIniMisfitFunction(self, func, *args, **kwargs):
        if not self.config['doMisfit']:
            raise RuntimeError("setIniMisfitFunction failed: -calc_misfit not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.finimis = context # this saves a pointer to context
            ierr = TMMSetIniMisfitFunction(self.state, iniMisfitCB, <void*>context)
        else:
            ierr = TMMSetIniMisfitFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetIniMisfitFunction failed with error code {ierr}")

    def setCalcMisfitFunction(self, func, *args, **kwargs):
        if not self.config['doMisfit']:
            raise RuntimeError("setCalcMisfitFunction failed: -calc_misfit not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fcalcmis = context # this saves a pointer to context
            ierr = TMMSetCalcMisfitFunction(self.state, calcMisfitCB, <void*>context)
        else:
            ierr = TMMSetCalcMisfitFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetCalcMisfitFunction failed with error code {ierr}")

    def setWriMisfitFunction(self, func, *args, **kwargs):
        if not self.config['doMisfit']:
            raise RuntimeError("setWriMisfitFunction failed: -calc_misfit not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.fwrimis = context # this saves a pointer to context
            ierr = TMMSetWriMisfitFunction(self.state, wriMisfitCB, <void*>context)
        else:
            ierr = TMMSetWriMisfitFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetWriMisfitFunction failed with error code {ierr}")

    def setFinMisfitFunction(self, func, *args, **kwargs):
        if not self.config['doMisfit']:
            raise RuntimeError("setFinMisfitFunction failed: -calc_misfit not specified")
        if func is not None:
            context = (func, args, kwargs, self)
            self.ffinmis = context # this saves a pointer to context
            ierr = TMMSetFinMisfitFunction(self.state, finMisfitCB, <void*>context)
        else:
            ierr = TMMSetFinMisfitFunction(self.state, NULL, NULL)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetFinMisfitFunction failed with error code {ierr}")

    def getConfig(self):
        if True: #toBool(self.state.isInitialized):
          config = _getStateConfig(self.state)
          return config
        else:
          raise RuntimeError("getConfig failed because state is uninitialized!")    
            
    def getNumTracers(self):
        return toInt(self.state.numTracers)

    def setFromOptions(self, prefix="", doOutput=True):
        cdef const char* pre
        cdef PetscScalar *data_pointer
        
        pre = PyBytes_AsString(prefix.encode())

        cdef PetscErrorCode ierr
        ierr = TMMSetFromOptions(self.state, pre, asBool(doOutput))
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMSetFromOptions failed with error code {ierr}")

        self.stateId=toInt(self.state.stateId)
        self.numTracers=toInt(self.state.numTracers)
        self.config=_getStateConfig(self.state)
        
        OptDBState=PETSc.Options(prefix)
        if OptDBState.hasName("tracer_names"):
          self.tracerNames=OptDBState.getString("tracer_names").split(",")
        else:
          self.tracerNames=None

        rank = petsc4py.PETSc.COMM_WORLD.getRank()
        c = [[] for itr in range(self.numTracers)]
        for itr in range(self.numTracers):
#           Get pointer to data
            ierr = VecGetArray(self.state.c[itr], &data_pointer)
            if ierr != 0: raise PETSc.Error(ierr)
            myId = f"Id={self.stateId}-c[{itr}], r={rank}"
            c[itr] = _ToNumpyRealArray(lSize, data_pointer, self.state.c[itr], myId)
        self.c=_ImmutableList(c)
        del c

        if self.config['useForcingFromFile']:
            qf = [[] for itr in range(self.numTracers)]
            for itr in range(self.numTracers):
#               Get pointer to data
                ierr = VecGetArray(self.state.qf[itr], &data_pointer)
                if ierr != 0: raise PETSc.Error(ierr)
                myId = f"Id={self.stateId}-qf[{itr}], r={rank}"
                qf[itr] = _ToNumpyRealArray(lSize, data_pointer, self.state.qf[itr], myId)
            self.qf=_ImmutableList(qf)
            del qf

        if self.config['useExternalForcing']:
            qef = [[] for itr in range(self.numTracers)]
            for itr in range(self.numTracers):
#               Get pointer to data
                ierr = VecGetArray(self.state.qef[itr], &data_pointer)
                if ierr != 0: raise PETSc.Error(ierr)
                myId = f"Id={self.stateId}-qef[{itr}], r={rank}"
                qef[itr] = _ToNumpyRealArray(lSize, data_pointer, self.state.qef[itr], myId)
            self.qef=_ImmutableList(qef)
            del qef

        if self.config['relaxTracer']:
            qrel = [[] for itr in range(self.numTracers)]
            for itr in range(self.numTracers):
#               Get pointer to data
                ierr = VecGetArray(self.state.qrel[itr], &data_pointer)
                if ierr != 0: raise PETSc.Error(ierr)
                myId = f"Id={self.stateId}-qrel[{itr}], r={rank}"
                self.qrel[itr] = _ToNumpyRealArray(lSize, data_pointer, self.state.qrel[itr], myId)
            self.qrel=_ImmutableList(qrel)
            del qrel

        if self.config['usePrescribedBC']:
            cbc = [[] for itr in range(self.numTracers)]
            for itr in range(self.numTracers):
#               Get pointer to data
                ierr = VecGetArray(self.state.cbc[itr],&data_pointer)
                if ierr != 0: raise PETSc.Error(ierr)
                myId = f"Id={self.stateId}-cbc[{itr}], r={rank}"
                cbc[itr] = _ToNumpyRealArray(lBCSize, data_pointer, self.state.cbc[itr], myId)
            self.cbc=_ImmutableList(cbc)
            del cbc

            cbf = [[] for itr in range(self.numTracers)]
            for itr in range(self.numTracers):
#               Get pointer to data
                ierr = VecGetArray(self.state.cbf[itr],&data_pointer)
                if ierr != 0: raise PETSc.Error(ierr)
                myId = f"Id={self.stateId}-cbf[{itr}], r={rank}"
                cbf[itr] = _ToNumpyRealArray(lBCSize, data_pointer, self.state.cbf[itr], myId)
            self.cbf=_ImmutableList(cbf)
            del cbf
        
    def forcingUpdate(self, tc, Iterc, iLoop):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscInt Itercc = asInt(Iterc)
        cdef PetscInt iLoopc = asInt(iLoop)
        cdef PetscErrorCode ierr
        ierr = TMMForcingUpdate(tcc, Itercc, iLoopc, self.state)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMForcingUpdate failed with error code {ierr}")

    def timeStep(self, tc, Iterc, iLoop):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscInt Itercc = asInt(Iterc)
        cdef PetscInt iLoopc = asInt(iLoop)
        cdef PetscErrorCode ierr
        ierr = TMMTimeStep(tcc, Itercc, iLoopc, self.state)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMTimeStep failed with error code {ierr}")

    def timeStepPost(self, tc, Iterc, iLoop):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscInt Itercc = asInt(Iterc)
        cdef PetscInt iLoopc = asInt(iLoop)
        cdef PetscErrorCode ierr
        # Careful; need to pass time at end of time step
        ierr = TMMTimeStepPost(tcc, Iterc, iLoop, self.state)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMTimeStepPost failed with error code {ierr}")

    def output(self, tc, Iterc, iLoop):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscInt Itercc = asInt(Iterc)
        cdef PetscInt iLoopc = asInt(iLoop)
        cdef PetscErrorCode ierr
        # Careful; need to pass time at end of time step
        ierr = TMMOutput(tcc, Itercc, iLoopc, self.state)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMOutput failed with error code {ierr}")

    def pushConfigVar(self,**kwargs):
        for key in kwargs:
          if key=='isInitializedExternalForcing':
            val=kwargs[key]
            if type(val) is not bool:
              raise RuntimeError(f"Incorrect type for {key}. Must be a bool!")
            self.state.isInitializedExternalForcing=asBool(val)
          elif key=='isInitializedCalcBC':
            val=kwargs[key]
            if type(val) is not bool:
              raise RuntimeError(f"Incorrect type for {key}. Must be a bool!")
            self.state.isInitializedCalcBC=asBool(val)
        # Update config
        self.config=_getStateConfig(self.state)
        
    def destroy(self, tc):
        cdef PetscErrorCode ierr

        if debug: print(f"Within destroy at time {tc}")

        if self.c:
            del self.c

        if self.qf:
            del self.qf

        if self.qef:
            del self.qef

        if self.qrel:
            del self.qrel

        if self.cbc:
            del self.cbc

        if self.cbf:
            del self.cbf

        ierr = TMMDestroy(asReal(tc), &self.state)
        if ierr != 0: raise PETSc.Error(ierr)
#         if ierr != 0:
#             raise RuntimeError(f"TMMDestroy failed with error code {ierr}")

    def setTracer(self, alpha, itr=0):
        cdef PetscScalar sval = asReal(alpha)
        cdef PetscErrorCode ierr
        ierr = VecSet(self.state.c[itr], sval)
        if ierr != 0: raise PETSc.Error(ierr)

    def loadTracerFromFile(self, typ, filename, itr, iRec=0):
#     Note: iRec is 0-based while VecLoadIntoVectorRandomAccessFromFile is 1-based
      cdef PetscErrorCode ierr
      cdef const char* c_filename
      cdef PetscInt iRecc = asInt(iRec)+1
      c_filename = PyBytes_AsString(filename.encode())
      
      if typ is 'tracer':
        nl = nb
        ierr = VecLoadIntoVectorRandomAccessFromFile(c_filename,self.state.c[itr], asInt(nl), iRecc)
        if ierr != 0: raise PETSc.Error(ierr)
      elif typ is 'cbc':
        nl = gBCSize
        ierr = VecLoadIntoVectorRandomAccessFromFile(c_filename,self.state.cbc[itr], asInt(nl), iRecc)
        if ierr != 0: raise PETSc.Error(ierr)
      elif typ is 'cbf':
        nl = gBCSize
        ierr = VecLoadIntoVectorRandomAccessFromFile(c_filename,self.state.cbf[itr], asInt(nl), iRecc)
        if ierr != 0: raise PETSc.Error(ierr)
      else:
        raise RuntimeError(f"Unknown type of vector {typ} when loading")

def _cleanUp():
    tconfig=getTimeConfig()
    tc = tconfig['time0'] + tconfig['deltaTClock']*tconfig['maxSteps']
    while len(_statescreated)>0:
      s=_statescreated.pop(0)
      if debug: print(f"Calling destroy at time {tc}")
      s.destroy(tc)
    if debug: print(f"Calling finalize at time {time0}, {maxSteps}, {tc}")
    finalize(tc)

atexit.register(_cleanUp)

def initialize():

    cdef PetscInt Iter0ret
    cdef PetscInt maxStepsret
    cdef PetscScalar time0ret
    cdef PetscScalar deltaTClockret
    cdef PetscErrorCode ierr

    # Call the C function
    ierr = TMMInitialize(&Iter0ret, &maxStepsret, &time0ret, &deltaTClockret)
    if ierr != 0: raise PETSc.Error(ierr)
#     if ierr != 0:
#         raise RuntimeError(f"TMMInitialize failed with error code {ierr}")

    # Return results as a Python tuple
    return toInt(Iter0ret), toInt(maxStepsret), toReal(time0ret), toReal(deltaTClockret)

def updateTMs(tc):
    cdef PetscErrorCode ierr
    cdef PetscScalar tcc = asReal(tc)
    ierr = TMMUpdateTMs(tcc)
    if ierr != 0: raise PETSc.Error(ierr)
#     if ierr != 0:
#         raise RuntimeError(f"TMMUpdateTMs failed with error code {ierr}")

def finalize(tc):
    cdef PetscErrorCode ierr
    cdef PetscScalar tcc = asReal(tc)
    ierr = TMMFinalize(tcc)
    if ierr != 0: raise PETSc.Error(ierr)
#     if ierr != 0:
#         raise RuntimeError(f"TMMFinalize failed with error code {ierr}")

def getProfileConfig():
    config = {key: False for key in profileVars}
    if toBool(useProfiles):
        config["useProfiles"] = toBool(useProfiles)
        for k in config.keys():
          if k == "lSize":
              config[k] = toInt(lSize)
          if k == "lNumProfiles":
              config[k] = toInt(lNumProfiles)
          if k == "totalNumProfiles":
              config[k] = toInt(totalNumProfiles)
          # We copy these as we don't want a user to accidentally change the values!
          if k == "lProfileLength":
              config[k] = copyNumpyIntArrayFromPointer(lProfileLength, lNumProfiles)
          if k == "lStartIndices":
              config[k] = copyNumpyIntArrayFromPointer(lStartIndices, lNumProfiles)
          if k == "lEndIndices":
              config[k] = copyNumpyIntArrayFromPointer(lEndIndices, lNumProfiles)
    else:
        config = None    
    return config

def getTimeConfig():
    config = {'Iter0': toInt(Iter0), 'time0': toReal(time0), 'maxSteps': toInt(maxSteps), 'deltaTClock': toReal(deltaTClock)}
    return config

def getSizes():
    return toInt(nb), toInt(lSize)

def getBCSizes():
    if toBool(prescribedBCInUse):
      return toInt(gBCSize), toInt(lBCSize)
    else:  
      raise RuntimeError("getBCSizes failed because no prescribed BCs are in use!")
      
cdef _getStateConfig(PetscTMMState state):
    if True: #toBool(state.isInitialized):
      config = {key: False for key in configVars}
      for k in config.keys():
        if k == "useExternalForcing": config[k] = toBool(state.useExternalForcing)
        if k == "useForcingFromFile": config[k] = toBool(state.useForcingFromFile)
        if k == "usePrescribedBC": config[k] = toBool(state.usePrescribedBC)
        if k == "applyExternalForcing": config[k] = toBool(state.applyExternalForcing)
        if k == "applyForcingFromFile": config[k] = toBool(state.applyForcingFromFile)
        if k == "applyBC": config[k] = toBool(state.applyBC)
        if k == "periodicForcing": config[k] = toBool(state.periodicForcing)
        if k == "timeDependentForcing": config[k] = toBool(state.timeDependentForcing)
        if k == "constantForcing": config[k] = toBool(state.constantForcing)
        if k == "periodicBC": config[k] = toBool(state.periodicBC)
        if k == "timeDependentBC": config[k] = toBool(state.timeDependentBC)
        if k == "constantBC": config[k] = toBool(state.constantBC)
        if k == "doCalcBC": config[k] = toBool(state.doCalcBC)
        if k == "useMonitor": config[k] = toBool(state.useMonitor)
        if k == "doMisfit": config[k] = toBool(state.doMisfit)
        if k == "relaxTracer": config[k] = toBool(state.relaxTracer)
        if k == "isInitializedExternalForcing": config[k] = toBool(state.isInitializedExternalForcing)
        if k == "isInitializedCalcBC": config[k] = toBool(state.isInitializedCalcBC)
        if k == "doOutput": config[k] = toBool(state.doOutput)
        if k == "appendOutput": config[k] = toBool(state.appendOutput)
        if k == "writePickup": config[k] = toBool(state.writePickup)
        if k == "doWriteBC": config[k] = toBool(state.doWriteBC)
        if k == "doWriteQF": config[k] = toBool(state.doWriteQF)
        if k == "doWriteQEF": config[k] = toBool(state.doWriteQEF)
        if k == "pickupFromFile": config[k] = toBool(state.pickupFromFile)
        if k == "doTimeAverage": config[k] = toBool(state.doTimeAverage)
        if k == "avgAppendOutput": config[k] = toBool(state.avgAppendOutput)
        if k == "doExtraWrite": config[k] = toBool(state.doExtraWrite)
        if k == "numTracers": config[k] = toInt(state.numTracers)
      return config
    else:
      raise RuntimeError("getConfig failed because state is uninitialized!")

cdef class StepTimer():

    def __cinit__(self):
        self.timer = NULL
      
    def __init__(self):
        pass
        
    def create(self, prefix="", prefix2prefix=None, PetscInt startTimeStep=0):
        cdef PetscErrorCode ierr 
        cdef const char* c_pre
        cdef const char* c_pre2pre
        cdef PetscStepTimer newtimer
        c_pre = PyBytes_AsString(prefix.encode())
        if prefix2prefix is None:
          c_pre2pre = NULL
        else:  
          c_pre2pre = PyBytes_AsString(prefix2prefix.encode())
        ierr = StepTimerCreate(&newtimer)
        if ierr != 0: raise PETSc.Error(ierr)
        self.timer=newtimer
        ierr = StepTimerIni(c_pre, c_pre2pre, startTimeStep, self.timer)
        if ierr != 0: raise PETSc.Error(ierr)
#SPK do I need a toInt here?
        self.count = self.timer.count
        self.startTimeStep = self.timer.startTimeStep
        self.numTimeSteps = self.timer.numTimeSteps

        return self

    def update(self, PetscInt Iter):
        cdef PetscErrorCode ierr 
        ierr = StepTimerUpdate(Iter, self.timer)
        if ierr != 0: raise PETSc.Error(ierr)
        self.count = self.timer.count
        self.startTimeStep = self.timer.startTimeStep
        self.numTimeSteps = self.timer.numTimeSteps

    def incr(self):
        self.timer.count = self.timer.count + 1
        self.count = self.timer.count

cdef class PeriodicTimer():

    def __cinit__(self):
        self.timer = NULL
      
    def __init__(self):
        pass
        
    def create(self, prefix="", prefix2prefix=None, np.ndarray[double, ndim=1, mode="c"] tdp=None, cyclePeriod=None):
        cdef PetscErrorCode ierr 
        cdef const char* c_pre
        cdef const char* c_pre2pre
        cdef PetscScalar *c_tdp
        cdef PetscPeriodicTimer newtimer

        OptDB = PETSc.Options() # This needs to be here; moving it below somehow corrupts c_pre

        if prefix2prefix is None:
          pre=prefix
        else:  
          pre=prefix2prefix+prefix

#       This needs to be here as well. OptDB.setValue somehow corrupts c_pre.
        if tdp is None:
          c_tdp = NULL
        else:
          try:
#             Check the size ...
              n=OptDB.getInt(pre+"num_per_period")
              if n != len(tdp):   # tdp.shape[0] gives a compilation error
                raise RuntimeError(f"Length of tdp array ({len(tdp)}) does not match size specified in options database ({n}) for {pre+'num_per_period'}")
          except KeyError:
              print(f"Warning: no options database value found for {pre+'num_per_period'}. Setting options database.")
              OptDB.setValue(pre+"num_per_period",len(tdp))
          c_tdp = &tdp[0]

        if cyclePeriod is not None:
            print(f"Setting options database for {pre+'cycle_period'}. This will override any previously set value.")
            OptDB.setValue(pre+"cycle_period",cyclePeriod)

        c_pre = PyBytes_AsString(prefix.encode())
        if prefix2prefix is None:
          c_pre2pre = NULL
        else:  
          c_pre2pre = PyBytes_AsString(prefix2prefix.encode())

        ierr = PeriodicTimerCreate(&newtimer)
        if ierr != 0: raise PETSc.Error(ierr)
        self.timer=newtimer
        ierr = PeriodicTimerIni(c_pre, c_pre2pre, c_tdp, self.timer)
        if ierr != 0: raise PETSc.Error(ierr)

        return self

cdef class TimeDependentTimer():

    def __cinit__(self):
        self.timer = NULL
      
    def __init__(self):
        pass
        
    def create(self, prefix="", prefix2prefix=None, np.ndarray[double, ndim=1, mode="c"] tdt=None):
        cdef PetscErrorCode ierr 
        cdef const char* c_pre
        cdef const char* c_pre2pre
        cdef PetscScalar *c_tdt
        cdef PetscTimeDependentTimer newtimer

        OptDB = PETSc.Options() # This needs to be here; moving it below somehow corrupts c_pre

        if prefix2prefix is None:
          pre=prefix
        else:  
          pre=prefix2prefix+prefix

#       This needs to be here as well. OptDB.setValue somehow corrupts c_pre.
        if tdt is None:
          c_tdt = NULL
        else:
          try:
#             Check the size ...
              n=OptDB.getInt(pre+"num_times")
              if n != len(tdt):   # tdt.shape[0] gives a compilation error
                raise RuntimeError(f"Length of tdt array ({len(tdt)}) does not match size specified in options database ({n}) for {pre+'num_times'}")
          except KeyError:
              print(f"Warning: no options database value found for {pre+'num_times'}. Setting options database.")
              OptDB.setValue(pre+"num_times",len(tdt))
          c_tdt = &tdt[0]

        c_pre = PyBytes_AsString(prefix.encode())
        if prefix2prefix is None:
          c_pre2pre = NULL
        else:  
          c_pre2pre = PyBytes_AsString(prefix2prefix.encode())

        ierr = TimeDependentTimerCreate(&newtimer)
        if ierr != 0: raise PETSc.Error(ierr)
        self.timer=newtimer
        ierr = TimeDependentTimerIni(c_pre, c_pre2pre, c_tdt, self.timer)
        if ierr != 0: raise PETSc.Error(ierr)

        return self

cdef class PeriodicVec():

    def __cinit__(self):
        self.pvec = NULL
      
    def __init__(self):
        pass
        
    def create(self, typ='tracer', np.ndarray[double, ndim=1, mode="c"] buf=None):
        cdef PetscErrorCode ierr 
        cdef PetscPeriodicVec newvec
        cdef PetscScalar *data_pointer

        ierr = PeriodicVecCreate(&newvec)
        if ierr != 0: raise PETSc.Error(ierr)
        self.pvec=newvec
        
        if typ is None or typ is 'tracer':
          self.lDim = lSize
        else:
          self.lDim = lBCSize

        if buf is None:
          self.frombuf = 0
          ierr = VecCreateFromLocalSize(self.lDim, &self.c)
          if ierr != 0: raise PETSc.Error(ierr)
          ierr = VecGetArray(self.c, &data_pointer)
          if ierr != 0: raise PETSc.Error(ierr)
#         Note: the memory is owned by self.c          
        else:
          if len(buf) != self.lDim:
            raise RuntimeError(f"Provided array length does not match expected length {self.lDim}")
          self.frombuf = 1
          data_pointer = &buf[0]
          ierr = VecCreateWithArrayFromLocalSize(self.lDim, data_pointer, &self.c)
          if ierr != 0: raise PETSc.Error(ierr)
#         Note: the memory is owned by the array buf (which may be owned by a vec)

        self.arr = _ToNumpyRealArray(self.lDim, data_pointer, self.c, "")

        return self

    def interp(self, tc, PeriodicTimer timer, filename):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscErrorCode ierr
        cdef const char* c_filename
        
        c_filename = PyBytes_AsString(filename.encode())
        ierr = PeriodicVecInterp(tcc, &self.c, timer.timer.cyclePeriod, timer.timer.numPerPeriod, timer.timer.tdp, self.pvec, c_filename)
        if ierr != 0: raise PETSc.Error(ierr)

    def destroy(self):
        cdef PetscErrorCode ierr
        
        del self.arr # calls VecRestoreArray (??????)
        ierr = VecDestroy(&self.c)
        if ierr != 0: raise PETSc.Error(ierr)
        ierr = PeriodicVecDestroy(&self.pvec)
        if ierr != 0: raise PETSc.Error(ierr)

cdef class TimeDependentVec():

    def __cinit__(self):
        self.pvec = NULL
      
    def __init__(self):
        pass
        
    def create(self, typ='tracer', np.ndarray[double, ndim=1, mode="c"] buf=None):
        cdef PetscErrorCode ierr 
        cdef PetscTimeDependentVec newvec
        cdef PetscScalar *data_pointer

        ierr = TimeDependentVecCreate(&newvec)
        if ierr != 0: raise PETSc.Error(ierr)
        self.pvec=newvec
        
        if typ is None or typ is 'tracer':
          self.lDim = lSize
        else:
          self.lDim = lBCSize

        if buf is None:
          self.frombuf = 0
          ierr = VecCreateFromLocalSize(self.lDim, &self.c)
          if ierr != 0: raise PETSc.Error(ierr)
          ierr = VecGetArray(self.c, &data_pointer)
          if ierr != 0: raise PETSc.Error(ierr)
#         Note: the memory is owned by self.c          
        else:
          if len(buf) != self.lDim:
            raise RuntimeError(f"Provided array length does not match expected length {self.lDim}")
          self.frombuf = 1
          data_pointer = &buf[0]
          ierr = VecCreateWithArrayFromLocalSize(self.lDim, data_pointer, &self.c)
          if ierr != 0: raise PETSc.Error(ierr)
#         Note: the memory is owned by the array buf (which may be owned by a vec)

        self.arr = _ToNumpyRealArray(self.lDim, data_pointer, self.c, "")

        return self

    def interp(self, tc, TimeDependentTimer timer, filename):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscErrorCode ierr
        cdef const char* c_filename

        c_filename = PyBytes_AsString(filename.encode())
        ierr = TimeDependentVecInterp(tcc, &self.c, timer.timer.numTimes, timer.timer.tdt, self.pvec, c_filename)
        if ierr != 0: raise PETSc.Error(ierr)

    def destroy(self):
        cdef PetscErrorCode ierr
        
        del self.arr # calls VecRestoreArray (??????)
        ierr = VecDestroy(&self.c)
        if ierr != 0: raise PETSc.Error(ierr)
        ierr = TimeDependentVecDestroy(&self.pvec)
        if ierr != 0: raise PETSc.Error(ierr)

cdef class PeriodicArray():

    def __cinit__(self):
        self.parr = NULL
      
    def __init__(self):
        pass
        
    def create(self, arrayLength, np.ndarray[double, ndim=1, mode="c"] buf=None):
        cdef PetscErrorCode ierr 
        cdef PetscPeriodicArray newarr
#         cdef PetscScalar *data_pointer

        ierr = PeriodicArrayCreate(&newarr, arrayLength)
        if ierr != 0: raise PETSc.Error(ierr)
        self.parr=newarr
        self.lDim = arrayLength

        if buf is None:
          self.frombuf = 0
          self.data_pointer = <PetscScalar *> PyMem_Malloc(self.lDim*sizeof(PetscScalar))
          self.arr = _ToNumpyRealArrayWithDel(self.lDim, self.data_pointer, NULL, "")
#         Note: the memory is owned by self.arr and will be freed when the destroy method is called
        else:
          if len(buf) != self.lDim:
            raise RuntimeError(f"Provided array length does not match expected length {self.lDim}")
          self.frombuf = 1
          self.data_pointer = &buf[0]
          self.arr = _ToNumpyRealArray(self.lDim, self.data_pointer, NULL, "")
#         Note: the memory is owned by the array buf

        return self

    def interp(self, tc, PeriodicTimer timer, filename):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscErrorCode ierr
        cdef const char* c_filename
        
        c_filename = PyBytes_AsString(filename.encode())
        ierr = interpPeriodicProfileSurfaceScalarData(tcc, self.data_pointer, timer.timer.cyclePeriod, timer.timer.numPerPeriod, timer.timer.tdp, self.parr, c_filename)
        if ierr != 0: raise PETSc.Error(ierr)
                                    
    def destroy(self):
        cdef PetscErrorCode ierr
        
        del self.arr
        ierr = PeriodicArrayDestroy(&self.parr)
        if ierr != 0: raise PETSc.Error(ierr)

cdef class TimeDependentArray():

    def __cinit__(self):
        self.parr = NULL
      
    def __init__(self):
        pass
        
    def create(self, arrayLength, np.ndarray[double, ndim=1, mode="c"] buf=None):
        cdef PetscErrorCode ierr 
        cdef PetscTimeDependentArray newarr
#         cdef PetscScalar *data_pointer

        ierr = TimeDependentArrayCreate(&newarr, arrayLength)
        if ierr != 0: raise PETSc.Error(ierr)
        self.parr=newarr
        self.lDim = arrayLength

        if buf is None:
          self.frombuf = 0
          self.data_pointer = <PetscScalar *> PyMem_Malloc(self.lDim*sizeof(PetscScalar))
          self.arr = _ToNumpyRealArrayWithDel(self.lDim, self.data_pointer, NULL, "")
#         Note: the memory is owned by self.arr and will be freed when the destroy method is called
        else:
          if len(buf) != self.lDim:
            raise RuntimeError(f"Provided array length does not match expected length {self.lDim}")
          self.frombuf = 1
          self.data_pointer = &buf[0]
          self.arr = _ToNumpyRealArray(self.lDim, self.data_pointer, NULL, "")
#         Note: the memory is owned by the array buf

        return self

    def interp(self, tc, TimeDependentTimer timer, filename):
        cdef PetscScalar tcc = asReal(tc)
        cdef PetscErrorCode ierr
        cdef const char* c_filename
        
        c_filename = PyBytes_AsString(filename.encode())
        ierr = interpTimeDependentProfileSurfaceScalarData(tcc, self.data_pointer, timer.timer.numTimes, timer.timer.tdt, self.parr, c_filename)
        if ierr != 0: raise PETSc.Error(ierr)
                                    
    def destroy(self):
        cdef PetscErrorCode ierr
        
        del self.arr
        ierr = TimeDependentArrayDestroy(&self.parr)
        if ierr != 0: raise PETSc.Error(ierr)


def loadVecIntoArray(filename, np.ndarray[double, ndim=1, mode="c"] buf=None):
    cdef PetscErrorCode ierr
    cdef PetscScalar *data_pointer
    cdef const char* c_filename
#  I don't think this is optimal. See https://cython.readthedocs.io/en/latest/src/tutorial/memory_allocation.html and 
#  https://cython.readthedocs.io/en/latest/src/tutorial/array.html#array-array

    c_filename = PyBytes_AsString(filename.encode())
    
    if buf is None:
        data_pointer = <PetscScalar *> PyMem_Malloc(lSize*sizeof(PetscScalar))
        ierr = VecLoadIntoArray(lSize, c_filename, data_pointer);
        if ierr != 0: raise PETSc.Error(ierr)
        nparr = _ToNumpyRealArrayWithDel(lSize, data_pointer, NULL, "")
#       Note: the memory will be freed when nparr is deleted
        return nparr
    else:
        data_pointer = &buf[0]
        ierr = VecLoadIntoArray(lSize, c_filename, data_pointer);
        if ierr != 0: raise PETSc.Error(ierr)
#       Note: the memory is owned by the array buf

def readProfileScalarData(filename, numValsPerProfile=1, np.ndarray[double, ndim=1, mode="c"] buf=None):
    cdef PetscErrorCode ierr
    cdef PetscScalar *data_pointer
    cdef const char* c_filename

    c_filename = PyBytes_AsString(filename.encode())
    
    lDim = lNumProfiles*numValsPerProfile
    if buf is None:
        data_pointer = <PetscScalar *> PyMem_Malloc(lDim*sizeof(PetscScalar))
        ierr = readProfileSurfaceScalarData(c_filename,data_pointer,numValsPerProfile);
        if ierr != 0: raise PETSc.Error(ierr)
        nparr = _ToNumpyRealArrayWithDel(lDim, data_pointer, NULL, "")
#       Note: the memory will be freed when nparr is deleted
        return nparr
    else:
        if len(buf) != lDim:
          raise RuntimeError(f"Provided array length does not match expected length {lDim}")
        data_pointer = &buf[0]
        ierr = readProfileSurfaceScalarData(c_filename,data_pointer,numValsPerProfile);
        if ierr != 0: raise PETSc.Error(ierr)
#       Note: the memory is owned by the array buf
  
def writeProfileScalarData(filename, np.ndarray[double, ndim=1, mode="c"] arr not None, numValsPerProfile=1, appendToFile=False):
    cdef PetscErrorCode ierr
    cdef PetscScalar *data_pointer
    cdef const char* c_filename

    c_filename = PyBytes_AsString(filename.encode())
    
    lDim = lNumProfiles*numValsPerProfile
    if len(arr) != lDim:
      raise RuntimeError(f"arr length does not match expected length {lDim}")
    data_pointer = &arr[0]
    ierr = writeProfileSurfaceScalarData(c_filename, data_pointer, asInt(numValsPerProfile), asBool(appendToFile))
    if ierr != 0: raise PETSc.Error(ierr)

def writeArrayToVec(filename, np.ndarray[double, ndim=1, mode="c"] arr not None, mode):
#         if   mode == 'r'  : return PETSC_FILE_MODE_READ
#         elif mode == 'w'  : return PETSC_FILE_MODE_WRITE
#         elif mode == 'a'  : return PETSC_FILE_MODE_APPEND
#         elif mode == 'r+' : return PETSC_FILE_MODE_UPDATE
#         elif mode == 'w+' : return PETSC_FILE_MODE_UPDATE
#         elif mode == 'a+' : return PETSC_FILE_MODE_APPEND_UPDATE
#         elif mode == 'u'  : return PETSC_FILE_MODE_UPDATE
#         elif mode == 'au' : return PETSC_FILE_MODE_APPEND_UPDATE
#         elif mode == 'ua' : return PETSC_FILE_MODE_APPEND_UPDATE

    c = petsc4py.PETSc.Vec() #.create(comm=petsc4py.PETSc.COMM_WORLD)
    c.createWithArray(arr,comm=petsc4py.PETSc.COMM_WORLD)
    
#     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,mode,&fd);CHKERRQ(ierr);
    viewer = PETSc.Viewer().createBinary(filename, mode=mode)
    c.view(viewer)
    viewer.destroy()
    c.destroy()

def writeBinaryArray(filename, np.ndarray[double, ndim=1, mode="c"] arr not None, appendToFile=False):
#   Note: writeBinaryScalarData only writes on rank 0. You can call this from all ranks
    cdef PetscErrorCode ierr
    cdef PetscScalar *data_pointer
    cdef const char* c_filename

    c_filename = PyBytes_AsString(filename.encode())
    N=len(arr)
    data_pointer = &arr[0]
    ierr = writeBinaryScalarData(c_filename, data_pointer, N, asBool(appendToFile))
    if ierr != 0: raise PETSc.Error(ierr)

def barrier():

    cdef PetscErrorCode ierr
    ierr = Barrier()
    if ierr != 0: raise PETSc.Error(ierr)

def globalScalarSum(x):

    cdef PetscErrorCode ierr
    cdef PetscScalar xc = asReal(x)
    cdef PetscScalar tot

    ierr = sumScalar(xc, &tot)
    if ierr != 0: raise PETSc.Error(ierr)

    return toReal(tot)

# This should be deprecated
# def globalSum(np.ndarray[double, ndim=1, mode="c"] xarr not None):
# 
#     x = petsc4py.PETSc.Vec()
#     x.createWithArray(xarr,comm=petsc4py.PETSc.COMM_WORLD)
# 
#     tot=x.sum()    
# 
#     return tot

def globalSum(np.ndarray[double, ndim=1, mode="c"] xarr not None):

    cdef PetscErrorCode ierr
    cdef PetscScalar *xdata_pointer
    cdef PetscScalar res

    nx = len(xarr)
    xdata_pointer = &xarr[0]
        
    ierr = sumArray(xdata_pointer, asInt(nx), &res)
    if ierr != 0: raise PETSc.Error(ierr)

    return toReal(res)

def dotProd(np.ndarray[double, ndim=1, mode="c"] xarr not None, np.ndarray[double, ndim=1, mode="c"] yarr not None):
    cdef PetscErrorCode ierr
    cdef PetscScalar *xdata_pointer
    cdef PetscScalar *ydata_pointer
    cdef PetscScalar res

    nx = len(xarr)
    ny = len(yarr)
    if nx != ny:
      raise RuntimeError(f"xarr length {nx} does not equal yarr length {ny}")

    xdata_pointer = &xarr[0]
    ydata_pointer = &yarr[0]
        
    ierr = dotProdArrays(xdata_pointer, ydata_pointer, asInt(nx), &res)
    if ierr != 0: raise PETSc.Error(ierr)

    return toReal(res)

# This should be deprecated
def dotProdProfileScalarData(np.ndarray[double, ndim=1, mode="c"] xarr not None, np.ndarray[double, ndim=1, mode="c"] yarr not None):
    cdef PetscErrorCode ierr
    cdef PetscScalar *xdata_pointer
    cdef PetscScalar *ydata_pointer
    cdef PetscScalar res

    if (len(xarr) != lNumProfiles):
      raise RuntimeError(f"xarr length does not match expected length {lNumProfiles}")

    if (len(yarr) != lNumProfiles):
      raise RuntimeError(f"yarr length does not match expected length {lNumProfiles}")

    xdata_pointer = &xarr[0]
    ydata_pointer = &yarr[0]
        
    ierr = dotProdProfileSurfaceScalarData(xdata_pointer, ydata_pointer, &res)
    if ierr != 0: raise PETSc.Error(ierr)

    return toReal(res)

# This should be deprecated
def dotProdVecs(np.ndarray[double, ndim=1, mode="c"] xarr not None, np.ndarray[double, ndim=1, mode="c"] yarr not None):

    x = petsc4py.PETSc.Vec()
    x.createWithArray(xarr,comm=petsc4py.PETSc.COMM_WORLD)
    y = petsc4py.PETSc.Vec()
    y.createWithArray(yarr,comm=petsc4py.PETSc.COMM_WORLD)

    res=x.dot(y)    

    x.destroy()
    y.destroy()

    return res

# These are public functions for use from other extension modules via cimport. Note that 
# no error checking is done for the validity of the pointer/length passed in.
cdef getNumpyRealArrayFromPointer(PetscScalar *data_pointer, PetscInt n):
  c=_ToNumpyRealArray(n, data_pointer, NULL, "")
  return c
 
cdef getNumpyIntArrayFromPointer(PetscInt *data_pointer, PetscInt n):
  c=_ToNumpyIntArray(n, data_pointer)
  return c

cdef copyNumpyRealArrayFromPointer(PetscScalar *data_pointer, PetscInt n):
  cdef PetscScalar *arr_pointer
  cdef int i;
  arr_pointer = <PetscScalar *> PyMem_Malloc(n*sizeof(PetscScalar))
  for i in range(n): arr_pointer[i]=data_pointer[i]
  arr = _ToNumpyRealArrayWithDel(n, arr_pointer, NULL, "")
  return arr

cdef copyNumpyIntArrayFromPointer(PetscInt *data_pointer, PetscInt n):
  cdef PetscInt *arr_pointer
  cdef int i;
  arr_pointer = <PetscInt *> PyMem_Malloc(n*sizeof(PetscInt))
  for i in range(n): arr_pointer[i]=data_pointer[i]
  arr = _ToNumpyIntArrayWithDel(n, arr_pointer)
  return arr

cdef class _RealArrayWrapper:
    cdef PetscScalar* data_ptr
    cdef int size
    cdef PetscVec pv
    cdef str myId
    cdef set_data(self, PetscInt size, PetscScalar* data_ptr, PetscVec pv, str myId):
        # Note: pass NULL as the argument for pv if data_ptr doesn't point to a PETSc Vec    
        self.data_ptr = data_ptr
        self.size = size
        self.pv = pv # we need to store a reference to this for the destructor
        self.myId = myId

    def __array__(self, dtype=None, copy=None):
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1-d array of length size
        ndarray = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        # This is called by Python when the array goes out of scope or all the 
        # references to it are gone which (I think) is what I want
        # Set debug to True to see this
        cdef PetscErrorCode ierr

        msg = f"Before GC {self.myId}"
        if debug: print(msg)
        if self.pv:
            msg = f"Doing GC {self.myId}"
            if debug: print(msg)
            ierr = VecRestoreArray(self.pv,&self.data_ptr)
            if ierr != 0: raise PETSc.Error(ierr)
# uncomment this for destructor but I don't think we want this as it will free the 
# memory and we want PETSc to manage the memory
#     def __dealloc__(self):
#          # This is called by Python when the array goes out of scope or all the 
#          # references to it are gone which (I think) is what I want
#         free(<void*>self.data_ptr)

cdef class _RealArrayWrapperWithDel:
    cdef PetscScalar* data_ptr
    cdef int size
    cdef PetscVec pv
    cdef str myId
    cdef set_data(self, PetscInt size, PetscScalar* data_ptr, PetscVec pv, str myId):
        # Note: pass NULL as the argument for pv if data_ptr doesn't point to a PETSc Vec
        self.data_ptr = data_ptr
        self.size = size
        self.pv = pv # we need to store a reference to this for the destructor
        self.myId = myId

    def __array__(self, dtype=None, copy=None):
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1-d array of length 'size
        ndarray = np.PyArray_SimpleNewFromData(1, shape, np.NPY_DOUBLE, self.data_ptr)
        return ndarray
    def __dealloc__(self):
        # This is called by Python when the array goes out of scope or all the 
        # references to it are gone which (I think) is what I want
        if debug: print("freeing data")
        PyMem_Free(self.data_ptr)

cdef class _IntArrayWrapper:
    cdef PetscInt* data_ptr
    cdef int size
    cdef set_data(self, int size, PetscInt* data_ptr):
        self.data_ptr = data_ptr
        self.size = size # we need to store a reference to this for the destructor

    def __array__(self, dtype=None, copy=None):
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, self.data_ptr)
        return ndarray

cdef class _IntArrayWrapperWithDel:
    cdef PetscInt* data_ptr
    cdef int size
    cdef set_data(self, int size, PetscInt* data_ptr):
        self.data_ptr = data_ptr
        self.size = size # we need to store a reference to this for the destructor

    def __array__(self, dtype=None, copy=None):
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, self.data_ptr)
        return ndarray
    def __dealloc__(self):
        # This is called by Python when the array goes out of scope or all the 
        # references to it are gone which (I think) is what I want
        if debug: print("freeing data")
        PyMem_Free(self.data_ptr)

# uncomment this for destructor but I don't think we want this as it will free the 
# memory and we want PETSc to manage the memory
#     def __dealloc__(self):
#         free(<void*>self.data_ptr)

cdef _ToNumpyRealArray(PetscInt size, PetscScalar *data_pointer, PetscVec pv, myId):
# Note: pass NULL as the argument for pv if data_pointer doesn't point to a PETSc Vec
#       pass "" for myId if not required
    cdef PetscErrorCode ierr
    cdef np.ndarray nparr
    varray = _RealArrayWrapper()
    varray.set_data(size, data_pointer,pv,myId)
    nparr = np.array(varray, copy=False)
    np.PyArray_SetBaseObject(nparr,varray)
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(varray)
#   We restore data_pointer here because it was created with a VecGetArray in the calling 
#   function in order to convert to a numpy array and is (should?) no longer needed in the 
#   calling function. Is that what I want?
    if pv:
        ierr = VecRestoreArray(pv, &data_pointer)
        if ierr != 0: raise PETSc.Error(ierr)
    return nparr    

cdef _ToNumpyIntArray(PetscInt size, PetscInt *data_pointer):
    cdef PetscErrorCode ierr
    cdef np.ndarray nparr
    varray = _IntArrayWrapper()
    varray.set_data(size, data_pointer)
    nparr = np.array(varray, copy=False)
    np.PyArray_SetBaseObject(nparr,varray)
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(varray)
    return nparr

cdef _ToNumpyRealArrayWithDel(PetscInt size, PetscScalar *data_pointer, PetscVec pv, myId):
# Note: pass NULL as the argument for pv if data_pointer doesn't point to a vector
    cdef PetscErrorCode ierr
    cdef np.ndarray nparr
    varray = _RealArrayWrapperWithDel()
    varray.set_data(size, data_pointer,pv,myId)
    nparr = np.array(varray, copy=False)
    np.PyArray_SetBaseObject(nparr,varray)
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(varray)
    if pv:
        ierr = VecRestoreArray(pv, &data_pointer)
        if ierr != 0: raise PETSc.Error(ierr)
    return nparr

cdef _ToNumpyIntArrayWithDel(PetscInt size, PetscInt *data_pointer):
    cdef PetscErrorCode ierr
    cdef np.ndarray nparr
    varray = _IntArrayWrapperWithDel()
    varray.set_data(size, data_pointer)
    nparr = np.array(varray, copy=False)
    np.PyArray_SetBaseObject(nparr,varray)
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(varray)
    return nparr

# 
# # Memoryview on a NumPy array
# narr = np.arange(27, dtype=np.dtype("i")).reshape((3, 3, 3))
# cdef int [:, :, :] narr_view = narr
# 
# # Memoryview on a C array
# cdef int[3][3][3] carr
# cdef int [:, :, :] carr_view = carr
# 
# # Memoryview on a Cython array
# cyarr = cvarray(shape=(3, 3, 3), itemsize=sizeof(int), format="i")
# cdef int [:, :, :] cyarr_view = cyarr
