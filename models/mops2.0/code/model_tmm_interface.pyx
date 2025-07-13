# Note: the asReal/asInt conversion below is not really necessary (as the C and python data types are 
# generally the same and the inline functions don't do any type conversion and simply return the same 
# value). Nor is it strictly correct since what we want is asScalar. But I couldn't get it to work.

from mops_forcing cimport ExternalForcingContext

from petsc4py.PETSc import Error

from tmm4py cimport PetscTMMState, TMMState, PetscErrorCode, PetscInt, PetscScalar, PetscReal, PetscBool
from tmm4py cimport getNumpyRealArrayFromPointer
from tmm4py cimport asBool, asInt, asReal, toBool, toInt, toReal

cdef extern from "tmm_share.h":
    PetscInt nb
    PetscInt lBCSize, gBCSize
    PetscInt lNumProfiles, lSize, numPrevProfiles, totalNumProfiles
    PetscBool useProfiles

cdef extern PetscErrorCode getExternalForcingContext(PetscTMMState state, void **ctx);

cdef extern from "tmm_external_forcing.h":
    PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);
    PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);
    PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscTMMState state, void *ctx);

cdef extern from "tmm_monitor.h":
    PetscErrorCode iniMonitor(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);
    PetscErrorCode calcMonitor(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode writeMonitor(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode finalizeMonitor(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);

cdef extern from "tmm_misfit.h":
    PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);
    PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, PetscTMMState state, void *ctx);
    PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, PetscTMMState state, void *ctx);

# Python gateway functions to C functions
def iniExternalForcingFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = iniExternalForcing(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def calcExternalForcingFn(tc, Iter, iLoop, TMMState State):
    cdef int ierr;
    ierr = calcExternalForcing(asReal(tc), asInt(Iter), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def writeExternalForcingFn(tc, Iter, iLoop, TMMState State):
    cdef int ierr;
    ierr = writeExternalForcing(asReal(tc), asInt(Iter), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def finalizeExternalForcingFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = finalizeExternalForcing(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def reInitializeExternalForcingFn(tc, Iter, iLoop, TMMState State):
    cdef int ierr;
    ierr = reInitializeExternalForcing(asReal(tc), asInt(Iter), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

#
def iniMonitorFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = iniMonitor(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def calcMonitorFn(tc, iLoop, TMMState State):
    cdef int ierr;
    ierr = calcMonitor(asReal(tc), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def writeMonitorFn(tc, iLoop, TMMState State):
    cdef int ierr;
    ierr = writeMonitor(asReal(tc), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def finalizeMonitorFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = finalizeMonitor(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

#
def iniMisfitFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = iniMisfit(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def calcMisfitFn(tc, iLoop, TMMState State):
    cdef int ierr;
    ierr = calcMisfit(asReal(tc), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def writeMisfitFn(tc, iLoop, TMMState State):
    cdef int ierr;
    ierr = writeMisfit(asReal(tc), asInt(iLoop), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def finalizeMisfitFn(tc, Iter, TMMState State):
    cdef int ierr;
    ierr = finalizeMisfit(asReal(tc), asInt(Iter), State.state, NULL);
    if ierr != 0: raise Error(ierr)

def pushToExternalForcingFn(TMMState State, **kwargs):
    cdef int ierr;
    cdef void *ctx = NULL;
    ierr = getExternalForcingContext(State.state, &ctx);
    if ierr != 0: raise Error(ierr)
    cdef ExternalForcingContext ef=<ExternalForcingContext>ctx;
#     print(f"Python: {toReal(ef.GRunoff)}")
#     print(f"Python: {toReal(ef.pCO2atm)}")
    for key in kwargs:
      if key == 'GRunoff':
          ef.GRunoff=asReal(kwargs[key])

def pullFromExternalForcingFn(TMMState State, *args):
    cdef int ierr;
    cdef void *ctx = NULL;
    ierr = getExternalForcingContext(State.state, &ctx);
    if ierr != 0: raise Error(ierr)
    cdef ExternalForcingContext ef=<ExternalForcingContext>ctx;
    n=len(args)
    out=[]
    for i in range(len(args)):
      if args[i] == 'GRunoff':
        out.append(toReal(ef.GRunoff))
      elif args[i] == 'ph':
        if toBool(useProfiles):
          v = getNumpyRealArrayFromPointer(ef.localph, lNumProfiles)
          out.append(v)
        else:
          raise RuntimeError("To pull localph requires useProfiles to be set to True")
    return out
 