#include "tmmimpl.h"

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMExtIniFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMExtCalcFunctionFn)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMExtWriFunctionFn)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMExtFinFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMExtReiFunctionFn)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMCalcBCIniFunctionFn)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMCalcBCCalcFunctionFn)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMCalcBCWriFunctionFn)(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMCalcBCFinFunctionFn)(PetscScalar tc, PetscInt Iterc, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMCalcBCReiFunctionFn)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMonitorIniFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMonitorCalcFunctionFn)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMonitorWriFunctionFn)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMonitorFinFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMisfitIniFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMisfitCalcFunctionFn)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMisfitWriFunctionFn)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode(TMMMisfitFinFunctionFn)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);

/* Public interfaces */
extern PetscErrorCode TMMCreate(TMMState *state);
extern PetscErrorCode TMMSetFromOptions(TMMState state, const char pre[], PetscBool doOutput);
extern PetscErrorCode TMMDestroy(PetscScalar tc, TMMState *state);
extern PetscErrorCode TMMSetOptionsPrefix(TMMState, const char[]);
extern PetscErrorCode TMMGetOptionsPrefix(TMMState, const char *[]);

// extern PetscErrorCode TMMSetExternalForcingFunctions(TMMState state, 
//   TMMExtIniFunctionFn *fini,
//   TMMExtCalcFunctionFn *fcalc,
//   TMMExtWriFunctionFn *fwri,
//   TMMExtFinFunctionFn *ffin,
//   TMMExtReiFunctionFn *frei
//   );

extern PetscErrorCode TMMSetIniExtForcFunction(TMMState state, TMMExtIniFunctionFn *fini, void *ctx);
extern PetscErrorCode TMMSetCalcExtForcFunction(TMMState state, TMMExtCalcFunctionFn *fcalc, void *ctx);
extern PetscErrorCode TMMSetWriExtForcFunction(TMMState state, TMMExtWriFunctionFn *fwri, void *ctx);
extern PetscErrorCode TMMSetFinExtForcFunction(TMMState state, TMMExtFinFunctionFn *ffin, void *ctx);
extern PetscErrorCode TMMSetReiExtForcFunction(TMMState state, TMMExtReiFunctionFn *frei, void *ctx);

extern PetscErrorCode TMMSetIniCalcBCFunction(TMMState state, TMMCalcBCIniFunctionFn *fini, void *ctx);
extern PetscErrorCode TMMSetCalcCalcBCFunction(TMMState state, TMMCalcBCCalcFunctionFn *fcalc, void *ctx);
extern PetscErrorCode TMMSetWriCalcBCFunction(TMMState state, TMMCalcBCWriFunctionFn *fwri, void *ctx);
extern PetscErrorCode TMMSetFinCalcBCFunction(TMMState state, TMMCalcBCFinFunctionFn *ffin, void *ctx);
extern PetscErrorCode TMMSetReiCalcBCFunction(TMMState state, TMMCalcBCReiFunctionFn *frei, void *ctx);

extern PetscErrorCode TMMSetIniMonitorFunction(TMMState state, TMMMonitorIniFunctionFn *fini, void *ctx);
extern PetscErrorCode TMMSetCalcMonitorFunction(TMMState state, TMMMonitorCalcFunctionFn *fcalc, void *ctx);
extern PetscErrorCode TMMSetWriMonitorFunction(TMMState state, TMMMonitorWriFunctionFn *fwri, void *ctx);
extern PetscErrorCode TMMSetFinMonitorFunction(TMMState state, TMMMonitorFinFunctionFn *ffin, void *ctx);

extern PetscErrorCode TMMSetIniMisfitFunction(TMMState state, TMMMisfitIniFunctionFn *fini, void *ctx);
extern PetscErrorCode TMMSetCalcMisfitFunction(TMMState state, TMMMisfitCalcFunctionFn *fcalc, void *ctx);
extern PetscErrorCode TMMSetWriMisfitFunction(TMMState state, TMMMisfitWriFunctionFn *fwri, void *ctx);
extern PetscErrorCode TMMSetFinMisfitFunction(TMMState state, TMMMisfitFinFunctionFn *ffin, void *ctx);

extern PetscErrorCode TMMComputeExtForcFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc);

// extern PetscErrorCode TMMSetCalcBCFunctions(TMMState state, 
//   TMMCalcBCIniFunctionFn *fini,
//   TMMCalcBCCalcFunctionFn *fcalc,
//   TMMCalcBCWriFunctionFn *fwri,
//   TMMCalcBCFinFunctionFn *ffin,
//   TMMCalcBCReiFunctionFn *frei,
//   void *ctx);

extern PetscErrorCode TMMComputeCalcBCFunction(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc);

extern PetscErrorCode TMMComputeMonitorFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc);

extern PetscErrorCode TMMComputeMisfitFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc);

// extern void* TMMGetFunction(TMMState state, TMMFUNCTYPE theFunc);
// extern TMMExtIniFunctionFn* TMMGetIniExternalForcingFunction(TMMState state);
// extern TMMExtCalcFunctionFn* TMMGetCalcExternalForcingFunction(TMMState state);
// extern TMMExtWriFunctionFn* TMMGetWriExternalForcingFunction(TMMState state);
// extern TMMExtFinFunctionFn* TMMGetFinExternalForcingFunction(TMMState state);
// extern TMMExtReiFunctionFn* TMMGetReiExternalForcingFunction(TMMState state);

// PETSC_EXTERN_TYPEDEF typedef void*(TMMExtCtxFunctionFn)(MPI_Comm comm, TMMState state);
// extern void* TMMCreateExternalForcingContext(MPI_Comm comm, TMMState state, TMMExtCtxFunctionFn *fctx);

// extern void* createExternalForcingContext(MPI_Comm comm, TMMState state);
// extern void* createCalcBCContext(MPI_Comm comm, TMMState state);

extern PetscErrorCode TMMUpdateTMs(PetscScalar tc);
extern PetscErrorCode TMMInitialize(PetscInt *Iter0, PetscInt *maxSteps, PetscScalar *time0, PetscScalar *deltaTClock);
extern PetscErrorCode TMMFinalize(PetscScalar tc);
extern PetscErrorCode TMMOutput(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state);
extern PetscErrorCode TMMForcingUpdate(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state);
extern PetscErrorCode TMMForcingReinitialize(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state);
extern PetscErrorCode TMMTimeStep(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state);
extern PetscErrorCode TMMTimeStepPost(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state);

extern PetscErrorCode doJacobianInitialize(Mat Q, TMMState state);
extern PetscErrorCode doJacobianCalc(PetscScalar tc, PetscScalar tf,PetscInt Iterc, PetscInt iLoop, Mat Q, TMMState state);
extern PetscErrorCode doJacobianFinalize();

extern PetscErrorCode spinupCalc(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscInt *neval, PetscReal *nm);
extern PetscErrorCode checkEquilibrium(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscBool doOutput, PetscBool runExternalModel, PetscReal *nm);
extern PetscErrorCode runOnExternalSignal(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state);
