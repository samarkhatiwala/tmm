import petsc4py
import sys
import numpy as np

petsc4py.init(sys.argv)

from petsc4py import PETSc

# Always abort on petsc errors
PETSc.Sys.pushErrorHandler("mpiabort")

import tmm4py as TMM
from time import time

OptDB = PETSc.Options()

startTime = time()

Iter0, maxSteps, time0, deltaTClock = TMM.initialize()

if OptDB.hasName("prefixes"):
    prefixes = OptDB.getString("-prefixes").split(",")
else:
    prefixes = [""]

numStates = len(prefixes)

states = [None] * numStates

if OptDB.hasName("module"):
  import importlib
  mod = OptDB.getString("module")
  model = importlib.import_module(mod)
  if OptDB.hasName("class"):
    cls = OptDB.getString("class")
    model = getattr(model, cls)
  if callable(model):
    models = [model(idx+1,prefixes[idx]) for idx in range(numStates)]
  else:
    models = [model for idx in range(numStates)]
else:
  models = [TMM.Stub(idx+1,prefixes[idx]) for idx in range(numStates)]

for i in range(numStates):
    states[i] = TMM.TMMState()
    states[i].create()
    states[i].setFromOptions(prefix=prefixes[i], doOutput=True)

    if states[i].config['useExternalForcing']:
      states[i].setIniExternalForcingFunction(models[i].iniExternalForcingFn)
      states[i].setCalcExternalForcingFunction(models[i].calcExternalForcingFn)
      states[i].setWriExternalForcingFunction(models[i].writeExternalForcingFn)
      states[i].setFinExternalForcingFunction(models[i].finalizeExternalForcingFn)
      states[i].setReiExternalForcingFunction(models[i].reInitializeExternalForcingFn)

    if states[i].config['doCalcBC']:
      states[i].setIniCalcBCFunction(models[i].iniCalcBCFn)
      states[i].setCalcCalcBCFunction(models[i].calcCalcBCFn)
      states[i].setWriCalcBCFunction(models[i].writeCalcBCFn)
      states[i].setFinCalcBCFunction(models[i].finalizeCalcBCFn)
      states[i].setReiCalcBCFunction(models[i].reInitializeCalcBCFn)

    if states[i].config['useMonitor']:
      states[i].setIniMonitorFunction(models[i].iniMonitorFn)
      states[i].setCalcMonitorFunction(models[i].calcMonitorFn)
      states[i].setWriMonitorFunction(models[i].writeMonitorFn)
      states[i].setFinMonitorFunction(models[i].finalizeMonitorFn)

    if states[i].config['doMisfit']:
      states[i].setIniMisfitFunction(models[i].iniMisfitFn)
      states[i].setCalcMisfitFunction(models[i].calcMisfitFn)
      states[i].setWriMisfitFunction(models[i].writeMisfitFn)
      states[i].setFinMisfitFunction(models[i].finalizeMisfitFn)

for iLoop in range(1, maxSteps + 1):
    tc = time0 + deltaTClock * (iLoop - 1)
    tf = time0 + deltaTClock * iLoop
    Iterc = Iter0 + iLoop - 1
    TMM.updateTMs(tc)
    for i in range(numStates):
        states[i].forcingUpdate(tc, Iterc, iLoop)
        states[i].timeStep(tc, Iterc, iLoop)
        states[i].timeStepPost(tf, Iterc, iLoop)
        states[i].output(tf, Iterc, iLoop)

endTime = time()

PETSc.Sys.Print("Wall time %10.5f" % (endTime-startTime))
