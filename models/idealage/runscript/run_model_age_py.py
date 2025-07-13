import petsc4py
import sys

petsc4py.init(sys.argv)

from petsc4py import PETSc

# Always abort on petsc errors
# PETSc.Sys.pushErrorHandler("mpiabort")

import tmm4py as TMM
from time import time

from Age import Age

OptDB = PETSc.Options()
OptDB.setValue('t0',0.0)
OptDB.setValue('iter0',0)
OptDB.setValue('me','Ae1')
OptDB.setValue('mi','Ai1')
OptDB.setValue('deltat_clock',0.0009132420091324)
OptDB.setValue('max_steps',10950)
OptDB.setValue('periodic_matrix',None)
OptDB.setValue('matrix_cycle_period',1.0)
OptDB.setValue('matrix_num_per_period',12)
OptDB.setValue('biogeochem_deltat',28800.0)
OptDB.setValue('matrix_periodic_times_file','periodic_times_365d.bin')
OptDB.setValue('prefixes','state1_,state2_')
OptDB.setValue('state1_numtracers',1)
OptDB.setValue('state1_external_forcing',None)
OptDB.setValue('state1_i','ageini.petsc')
OptDB.setValue('state1_write_time_steps',1095)
OptDB.setValue('state1_o','age_1.petsc')
OptDB.setValue('state1_pickup_out','pickup_1.petsc')
OptDB.setValue('state1_time_file','output_time_1.txt')
OptDB.setValue('state2_numtracers',1)
OptDB.setValue('state2_external_forcing',None)
OptDB.setValue('state2_i','ageini.petsc')
OptDB.setValue('state2_write_time_steps',1095)
OptDB.setValue('state2_o','age_2.petsc')
OptDB.setValue('state2_pickup_out','pickup_2.petsc')
OptDB.setValue('state2_time_file','output_time_2.txt')

startTime = time()

Iter0, maxSteps, time0, deltaTClock = TMM.initialize()

if OptDB.hasName("prefixes"):
    prefixes = OptDB.getString("-prefixes").split(",")
else:
    prefixes = [""]

numStates = len(prefixes)

states = [None] * numStates

models = [Age(idx+1, prefixes[idx]) for idx in range(numStates)]

for i in range(numStates):
    states[i] = TMM.TMMState()
    states[i].create()
    states[i].setFromOptions(prefix=prefixes[i], doOutput=True)
    if states[i].config['useExternalForcing']:
      states[i].setIniExternalForcingFunction(models[i].iniAgeTracer)
      states[i].setCalcExternalForcingFunction(models[i].calcAgeTracer)
      states[i].setFinExternalForcingFunction(models[i].finalizeAgeTracer)
    
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
