import numpy as np
from petsc4py import PETSc
import tmm4py as TMM

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        return self.__dict__.update(d)

# This example illustrates the use of -calc_bc and creation/use of objects such 
# as PeriodicVec. We're not actually going to calculate BCs in this example but 
# read them from file and interpolate in time. (This is purely for illustration. 
# The TMM driver can directly handle this situation via the -prescribed_bc/-bc_files 
# options.)

class BC:

    p = DotDict()

    def __init__(self, efctxId, prefix):
      self.efctxId = efctxId
      self.prefix = prefix
      self.bc = DotDict()
      
    def iniCalcBCFn(self, tc, Iterc, tf, Iterf, state, *args, **kwargs):
      bc=self.bc
#     Set some data for this model instance
      bc.config=state.getConfig()
      numTracers=bc.config['numTracers']
#     Create some objects to facilitate 'calculating' BCs.
#     Instantiate timer object to facilitate 'calculation' of BCs
      bc.t=TMM.PeriodicTimer()
#     Read timer data specified on the command line with prefix 'bc_'
      bc.t.create(prefix="bc_")
#     Instantiate PeriodicVecs to facilitate reading BC data from file and 
#     interpolating in time
      bc.cbc=[TMM.PeriodicVec() for _ in range(numTracers)]
      bc.cbf=[TMM.PeriodicVec() for _ in range(numTracers)]
#     We create PeriodicVec's of type 'bc'. This object internally holds data (typically 
#     read from file) that are periodic in time and given at discrete intervals. These 
#     data are not exposed to the user but can be interpolated in time with the 'interp' 
#     method. The results are exposed via the numpy array PeriodicVec.arr (e.g., bc.cbc[0].arr). 
#     The 'create' method will automatically allocate memory for these arrays and subsequently, 
#     after interpolation, you can copy the data over to the corresponding numpy arrays 
#     state.cbc[0] which hold the BCs that are actually propagated by the TMM. (state.cbc[0] 
#     etc in turn point to the actual PETSc Vec objects internally used by the TMM.) 
#     We avoid this copy step by not allocating new memory for bc.cbc[0].arr etc but 
#     simply pointing their memory locations to the memory locations of state.cbc[0] etc 
#     via the optional 'buf' argument to the create method. 
      for itr in range(numTracers):
        bc.cbc[itr].create(typ='bc',buf=state.cbc[itr])
        bc.cbf[itr].create(typ='bc',buf=state.cbf[itr])

    def calcCalcBCFn(self, tc, Iterc, tf, Iterf, iLoop, state, *args, **kwargs):
      bc=self.bc
      numTracers=bc.config['numTracers']
#     Interpolate periodic data (read from files Cbc_01_XX etc) with timer object bc.t 
#     to current and next time steps tc and tf, respectively.
      for itr in range(numTracers):
        bc.cbc[itr].interp(tc,bc.t,f"Cbc_{itr+1:02d}_")
        bc.cbf[itr].interp(tf,bc.t,f"Cbc_{itr+1:02d}_")

    def writeCalcBCFn(self, tc, Iter, iLoop, state, *args, **kwargs):
      pass

    def finalizeCalcBCFn(self, tc, Iter, state, *args, **kwargs):
      bc=self.bc
      numTracers=bc.config['numTracers']
#     Destroy PeriodicVec's by calling the destroy method. Note: if the create method 
#     was called with the buf argument this won't free the underlying memory being 
#     pointed to by bc.cbc[0].arr. That memory is owned by state.cbc[0] and will 
#     be freed when state.destroy is called.
      for itr in range(numTracers):
        bc.cbc[itr].destroy()
        bc.cbf[itr].destroy()

    def reInitializeCalcBCFn(self, tc, Iter, iLoop, state):
      pass
