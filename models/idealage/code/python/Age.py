import numpy as np
from petsc4py import PETSc

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        return self.__dict__.update(d)

class Age:

    p = DotDict()

    def __init__(self, efctxId, prefix):
      self.efctxId = efctxId
      self.prefix = prefix
      self.ef = DotDict()
      
    def iniAgeTracer(self, tc, Iter, state, *args, **kwargs): 
      p=Age.p
      PETSc.Sys.Print("Initializing instance =%d" % (self.efctxId))
#     First set data common to all instances of Age
      if self.efctxId==1:
        OptDB=PETSc.Options()
        try:
          p.DeltaT=OptDB.getReal('biogeochem_deltat')
        except:
          PETSc.Error("Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option")

    def calcAgeTracer(self, tc, Iter, iLoop, state, *args, **kwargs):
        state.qef[0][:]=(self.efctxId)*28800.0/(86400.*365.)

    def finalizeAgeTracer(self, tc, Iter, state, *args, **kwargs):
        pass
