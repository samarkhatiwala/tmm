include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# Set compiler-dependent flags
FPPFLAG=
COMPILER=$(word 1,$(MPIFC_SHOW))
ifneq (,$(findstring ifort,$(COMPILER)))
	FFLAGS+=-r8
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FFLAGS+=-fdefault-real-8
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FFLAGS+=-qrealsize=8 -qzerosize
endif

SRCTMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c petsc_signal_utils.c petsc_matvec_utils.c        
SRCTRACEGASES = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_tracegases.c
SRFTRACEGASES = tracegases_model.F tracegases_chem.F

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

OBJTMM = $(SRCTMM:.c=.o)
OBJCTRACEGASES = $(SRCTRACEGASES:.c=.o)
OBJFTRACEGASES = $(SRFTRACEGASES:.F=.o)

OBJALL = $(OBJTMM) $(OBJCTRACEGASES) $(OBJFTRACEGASES)

$(OBJALL) : TRACEGASES_OPTIONS.h

tmmtracegases:
tmmtracegases: $(OBJALL) chkopts
	-${CLINKER} -o $@ $(OBJALL) ${PETSC_MAT_LIB}

tmmtracegasesspinup: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmmtracegasesspinup: CPPFLAGS+=-DFORSPINUP
tmmtracegasesspinup: $(OBJALL) chkopts
	-${CLINKER} -o $@ $(OBJALL) ${PETSC_MAT_LIB}

tmmtracegasesjacobian: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmmtracegasesjacobian: CPPFLAGS+=-DFORJACOBIAN
tmmtracegasesjacobian: $(OBJALL) chkopts
	-${CLINKER} -o $@ $(OBJALL) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
