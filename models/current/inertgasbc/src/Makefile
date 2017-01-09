include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/petscvariables

# Set compiler-dependent flags
FPPFLAG=
ifeq ($(strip $(MPIFC_SHOW)),)
	COMPILER=$(FC)
else
	COMPILER=$(word 1,$(MPIFC_SHOW))
endif

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

STMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCINERTGASBC = $(STMM) tmm_external_forcing.c tmm_monitor.c tmm_misfit.c external_bc_inert_gas.c
SRFINERTGASBC = inert_gas_bc.F inert_gas_chem.F inert_gas_diagnostics.F \
        seawater.F

OBJINERTGASBC = $(SRCINERTGASBC:.c=.o) $(SRFINERTGASBC:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

tmminertgasbc:
tmminertgasbc: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}

tmminertgasbcspinup: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmminertgasbcspinup: CPPFLAGS+=-DFORSPINUP
tmminertgasbcspinup: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}
	
tmminertgasbcjacobian: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmminertgasbcjacobian: CPPFLAGS+=-DFORJACOBIAN
tmminertgasbcjacobian: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
