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
	FC_FLAGS+=
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qzerosize
endif

STMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCINERTGAS = $(STMM) tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_inert_gas.c
SRFINERTGAS = inert_gas_fluxes.F inert_gas_chem.F inert_gas_diagnostics.F \
        seawater.F

OBJINERTGAS = $(SRCINERTGAS:.c=.o) $(SRFINERTGAS:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

tmminertgasgasex:
tmminertgasgasex: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}

tmminertgasgasexspinup: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmminertgasgasexspinup: CPPFLAGS+=-DFORSPINUP
tmminertgasgasexspinup: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}
	
tmminertgasgasexjacobian: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmminertgasgasexjacobian: CPPFLAGS+=-DFORJACOBIAN
tmminertgasgasexjacobian: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
