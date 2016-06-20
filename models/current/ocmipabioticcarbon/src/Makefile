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

STMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCOCMIPLANDATM = $(STMM) tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_ocmip_abiotic_carbon_landatm.c
SRFOCMIPLANDATM = ocmip_abiotic_carbon_model.F carbon_chem.F ocmip_abiotic_carbon_ini.F \
        landsource.F

OBJOCMIPLANDATM = $(SRCOCMIPLANDATM:.c=.o) $(SRFOCMIPLANDATM:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

tmmocmipabiolandatm:
tmmocmipabiolandatm: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmspinup: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmmocmipabiolandatmspinup: CPPFLAGS+=-DFORSPINUP
tmmocmipabiolandatmspinup: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmjacobian: FPPFLAGS+=$(FPPFLAG)-DAD_SAFE
tmmocmipabiolandatmjacobian: CPPFLAGS+=-DFORJACOBIAN
tmmocmipabiolandatmjacobian: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
