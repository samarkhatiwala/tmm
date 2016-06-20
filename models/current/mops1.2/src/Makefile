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

PROGRAM = mops
RPROGRAM = rmops
IPROGRAM = imops
IRPROGRAM = irmops

SRC   = \
	tmm_main.c \
	tmm_forward_step.c \
	tmm_forcing_utils.c \
	petsc_matvec_utils.c \
	petsc_signal_utils.c \
	external_forcing_mops_biogeochem.c \
	tmm_external_bc.c \
	tmm_monitor.c \
	tmm_misfit.c \
	tmm_profile_utils.c
SRF   = \
	BGC_MODEL.F \
	BGC_INI.F \
        CAR_CHEM.F \
        CAR_INI.F \
	mops_biogeochem_ini.F \
	mops_biogeochem_model.F \
	mops_biogeochem_diagnostics.F \
	mops_biogeochem_set_params.F \
	mops_biogeochem_copy_data.F \
	insolation.F 
	
OBJECTS = $(SRC:.c=.o) $(SRF:.F=.o)

# option -DCARBON: always couple carbon to the P-core

# use surface supply of buried P
$(PROGRAM): FPPFLAGS+=$(FPPFLAG)-DCARBON
$(PROGRAM): CPPFLAGS+=-DCARBON

# use mops with surface supply of buried P and implicit profiles
$(IPROGRAM): FPPFLAGS+=$(FPPFLAG)-DCARBON $(FPPFLAG)-DIMPRO
$(IPROGRAM): CPPFLAGS+=-DCARBON -DIMPRO

# use the original mops with river runoff of buried P
$(RPROGRAM): FPPFLAGS+=$(FPPFLAG)-DCARBON $(FPPFLAG)-DRUNOFF
$(RPROGRAM): CPPFLAGS+=-DCARBON -DRUNOFF

# use the original mops with river runoff of buried P and implicit profiles
$(IRPROGRAM): FPPFLAGS+=$(FPPFLAG)-DCARBON $(FPPFLAG)-DRUNOFF $(FPPFLAG)-DIMPRO
$(IRPROGRAM): CPPFLAGS+=-DCARBON -DRUNOFF -DIMPRO

$(PROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(IPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(RPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(IRPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}


