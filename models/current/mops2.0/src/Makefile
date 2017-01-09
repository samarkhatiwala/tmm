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
	misfit_mops_biogeochem.c \
	tmm_profile_utils.c \
	tmm_timer.c
SRF   = \
	BGC_MODEL.F \
	BGC_INI.F \
        CAR_CHEM.F \
        CAR_INI.F \
	mops_biogeochem_ini.F \
	mops_biogeochem_model.F \
	mops_biogeochem_diagnostics.F \
	mops_biogeochem_set_params.F \
	mops_biogeochem_misfit.F \
	mops_biogeochem_copy_data.F \
	insolation.F 
	
OBJECTS = $(SRC:.c=.o) $(SRF:.F=.o)

# option -DCARBON to couple carbon to the P-core
# option -DASCIIPARAMS to use ASCII parameter file as default

# use surface supply of buried P
$(PROGRAM): FPPFLAGS+=$(FPPFLAG)-DASCIIPARAMS
$(PROGRAM): CPPFLAGS+=-DASCIIPARAMS

# use mops with surface supply of buried P and implicit profiles
$(IPROGRAM): FPPFLAGS+=$(FPPFLAG)-DIMPRO -DASCIIPARAMS
$(IPROGRAM): CPPFLAGS+=-DIMPRO -DASCIIPARAMS

# use the original mops with river runoff of buried P
$(RPROGRAM): FPPFLAGS+=$(FPPFLAG)-DRUNOFF -DASCIIPARAMS
$(RPROGRAM): CPPFLAGS+=-DRUNOFF -DASCIIPARAMS

# use the original mops with river runoff of buried P and implicit profiles
$(IRPROGRAM): FPPFLAGS+=$(FPPFLAG)-DRUNOFF -DIMPRO -DASCIIPARAMS
$(IRPROGRAM): CPPFLAGS+=-DRUNOFF -DIMPRO -DASCIIPARAMS

$(PROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(IPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(RPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(IRPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

