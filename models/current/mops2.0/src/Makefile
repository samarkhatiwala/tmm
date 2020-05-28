# Set the environment variable TMMROOT to point to the top level of the TMM code
#
TMMBASE = ${TMMROOT}/driver/current
MODELBASE = ${TMMROOT}/models/current/mops2.0/src

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Set compiler-dependent flags
FPPFLAG=
ifeq ($(strip $(MPIFC_SHOW)),)
	COMPILER=$(FC)
else
	COMPILER=$(word 1,$(MPIFC_SHOW))
endif

ifneq (,$(findstring ifort,$(COMPILER)))
	FC_FLAGS+=-r8
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=-fdefault-real-8
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qrealsize=8 -qzerosize
endif
ifneq (,$(findstring ftn,$(COMPILER)))
	FC_FLAGS+=-s real64
endif

MYPATHS = $(TMMBASE):$(MODELBASE)

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

PROGRAM = mops
RPROGRAM = rmops
IPROGRAM = imops
IRPROGRAM = irmops

SRCTMM = tmm_main.c tmm_forward_step.c tmm_write.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c

SRCMODEL = tmm_external_bc.c tmm_monitor.c misfit_mops_biogeochem.c external_forcing_mops_biogeochem.c
SRFMODEL = 	BGC_MODEL.F \
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

OPTFILE =

ifneq ($(strip $(OPTFILE)),)
# Make sure OPTFILE is present in the current directory
  ifeq (,$(wildcard $(OPTFILE)))
    $(error $(OPTFILE) does not exist in the local directory!)
  endif
  CPPFLAGS+=$(shell awk '/^\#define/ {print "-D"$$2}' $(OPTFILE))
  FPPFLAGS+=$(shell awk '/^\#define/ {print "$(FPPFLAG)-D"$$2}' $(OPTFILE))
endif

# option -DCARBON to couple carbon to the P-core
# option -DASCIIPARAMS to use ASCII parameter file as default

# use surface supply of buried P
$(PROGRAM): FPPFLAGS+=$(FPPFLAG)-DASCIIPARAMS $(FPPFLAG)-DCARBON
$(PROGRAM): CPPFLAGS+=-DASCIIPARAMS -DCARBON

# use mops with surface supply of buried P and implicit profiles
$(IPROGRAM): FPPFLAGS+=$(FPPFLAG)-DIMPRO $(FPPFLAG)-DASCIIPARAMS $(FPPFLAG)-DCARBON
$(IPROGRAM): CPPFLAGS+=-DIMPRO -DASCIIPARAMS -DCARBON

# use the original mops with river runoff of buried P
$(RPROGRAM): FPPFLAGS+=$(FPPFLAG)-DRUNOFF $(FPPFLAG)-DASCIIPARAMS $(FPPFLAG)-DCARBON
$(RPROGRAM): CPPFLAGS+=-DRUNOFF -DASCIIPARAMS -DCARBON

# use the original mops with river runoff of buried P and implicit profiles
$(IRPROGRAM): FPPFLAGS+=$(FPPFLAG)-DRUNOFF $(FPPFLAG)-DIMPRO $(FPPFLAG)-DASCIIPARAMS $(FPPFLAG)-DCARBON
$(IRPROGRAM): CPPFLAGS+=-DRUNOFF -DIMPRO -DASCIIPARAMS -DCARBON

OBJTMM = $(SRCTMM:.c=.o)
OBJCMODEL = $(SRCMODEL:.c=.o)
OBJFMODEL = $(SRF90MODEL:.F90=.o) $(SRFMODEL:.F=.o)

$(OBJCMODEL) $(OBJFMODEL): $(OPTFILE)

%.o %.mod: %.F90
	$(FC) $(INCLUDEPATHS) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

%.o: %.F
	$(FC) $(INCLUDEPATHS) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -c -o $@ $<

$(PROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

$(IPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

$(RPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

$(IRPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

cleanall:
	make clean
	rm -f *.mod *.i *.i90
