# Set the environment variable TMM_DIR to point to the top level of the TMM code
#

MODEL_DIR ?= $(shell pwd)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SRCTMM = tmm_external_bc.c tmm_monitor.c tmm_misfit.c

SRCMODEL = external_forcing_mops_biogeochem.c
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

BUILDDIR := build
LIB := libmops.so
LIBTMM := ${TMM_DIR}/lib/libtmm.so

#
OBJTMM = $(addprefix $(BUILDDIR)/,$(SRCTMM:.c=.o))
OBJCMODEL = $(addprefix $(BUILDDIR)/,$(SRCMODEL:.c=.o))
OBJFMODEL = $(addprefix $(BUILDDIR)/,$(SRF90MODEL:.F90=.o) $(SRFMODEL:.F=.o))

INCLUDEPATHS = $(addprefix -I,$(MODEL_DIR)/src $(TMM_DIR)/include)

VPATH = $(MODEL_DIR)/src:$(TMM_DIR)/src

PROGRAM = mops
RPROGRAM = rmops
IPROGRAM = imops
IRPROGRAM = irmops

OPTFILE =

ifneq ($(strip $(OPTFILE)),)
# Make sure OPTFILE is present in the current directory
  ifeq (,$(wildcard $(OPTFILE)))
    $(error $(OPTFILE) does not exist in the local directory!)
  endif
  CPPFLAGS+=$(shell awk '/^\#define/ {print "-D"$$2}' $(OPTFILE))
  FPPFLAGS+=$(shell awk '/^\#define/ {print "$(FPPFLAG)-D"$$2}' $(OPTFILE))
endif

# Set extra linker flag for clang (MacOS)
EXTRALINKERFLAG=
ifneq (,$(findstring clang,$(C_VERSION)))
	EXTRALINKERFLAG+=-install_name @rpath/$(LIB)
endif

CFLAGS += -fvisibility=default
CFLAGS += -Wl,-flat_namespace,-undefined,suppress

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

$(OBJCMODEL) $(OBJFMODEL): $(OPTFILE)

$(BUILDDIR)/%.o $(BUILDDIR)/%.mod: %.F90 | $(BUILDDIR)
	$(FC) $(INCLUDEPATHS) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -fPIC -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

$(BUILDDIR)/%.o: %.F | $(BUILDDIR)
	$(FC) $(INCLUDEPATHS) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -fPIC -c -o $@ $<

$(BUILDDIR)/%.o: %.c | $(BUILDDIR)
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -fPIC -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# -undefined dynamic_lookup
#https://github.com/ziglang/zig/issues/3085
# https://github.com/ziglang/zig/issues/3085
$(PROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL)
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $(LIB) $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) -L$(TMM_DIR)/lib -ltmm

$(IPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(LIBTMM)
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $(LIB) $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) -L$(TMM_DIR)/lib -ltmm

$(RPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(LIBTMM)
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $(LIB) $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) -L$(TMM_DIR)/lib -ltmm

$(IRPROGRAM): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(LIBTMM)
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $(LIB) $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) -L$(TMM_DIR)/lib -ltmm

cleanall:
	rm -f $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL)
	rm -f $(BUILDDIR)/*.mod $(BUILDDIR)/*.i $(BUILDDIR)/*.i90
