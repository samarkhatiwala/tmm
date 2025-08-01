# (1) Set the environment variables TMM_DIR and UVICESCMROOT to point to the top level 
# of the TMM and UVic ECSM 2.9 code, respectively.
# (2) Set UVICUPDATESLEVEL below to the updates 'level' of the UVic ESCM compatible with MOBI.
#

default: all

UVICUPDATESLEVEL = 10
#
MODEL_DIR ?= $(shell pwd)
UVICBASE1 = ${UVICESCMROOT}/updates/${UVICUPDATESLEVEL}/source
UVICBASE2 = ${UVICESCMROOT}/source

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SRCTMM = tmm_external_bc.c tmm_monitor.c tmm_misfit.c
        
SRCMODEL = external_forcing_mobi.c 
SRFMODEL = mobi_copy_data.F mobi_sed_copy_data.F mobi_ini.F mobi_calc.F mobi_diags.F mobi_stubs.F \
        co2calc.F DIC_sat_calc.F file_names.F gasbc.F gosbc.F iomngr.F mobi.F setvbc.F \
        tracer.F set_sbc.F UVic_ESCM.F setmom.F sed.F setsed.F sediment.F \
        protac_thor.F
SRF90MODEL = mobi_diags_mod.F90

BUILDDIR := build
LIB := libmobi.so
LIBTMM := ${TMM_DIR}/lib/libtmm.so

#
SRFPPMODEL = $(addprefix $(BUILDDIR)/,$(SRFMODEL:.F=.f))
SRF90PPMODEL = $(addprefix $(BUILDDIR)/,$(SRF90MODEL:.F90=.f90))
OBJTMM = $(addprefix $(BUILDDIR)/,$(SRCTMM:.c=.o))
OBJCMODEL = $(addprefix $(BUILDDIR)/,$(SRCMODEL:.c=.o))
OBJFMODEL = $(SRF90PPMODEL:.f90=.o) $(SRFPPMODEL:.f=.o)

INCLUDEPATHS = $(addprefix -I,$(MODEL_DIR)/src $(TMM_DIR)/include $(UVICBASE1)/common $(UVICBASE1)/mom $(UVICBASE1)/embm $(UVICBASE1)/sed $(UVICBASE2)/common $(UVICBASE2)/mom $(UVICBASE2)/embm $(UVICBASE2)/sed)

VPATH = $(MODEL_DIR)/src:$(TMM_DIR)/src:$(UVICBASE1)/common:$(UVICBASE1)/mom:$(UVICBASE1)/embm:$(UVICBASE1)/sed:$(UVICBASE2)/common:$(UVICBASE2)/mom:$(UVICBASE2)/embm:$(UVICBASE2)/sed

# Set compiler-dependent flags
findany = $(if $1,$(or $(findstring $(firstword $1),$2),$(call findany,$(wordlist 2,$(words $1),$1),$2)))
FPPFLAG=
ifneq (,$(call findany,ifort intel Intel,$(FC_VERSION)))
	FC_FLAGS+=-convert big_endian -assume byterecl -r8
endif
ifneq (,$(call findany,gfortran gcc GCC,$(FC_VERSION)))
	FC_FLAGS+=-fconvert=big-endian -fdefault-real-8
endif
ifneq (,$(call findany,xlf IBM,$(FC_VERSION)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qrealsize=8 -qzerosize
endif
ifneq (,$(call findany,ftn cray Cray CRAY,$(FC_VERSION)))
	FC_FLAGS+=-h byteswapio -s real64
endif

OPTFILE = MOBI_TMM_OPTIONS.h

ifneq ($(strip $(OPTFILE)),)
# Make sure OPTFILE is present in the current directory
  ifeq (,$(wildcard $(OPTFILE)))
    $(error $(OPTFILE) does not exist in the local directory!)
  endif
  CPPKEYS=$(shell awk '/^\#define/ {print "-D"$$2}' $(OPTFILE))
  CPPFLAGS+=$(CPPKEYS)
endif

# Set extra linker flag for clang (MacOS)
EXTRALINKERFLAG=
ifneq (,$(findstring clang,$(C_VERSION)))
	EXTRALINKERFLAG+=-install_name @rpath/$(LIB)
endif

CFLAGS += -fvisibility=default
CFLAGS += -Wl,-flat_namespace,-undefined,suppress

all: $(LIB)

$(OBJCMODEL) $(OBJFMODEL) $(SRF90PPMODEL) $(SRFPPMODEL): $(OPTFILE)

smallf: $(SRF90PPMODEL) $(SRFPPMODEL)
CPPCMD = cat $< |  cpp  -traditional -P $(INCLUDEPATHS) $(CPPKEYS)

$(BUILDDIR)/%.f: %.F
	$(CPPCMD)  > $@

$(BUILDDIR)/%.f90: %.F90
	$(CPPCMD)  > $@

$(BUILDDIR)/%.o $(BUILDDIR)/%.mod: $(BUILDDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(INCLUDEPATHS) -I$(BUILDDIR) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -fPIC -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

$(BUILDDIR)/%.o: $(BUILDDIR)/%.f | $(BUILDDIR)
	$(FC) $(INCLUDEPATHS) -I$(BUILDDIR) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -fPIC -c -o $@ $<

$(BUILDDIR)/%.o: %.c | $(BUILDDIR)
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -fPIC -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(LIB): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(LIBTMM)
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) -L$(TMM_DIR)/lib -ltmm

cleanall:
	rm -f $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(SRF90PPMODEL) $(SRFPPMODEL)
	rm -f $(BUILDDIR)/*.mod $(BUILDDIR)/*.i $(BUILDDIR)/*.i90

