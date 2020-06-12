# (1) Set the environment variables TMMROOT and UVICESCMROOT to point to the top level 
#     of the TMM and UVic ECSM 2.9 code, respectively.
# (2) Set UVICESCMKIELMODS below to the top level of the Kiel modifications
# (3) Set UVICUPDATESLEVEL below to the updates 'level' of the UVic ESCM compatible with UVOK.
#
UVICESCMKIELMODS = /data2/spk/UVic_Kiel_Matrix/Kiel_Jan_2017_updates_to_UVIC2.9/source
UVICUPDATESLEVEL = 02
#
TMMBASE = ${TMMROOT}/driver/current
MODELBASE = ${TMMROOT}/models/current/uvok1.0/src
UVICBASE1 = ${UVICESCMROOT}/updates/${UVICUPDATESLEVEL}/source
UVICBASE2 = ${UVICESCMROOT}/source

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

MYPATHS = $(TMMBASE):$(MODELBASE):$(UVICESCMKIELMODS)/common:$(UVICESCMKIELMODS)/mom:$(UVICESCMKIELMODS)/embm:$(UVICBASE1)/common:$(UVICBASE1)/mom:$(UVICBASE1)/embm:$(UVICBASE2)/common:$(UVICBASE2)/mom:$(UVICBASE2)/embm

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

EXE = tmmuvok

SRCTMM = tmm_main.c tmm_forward_step.c tmm_write.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCMODEL = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_uvok.c 
SRFMODEL = uvok_copy_data.F uvok_ini.F uvok_calc.F uvok_diags.F uvok_stubs.F \
        co2calc.F file_names.F gasbc.F gosbc.F iomngr.F npzd_src.F setvbc.F \
        tracer.F UVic_ESCM.F
SRF90MODEL = uvok_diags_mod.F90

OPTFILE = UVOK_TMM_OPTIONS.h

ifneq ($(strip $(OPTFILE)),)
# Make sure OPTFILE is present in the current directory
  ifeq (,$(wildcard $(OPTFILE)))
    $(error $(OPTFILE) does not exist in the local directory!)
  endif
  CPPKEYS=$(shell awk '/^\#define/ {print "-D"$$2}' $(OPTFILE))
  CPPFLAGS+=$(CPPKEYS)
endif

SRFPPMODEL = $(SRFMODEL:.F=.f)
SRF90PPMODEL = $(SRF90MODEL:.F90=.f90)
OBJTMM = $(SRCTMM:.c=.o)
OBJCMODEL = $(SRCMODEL:.c=.o)
OBJFMODEL = $(SRF90PPMODEL:.f90=.o) $(SRFPPMODEL:.f=.o)

$(OBJCMODEL) $(OBJFMODEL) $(SRF90PPMODEL) $(SRFPPMODEL): $(OPTFILE)

smallf: $(SRF90PPMODEL) $(SRFPPMODEL) $(SRFPPMODELWRITE)
CPPCMD = cat $< |  cpp  -traditional -P $(INCLUDEPATHS) $(CPPKEYS)

.F.f:
	$(CPPCMD)  > $@

.F90.f90:
	$(CPPCMD)  > $@

%.o %.mod: %.f90
	$(FC) $(FC_FLAGS) $(FFLAGS) -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

# %.o: %.F
# 	$(FC) $(INCLUDEPATHS) $(FC_FLAGS) $(FFLAGS) $(FCPPFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -c -o $@ $<

$(EXE): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL)
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

cleanall:
	make clean
	rm -f *.f *.f90 *.mod *.i *.i90
