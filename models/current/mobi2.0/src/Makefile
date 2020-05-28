# (1) Set the environment variables TMMROOT and UVICESCMROOT to point to the top level 
# of the TMM and UVic ECSM 2.9 code, respectively.
# (2) Set UVICUPDATESLEVEL below to the updates 'level' of the UVic ESCM compatible with MOBI.
#
UVICUPDATESLEVEL = 07
#
TMMBASE = ${TMMROOT}/driver/current
MODELBASE = ${TMMROOT}/models/current/mobi2.0/src
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
	FC_FLAGS+=-convert big_endian -assume byterecl -r8
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=-fconvert=big-endian -fdefault-real-8
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qrealsize=8 -qzerosize
endif
ifneq (,$(findstring ftn,$(COMPILER)))
	FC_FLAGS+=-h byteswapio -s real64
endif

MYPATHS = $(TMMBASE):$(MODELBASE):$(UVICBASE1)/common:$(UVICBASE1)/mom:$(UVICBASE1)/embm:$(UVICBASE1)/sed:$(UVICBASE2)/common:$(UVICBASE2)/mom:$(UVICBASE2)/embm:$(UVICBASE2)/sed

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

EXE = tmmmobi
EXEWRITE = tmmmobiwrite

SRCTMM = tmm_main.c tmm_forward_step.c tmm_write.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCMODEL = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_mobi.c 
SRFMODEL = mobi_copy_data.F mobi_sed_copy_data.F mobi_ini.F mobi_calc.F mobi_diags.F mobi_stubs.F \
        co2calc.F file_names.F gasbc.F gosbc.F iomngr.F npzd_src.F setvbc.F \
        tracer.F set_sbc.F UVic_ESCM.F setmom.F sed.F setsed.F sediment.F
SRF90MODEL = mobi_diags_mod.F90

SRFMODELWRITE = mobi_write_ic.F mobi_ini.F iomngr.F npzd_src.F set_sbc.F \
               co2calc.F file_names.F mobi_stubs.F UVic_ESCM.F setmom.F \
               sed.F setsed.F sediment.F io.F

OPTFILE = MOBI_TMM_OPTIONS.h

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
SRFPPMODELWRITE = $(SRFMODELWRITE:.F=.f)
OBJTMM = $(SRCTMM:.c=.o)
OBJCMODEL = $(SRCMODEL:.c=.o)
OBJFMODEL = $(SRF90PPMODEL:.f90=.o) $(SRFPPMODEL:.f=.o)
OBJFMODELWRITE = $(SRF90PPMODEL:.f90=.o) $(SRFPPMODELWRITE:.f=.o)

$(OBJCMODEL) $(OBJFMODEL) $(SRF90PPMODEL) $(SRFPPMODEL) $(SRFPPMODELWRITE): $(OPTFILE)

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

$(EXE): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

$(EXEWRITE): $(OBJFMODELWRITE) chkopts
	-$(FLINKER) -o $@ $(OBJFMODELWRITE)

cleanall:
	make clean
	rm -f *.f *.f90 *.mod *.i *.i90
