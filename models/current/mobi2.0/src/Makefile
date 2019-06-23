# (1) Set the environment variables TMMROOT and UVICESCMROOT to point to the top level 
# of the TMM and UVic ECSM 2.9 code, respectively.
# (2) Set MOBIBASE below to the MOBI-TMM interface code directory. This can be different 
# from the directory from which make is invoked, which can contain modifications 
# to the main interface code.
# (3) Set UVICUPDATESLEVEL below to the updates 'level' of the UVic ESCM compatible 
# with MOBI.
MOBIBASE = ./
UVICUPDATESLEVEL = 06
#
TMMBASE = ${TMMROOT}/driver/current
BASE1 = ${UVICESCMROOT}/updates/${UVICUPDATESLEVEL}/source
BASE2 = ${UVICESCMROOT}/source

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
	FC_FLAGS+=-convert big_endian -assume byterecl -r8
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=-fconvert=big-endian -fdefault-real-8
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qrealsize=8 -qzerosize
endif

MYPATHS = $(TMMBASE):$(MOBIBASE):$(BASE1)/common:$(BASE1)/mom:$(BASE1)/embm:$(BASE2)/common:$(BASE2)/mom:$(BASE2)/embm

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

SRCTMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCMOBI = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_mobi.c 
SRFMOBI = mobi_copy_data.F mobi_ini.F mobi_calc.F mobi_diags.F mobi_stubs.F \
        co2calc.F file_names.F gasbc.F gosbc.F iomngr.F npzd_src.F setvbc.F \
        tracer.F set_sbc.F UVic_ESCM.F setmom.F
SRF90MOBI = mobi_diags_mod.F90

SRFMOBIWRITE = mobi_write_ic.F mobi_ini.F iomngr.F npzd_src.F set_sbc.F \
               co2calc.F file_names.F mobi_stubs.F UVic_ESCM.F setmom.F io.F

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

SRFPPMOBI = $(SRFMOBI:.F=.f) 
OBJTMM = $(SRCTMM:.c=.o)
OBJCMOBI = $(SRCMOBI:.c=.o)
OBJFMOBI = $(SRF90MOBI:.F90=.o) $(SRFPPMOBI:.f=.o)
SRFPPMOBIWRITE = $(SRFMOBIWRITE:.F=.f)
OBJFMOBIWRITE = $(SRF90MOBI:.F90=.o) $(SRFPPMOBIWRITE:.f=.o)

CPPKEYS=$(shell awk '/^\#define/ {print "-D"$$2}' MOBI_TMM_OPTIONS.h)
FPPKEYS=$(shell awk '/^\#define/ {print "$(FPPFLAG)-D"$$2}' MOBI_TMM_OPTIONS.h)

$(OBJCMOBI) $(OBJFMOBI) $(SRFPPMOBI) $(SRFPPMOBIWRITE): MOBI_TMM_OPTIONS.h

smallf: $(SRFPPMOBI) $(SRFPPMOBIWRITE)
CPPCMD = cat $< |  cpp  -traditional -P $(INCLUDEPATHS) $(CPPKEYS)
.F.f:
	$(CPPCMD)  > $@

%.o %.mod: %.F90
	$(FC) $(INCLUDEPATHS) $(FPPKEYS) $(FC_FLAGS) -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

%.o: %.c
	$(CC) $(INCLUDEPATHS) ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} -c -o $@ $<

tmmmobi: $(OBJTMM) $(OBJCMOBI) $(OBJFMOBI) chkopts
	-${CLINKER} -o $@ $(OBJTMM) $(OBJCMOBI) $(OBJFMOBI) ${PETSC_MAT_LIB}

tmmmobiwrite: $(OBJFMOBIWRITE) chkopts
	-${FLINKER} -o $@ $(OBJFMOBIWRITE)

cleanall:
	make clean
	rm -f *.f *.mod

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
