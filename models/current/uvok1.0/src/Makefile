# Set this to the top level path of the Kiel and UVic code base
BASE1 = /data2/spk/UVic_Kiel_Matrix/Kiel_Jan_2017_updates_to_UVIC2.9/source
BASE2 = /data2/spk/UVic_ESCM/2.9/updates/02/source
BASE3 = /data2/spk/UVic_ESCM/2.9/source

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

MYPATHS = $(BASE1)/common:$(BASE1)/mom:$(BASE1)/embm:$(BASE2)/common:$(BASE2)/mom:$(BASE2)/embm:$(BASE3)/common:$(BASE3)/mom:$(BASE3)/embm
   
INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

SRCTMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCUVOK = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_uvok.c 
SRFUVOK = uvok_copy_data.F uvok_ini.F uvok_calc.F uvok_diags.F uvok_stubs.F \
        co2calc.F file_names.F gasbc.F gosbc.F iomngr.F npzd_src.F setvbc.F \
        tracer.F UVic_ESCM.F
SRF90UVOK = uvok_diags_mod.F90

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

SRFPPUVOK = $(SRFUVOK:.F=.f) 
OBJTMM = $(SRCTMM:.c=.o)
OBJCUVOK = $(SRCUVOK:.c=.o)
OBJFUVOK = $(SRF90UVOK:.F90=.o) $(SRFPPUVOK:.f=.o)

CPPKEYS=$(shell awk '/^\#define/ {print "-D"$$2}' UVOK_TMM_OPTIONS.h)
FPPKEYS=$(shell awk '/^\#define/ {print "$(FPPFLAG)-D"$$2}' UVOK_TMM_OPTIONS.h)

$(OBJCUVOK) $(OBJFUVOK) $(SRFPPUVOK): UVOK_TMM_OPTIONS.h

smallf: $(SRFPPUVOK)
CPPCMD = cat $< |  cpp  -traditional -P $(INCLUDEPATHS) $(CPPKEYS)
.F.f:
	$(CPPCMD)  > $@

%.o %.mod: %.F90
	$(FC) $(INCLUDEPATHS) $(FPPKEYS) $(FC_FLAGS) -c -o $$(echo $@ | sed -e 's/mod$$/o/') $(FC_MODULE_OUTPUT_FLAG) $(dir $@) $<

tmmuvok: $(OBJTMM) $(OBJCUVOK) $(OBJFUVOK) chkopts
	-${CLINKER} -o $@ $(OBJTMM) $(OBJCUVOK) $(OBJFUVOK) ${PETSC_MAT_LIB}

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
