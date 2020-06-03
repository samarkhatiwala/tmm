# Set the environment variable TMMROOT to point to the top level of the TMM code
#
TMMBASE = ${TMMROOT}/driver/current

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
	FC_FLAGS+=
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qzerosize
endif
ifneq (,$(findstring ftn,$(COMPILER)))
	FC_FLAGS+=
endif

MYPATHS = $(TMMBASE)

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

EXE = tmm

STMM = tmm_main.c tmm_forward_step.c tmm_write.c tmm_forcing_utils.c \
       tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCTMM = $(STMM) tmm_external_forcing.c tmm_external_bc.c tmm_monitor.c tmm_misfit.c

OBJTMM = $(SRCTMM:.c=.o)

%.o: %.c
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -c -o $@ $<

$(EXE): $(OBJTMM) chkopts
	-$(CLINKER) -o $@ $(OBJTMM) $(PETSC_MAT_LIB)

cleanall:
	make clean
	rm -f *.mod *.i *.i90
