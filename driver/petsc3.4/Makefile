include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# Modify these as necessary
XLFCPP=
XLFFLAGS=
R8FLAG=-r8
# if using the IBM XLF compiler, uncomment the following lines
#XLFCPP=-WF,
#R8FLAG=-qrealsize=8
#XLFFLAGS=-qzerosize

STMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCTMM = $(STMM) tmm_external_forcing.c tmm_external_bc.c

OBJTMM = $(SRCTMM:.c=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

ARCH = $(shell uname)
NODE = $(shell uname -n)

tmm: $(OBJTMM) chkopts
	-${CLINKER} -o $@ $(OBJTMM) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
