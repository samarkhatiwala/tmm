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
        
SRCINERTGASBC = $(STMM) tmm_external_forcing.c tmm_monitor.c external_bc_inert_gas.c
SRFINERTGASBC = inert_gas_bc.F inert_gas_chem.F inert_gas_diagnostics.F \
        seawater.F

OBJINERTGASBC = $(SRCINERTGASBC:.c=.o) $(SRFINERTGASBC:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

ARCH = $(shell uname)
NODE = $(shell uname -n)

tmminertgasbc:
tmminertgasbc: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}

tmminertgasbcspinup: FPPFLAGS=-DAD_SAFE
tmminertgasbcspinup: CPPFLAGS=-DFORSPINUP
tmminertgasbcspinup: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}
	
tmminertgasbcjacobian: FPPFLAGS=-DAD_SAFE
tmminertgasbcjacobian: CPPFLAGS=-DFORJACOBIAN
tmminertgasbcjacobian: $(OBJINERTGASBC) chkopts
	-${CLINKER} -o $@ $(OBJINERTGASBC) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
