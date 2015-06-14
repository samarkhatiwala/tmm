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
        
SRCOCMIPLANDATM = $(STMM) tmm_external_bc.c tmm_monitor.c external_forcing_ocmip_abiotic_carbon_landatm.c
SRFOCMIPLANDATM = ocmip_abiotic_carbon_model.F carbon_chem.F ocmip_abiotic_carbon_ini.F \
        landsource.F

OBJOCMIPLANDATM = $(SRCOCMIPLANDATM:.c=.o) $(SRFOCMIPLANDATM:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

ARCH = $(shell uname)
NODE = $(shell uname -n)

tmmocmipabiolandatm:
tmmocmipabiolandatm: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmspinup: FPPFLAGS=$(XLFCPP)-DAD_SAFE
tmmocmipabiolandatmspinup: CPPFLAGS=-DFORSPINUP
tmmocmipabiolandatmspinup: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmjacobian: FPPFLAGS=$(XLFCPP)-DAD_SAFE
tmmocmipabiolandatmjacobian: CPPFLAGS=-DFORJACOBIAN
tmmocmipabiolandatmjacobian: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
