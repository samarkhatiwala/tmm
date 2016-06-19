include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

PROGRAM = tmmbgc
SPROGRAM = stmmbgc
IPROGRAM = itmmbgc
CPROGRAM = ctmmbgc
SIPROGRAM = sitmmbgc
CIPROGRAM = citmmbgc
SCPROGRAM = sctmmbgc
SCIPROGRAM = scitmmbgc

SRC   = \
	tmm_main.c \
	tmm_forward_step.c \
	tmm_forcing_utils.c \
	petsc_matvec_utils.c \
	petsc_signal_utils.c \
	external_forcing_kiel_biogeochem.c \
	tmm_external_bc.c \
	tmm_monitor.c \
	tmm_misfit.c \
	tmm_profile_utils.c
SRF   = \
	BGC_MODEL.F \
	BGC_INI.F \
        CAR_CHEM.F \
        CAR_INI.F \
	kiel_biogeochem_ini.F \
	kiel_biogeochem_model.F \
	kiel_biogeochem_diagnostics.F \
	kiel_biogeochem_set_params.F \
	insolation.F  
OBJECTS = $(SRC:.c=.o) $(SRF:.F=.o)

# use code without any options
$(PROGRAM): FPPFLAGS=
$(PROGRAM): CPPFLAGS=

# couple carbon to the P-core
$(CPROGRAM): FPPFLAGS=-DCARBON
$(CPROGRAM): CPPFLAGS=-DCARBON

# couple sediment to the P-core
$(SPROGRAM): FPPFLAGS=-DSEDIMENT
$(SPROGRAM): CPPFLAGS=-DSEDIMENT

# couple implicit profiles by Kriest & Oschlies, 2011
$(IPROGRAM): FPPFLAGS=-DIMPRO
$(IPROGRAM): FPPFLAGS=-DIMPRO

# couple both sediment and impro
$(SIPROGRAM): FPPFLAGS=-DSEDIMENT -DIMPRO
$(SIPROGRAM): CPPFLAGS=-DSEDIMENT -DIMPRO

# couple both sediment and carbon
$(SCPROGRAM): FPPFLAGS=-DSEDIMENT -DCARBON
$(SCPROGRAM): CPPFLAGS=-DSEDIMENT -DCARBON

# couple both carbon and impro
$(CIPROGRAM): FPPFLAGS=-DCARBON -DIMPRO
$(CIPROGRAM): CPPFLAGS=-DCARBON -DIMPRO

# couple all
$(SCIPROGRAM): FPPFLAGS=-DSEDIMENT -DCARBON -DIMPRO
$(SCIPROGRAM): CPPFLAGS=-DSEDIMENT -DCARBON -DIMPRO

$(PROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(CPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(SPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(IPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(SIPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(CIPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(SCPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}

$(SCIPROGRAM): $(OBJECTS) chkopts
	-${CLINKER} -o $@ $(OBJECTS) ${PETSC_MAT_LIB}




