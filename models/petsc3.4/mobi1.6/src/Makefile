include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

STMM = tmm_main.c tmm_forward_step.c tmm_forcing_utils.c \
        tmm_profile_utils.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCLIN = linearsolve.c

SRCLINCUS = linearsolve_custom.c

SRCLINMF = linearsolve_matrixfree.c mymatmult.c petsc_matvec_utils.c

SRCIPC = apply_direct_inverse_with_implicit_model_pc.c petsc_matvec_utils.c \
        petsc_signal_utils.c
       
SRCTMM = $(STMM) tmm_external_forcing.c tmm_external_bc.c

SRCPATH = $(STMM) tmm_external_bc.c external_forcing_Pa_Th.c
       
SRCRDECAY = $(STMM) tmm_external_bc.c external_forcing_radioactive_decay.c

SRCCFC = $(STMM) tmm_external_bc.c external_forcing_ocmip_cfcs.c 
SRFCFC = ocmip_cfcs_model.F ocmip_cfcs_diagnostics.F
       
SRCKIELBGC = $(STMM) tmm_external_bc.c external_forcing_kiel_biogeochem.c 
SRFKIELBGC = BGC_MODEL.F BGC_INI.F \
        kiel_biogeochem_ini.F kiel_biogeochem_model.F \
        kiel_biogeochem_diagnostics.F insolation.F \
        kiel_biogeochem_set_params.F \
        CAR_CHEM.F CAR_INI.F

SRCINERTGAS = $(STMM) tmm_external_bc.c external_forcing_inert_gas.c
SRFINERTGAS = inert_gas_fluxes.F inert_gas_chem.F inert_gas_diagnostics.F \
        seawater.F

SRCINERTGASBC = $(STMM) tmm_external_forcing.c external_bc_inert_gas.c
SRFINERTGASBC = inert_gas_bc.F inert_gas_chem.F inert_gas_diagnostics.F \
        seawater.F

SRCOCMIPLANDATM = $(STMM) tmm_external_bc.c external_forcing_ocmip_abiotic_landatm.c
SRFOCMIPLANDATM = ocmip_abiotic_model.F carbon_chem.F ini_ocmip_abiotic_model.F \
        ocmip_abiotic_diagnostics.F landsource.F

SRCMITBGC = $(STMM) tmm_external_bc.c external_forcing_mitgcm_biogeochem.c
SRFMITBGC = alk_surfforcing.F carbon_chem.F fe_chem.F bio_export.F \
        dic_abiotic_param.F ini_mitgcm_biogeochem.F mitgcm_biogeochem_model.F \
        dic_biotic_forcing.F calcite_saturation.F \
        dic_biotic_param.F insol.F car_flux.F dic_surfforcing.F \
        o2_surfforcing.F car_flux_omega_top.F dic_surfforcing_init.F \
        phos_flux.F landsource.F \
        mitgcm_biogeochem_diagnostics.F

SRCMOBI = $(STMM) tmm_external_bc.c external_forcing_mobi.c 
SRFMOBI = mobi_copy_data.F mobi_ini.F mobi_calc.F mobi_util.F \
        co2calc.F gosbc.F file_names.F gasbc.F iomngr.F npzd_src.F

SRCTEST = test.c
SRFTEST = forttest.F

OBJTEST = $(SRCTEST:.c=.o) $(SRFTEST:.F=.o)
OBJLIN = $(SRCLIN:.c=.o)
OBJLINCUS = $(SRCCUS:.c=.o)
OBJLINMF = $(SRCLINMF:.c=.o)
OBJIPC = $(SRCIPC:.c=.o)
OBJTMM = $(SRCTMM:.c=.o)
OBJPATH = $(SRCPATH:.c=.o)
OBJRDECAY = $(SRCRDECAY:.c=.o)
OBJCFC = $(SRCCFC:.c=.o) $(SRFCFC:.F=.o)
OBJKIELBGC = $(SRCKIELBGC:.c=.o) $(SRFKIELBGC:.F=.o)
OBJINERTGAS = $(SRCINERTGAS:.c=.o) $(SRFINERTGAS:.F=.o)
OBJINERTGASBC = $(SRCINERTGASBC:.c=.o) $(SRFINERTGASBC:.F=.o)
OBJOCMIPLANDATM = $(SRCOCMIPLANDATM:.c=.o) $(SRFOCMIPLANDATM:.F=.o)
OBJMITBGC = $(SRCMITBGC:.c=.o) $(SRFMITBGC:.F=.o)
OBJMOBI = $(SRCMOBI:.c=.o) $(SRFMOBI:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

ARCH = $(shell uname)
NODE = $(shell uname -n)

FFLAGS=-extend_source
ifeq ($(NODE),rzcluster.rz.uni-kiel.de)
FFLAGS=-Mextend -r8
endif
ifeq ($(NODE),fats)
FFLAGS=-Mextend -r8
endif
ifeq ($(findstring ranger,$(NODE)),ranger)
FFLAGS=
endif

test: $(OBJTEST) chkopts
	-${CLINKER} -o $@ $(OBJTEST) ${PETSC_KSP_LIB}

linearsolve: $(OBJLIN) chkopts
	-${CLINKER} -o $@ $(OBJLIN) ${PETSC_KSP_LIB}

linearsolvemf: $(OBJLINMF) chkopts
	-${CLINKER} -o $@ $(OBJLINMF) ${PETSC_KSP_LIB}

linearsolvecustom: $(OBJLINCUS) chkopts
	-${CLINKER} -o $@ $(OBJLINCUS) ${PETSC_KSP_LIB}

tmm: $(OBJTMM) chkopts
	-${CLINKER} -o $@ $(OBJTMM) ${PETSC_MAT_LIB}

tmmpath: $(OBJPATH) chkopts
	-${CLINKER} -o $@ $(OBJPATH) ${PETSC_MAT_LIB}

tmmradioactivedecay: $(OBJRDECAY) chkopts
	-${CLINKER} -o $@ $(OBJRDECAY) ${PETSC_MAT_LIB}

tmmocmipcfcs: $(OBJCFC) chkopts
	-${CLINKER} -o $@ $(OBJCFC) ${PETSC_MAT_LIB}

directimplicitpc: $(OBJIPC) chkopts
	-${CLINKER} -o $@ $(OBJIPC) ${PETSC_KSP_LIB}

tmmkielbgc: FPPFLAGS=-DCARBON
tmmkielbgc: CPPFLAGS=-DCARBON
tmmkielbgc: $(OBJKIELBGC) chkopts
	-${CLINKER} -o $@ $(OBJKIELBGC) ${PETSC_MAT_LIB}

tmmkielbgcspinup: FPPFLAGS=-DCARBON -DAD_SAFE
tmmkielbgcspinup: CPPFLAGS=-DCARBON -DFORSPINUP
tmmkielbgcspinup: $(OBJKIELBGC) chkopts
	-${CLINKER} -o $@ $(OBJKIELBGC) ${PETSC_MAT_LIB}

tmmkielbgcjacobian: FPPFLAGS=-DCARBON -DAD_SAFE
tmmkielbgcjacobian: CPPFLAGS=-DCARBON -DFORSPINUP
tmmkielbgcjacobian: $(OBJKIELBGC) chkopts
	-${CLINKER} -o $@ $(OBJKIELBGC) ${PETSC_MAT_LIB}

tmminertgas:
tmminertgas: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}

tmminertgasspinup: FPPFLAGS=-DAD_SAFE
tmminertgasspinup: CPPFLAGS=-DFORSPINUP
tmminertgasspinup: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}
	
tmminertgasjacobian: FPPFLAGS=-DAD_SAFE
tmminertgasjacobian: CPPFLAGS=-DFORJACOBIAN
tmminertgasjacobian: $(OBJINERTGAS) chkopts
	-${CLINKER} -o $@ $(OBJINERTGAS) ${PETSC_MAT_LIB}

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

tmmocmipabiolandatm:
tmmocmipabiolandatm: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmspinup: FPPFLAGS=-DAD_SAFE
tmmocmipabiolandatmspinup: CPPFLAGS=-DFORSPINUP
tmmocmipabiolandatmspinup: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmocmipabiolandatmjacobian: FPPFLAGS=-DAD_SAFE
tmmocmipabiolandatmjacobian: CPPFLAGS=-DFORJACOBIAN
tmmocmipabiolandatmjacobian: $(OBJOCMIPLANDATM) chkopts
	-${CLINKER} -o $@ $(OBJOCMIPLANDATM) ${PETSC_MAT_LIB}

tmmmitbgc:
tmmmitbgc: $(OBJMITBGC) chkopts
	-${CLINKER} -o $@ $(OBJMITBGC) ${PETSC_MAT_LIB}

tmmmobi: FPPFLAGS=-DO_TMM -r8
tmmmobi: $(OBJMOBI) chkopts
	-${CLINKER} -o $@ $(OBJMOBI) ${PETSC_MAT_LIB}

#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
