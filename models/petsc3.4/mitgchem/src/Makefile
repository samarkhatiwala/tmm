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
        
SRCMITGCHEM = $(STMM) tmm_external_bc.c external_forcing_mitgchem_dic.c
SRFMITGCHEM = gchem_forcing_sep.F gchem_readparms.F gchem_init_fixed.F \
            gchem_init_vari.F dic_atmos.F dic_ini_forcing.F  dic_read_pickup.F \
            dic_surfforcing.F dic_biotic_diags.F  dic_diagnostics_init.F dic_init_fixed.F \
            dic_readparms.F dic_surfforcing_init.F dic_biotic_forcing.F \
            dic_init_varia.F dic_biotic_init.F dic_ini_atmos.F dic_read_co2_pickup.F \
            alk_surfforcing.F carbon_chem.F fe_chem.F bio_export.F \
            calcite_saturation.F insol.F car_flux.F o2_surfforcing.F car_flux_omega_top.F \
            phos_flux.F \
            different_multiple.F print.F utils.F write_utils.F open_copy_data_file.F \
            nml_change_syntax.F lef_zero.F mds_reclen.F mdsfindunit.F mds_flush.F \
            packages_unused_msg.F get_periodic_interval.F mitgcm_stubs.F \
            mitgchem_ini.F mitgchem_model.F mitgchem_copy_data.F mitgchem_diagnostics.F landsource.F

#OBJMITGCHEM = $(SRCMITGCHEM:.c=.o) $(SRFMITGCHEM:.F=.o)

.SUFFIXES: 
.SUFFIXES: .o .f .p .F .c

ARCH = $(shell uname)
NODE = $(shell uname -n)

SRFPPMITGCHEM = $(SRFMITGCHEM:.F=.f)
OBJMITGCHEM = $(SRCMITGCHEM:.c=.o) $(SRFPPMITGCHEM:.f=.o)

external_forcing_mitgchem_dic.o: external_forcing_mitgchem_dic.c DIC_OPTIONS.h DIC_OPTIONS_TMM.h
DIC_OPTIONS_TMM.h:
	cat DIC_OPTIONS.h | ./stripFortranComments.sh > $@

smallf: $(SRFPPMITGCHEM)
CPPCMD = cat $< |  cpp  -traditional -P -DWORDLENGTH=4 -DALLOW_TMM | ./set64bitConst.sh
.F.f:
	$(CPPCMD)  > $@

tmmmitgchemdic: 
tmmmitgchemdic: $(OBJMITGCHEM) chkopts
	-${CLINKER} -o $@ $(OBJMITGCHEM) ${PETSC_MAT_LIB}


#.F.o:
#	-${FLINKER} -c $<

#.c.o:
#	-${CLINKER} -c ${PETSC_MAT_LIB} $<  

depend:
	/usr/X11R6/bin/makedepend $(SRCF) $(SRCC)

#clean: 
#	rm -f core $(OBJ)
