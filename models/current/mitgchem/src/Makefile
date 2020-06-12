# (1) Set the environment variables TMMROOT and MITGCMROOT to point to the top level 
#     of the TMM and MITgcm code, respectively.
#
MITGCMROOT = /data2/spk/MITgcm_latest/MITgcm
#
TMMBASE = ${TMMROOT}/driver/current
MODELBASE = ${TMMROOT}/models/current/mitgchem/src

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
	FC_FLAGS+=-r8
endif
ifneq (,$(findstring gfortran,$(COMPILER)))
	FC_FLAGS+=-fdefault-real-8
endif
ifneq (,$(findstring xlf,$(COMPILER)))
	FPPFLAG+=-WF,
	FC_FLAGS+=-qrealsize=8 -qzerosize
endif
ifneq (,$(findstring ftn,$(COMPILER)))
	FC_FLAGS+=-s real64
endif

MYPATHS = $(TMMBASE):$(MODELBASE):$(MITGCMROOT)/model/src:$(MITGCMROOT)/model/inc:$(MITGCMROOT)/eesupp/src:$(MITGCMROOT)/eesupp/inc:$(MITGCMROOT)/pkg/ptracers:$(MITGCMROOT)/pkg/gchem:$(MITGCMROOT)/pkg/dic

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(MYPATHS)))

VPATH = $(MYPATHS)

EXE = tmmmitgchemdic

SRCTMM = tmm_main.c tmm_forward_step.c tmm_write.c tmm_forcing_utils.c \
        tmm_profile_utils.c tmm_timer.c petsc_signal_utils.c petsc_matvec_utils.c
        
SRCMODEL = tmm_external_bc.c tmm_monitor.c tmm_misfit.c external_forcing_mitgchem_dic.c
SRFMODEL = different_multiple.F print.F utils.F write_utils.F open_copy_data_file.F \
            nml_change_syntax.F lef_zero.F mds_reclen.F mdsfindunit.F mds_flush.F \
            packages_unused_msg.F get_periodic_interval.F \
            gchem_forcing_sep.F gchem_readparms.F gchem_init_fixed.F gchem_init_vari.F dic_atmos.F \
            dic_ini_forcing.F  dic_surfforcing.F dic_biotic_diags.F \
            dic_diagnostics_init.F dic_init_fixed.F dic_readparms.F dic_surfforcing_init.F \
            dic_biotic_forcing.F dic_init_varia.F dic_biotic_init.F dic_ini_atmos.F \
            dic_read_co2_pickup.F alk_surfforcing.F carbon_chem.F fe_chem.F bio_export.F \
            calcite_saturation.F insol.F car_flux.F o2_surfforcing.F car_flux_omega_top.F phos_flux.F \
            mitgcm_stubs.F mitgchem_ini.F mitgchem_model.F mitgchem_copy_data.F \
            mitgchem_diagnostics.F landsource.F

OPTFILE = DIC_OPTIONS.h

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
OBJTMM = $(SRCTMM:.c=.o)
OBJCMODEL = $(SRCMODEL:.c=.o)
OBJFMODEL = $(SRF90PPMODEL:.f90=.o) $(SRFPPMODEL:.f=.o)

$(OBJCMODEL) $(OBJFMODEL) $(SRF90PPMODEL) $(SRFPPMODEL): $(OPTFILE)

smallf: $(SRF90PPMODEL) $(SRFPPMODEL)
CPPCMD = cat $< |  cpp  -traditional -P $(INCLUDEPATHS) -DWORDLENGTH=4 -DALLOW_TMM | $(MITGCMROOT)/tools/set64bitConst.sh

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

$(EXE): $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL)
	-$(CLINKER) -o $@ $(OBJTMM) $(OBJCMODEL) $(OBJFMODEL) $(PETSC_MAT_LIB)

cleanall:
	make clean
	rm -f *.f *.f90 *.mod *.i *.i90
