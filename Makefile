default: all

TMMBASE = .

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

STMM = tmm.c tmm_forward_step.c tmm_transport_matrices.c tmm_forcing.c tmm_output.c \
       tmm_write_extra.c tmm_forcing_utils.c \
       tmm_profile_utils.c tmm_timer.c tmm_petsc_matvec_utils.c \
	      tmm_initialize.c tmm_timestep.c tmm_finalize.c

# PREFIX ?= $(abspath TMM)
LIB := libtmm.so
SRCDIR := $(TMMBASE)/src/tmm
INCDIR := $(TMMBASE)/include/tmm
BUILDDIR := build
EXE := $(BUILDDIR)/$(LIB)

#
SRCTMM = $(addprefix $(SRCDIR)/,$(STMM))

OBJTMM := $(subst $(SRCDIR), $(BUILDDIR), $(SRCTMM:.c=.o))
# $(info $(OBJTMM))

INCLUDEPATHS = $(addprefix -I,$(subst :, ,$(INCDIR)))

CFLAGS += -fvisibility=default

# Set extra linker flag for clang (MacOS)
EXTRALINKERFLAG=
ifneq (,$(findstring clang,$(C_VERSION)))
	EXTRALINKERFLAG+=-install_name @rpath/$(LIB)
endif

all: $(OBJTMM) $(EXE)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(INCLUDEPATHS) $(PCC_FLAGS) $(CFLAGS) $(CCPPFLAGS) -fPIC -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)
	
$(EXE): $(OBJTMM)
 ifeq (${PETSCPYTHONSITE}, )
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $@ $(OBJTMM) $(PETSC_MAT_LIB)
 else
	-$(CLINKER) -shared $(EXTRALINKERFLAG) -o $@ $(OBJTMM) $(PETSC_MAT_LIB) -Wl,-rpath,${PETSCPYTHONSITE}/petsc/lib -L${PETSCPYTHONSITE}/petsc/lib
 endif

install: PREFIX ?= $(abspath TMM)
install: $(EXE)
	@echo "Installing TMM in ${PREFIX}"
	mkdir -p ${PREFIX}/lib
	mv $< ${PREFIX}/lib
	cp -p lib/tmmversion ${PREFIX}/lib/
	mkdir -p ${PREFIX}/include
	cp -p ${INCDIR}/* ${PREFIX}/include/
# Stubs
	mkdir -p ${PREFIX}/src
	cp -p ${SRCDIR}/tmm_main.c ${SRCDIR}/tmm_external_forcing.c ${SRCDIR}/tmm_external_bc.c ${SRCDIR}/tmm_monitor.c ${SRCDIR}/tmm_misfit.c ${PREFIX}/src/
	@echo "Now set the TMM_DIR environment variable to ${PREFIX}"
	@echo "To build and install tmm4py do: make tmm4py to install in site packages or "
	@echo "                                make tmm4py PREFIX=XXX to install in XXX"

tmm4py: PIP ?= pip3
tmm4py:
	rm -f src/binding/tmm4py/tmm4py/tmm4py_core.c
	rm -rf src/binding/tmm4py/build
	rm -f src/binding/tmm4py/_build_info.txt
 ifeq (${PREFIX}, )
	@echo "Installing tmm4py in site packages" 
	cd src/binding/tmm4py && ${PIP} install --no-build-isolation .
 else
	@echo "Installing tmm4py in ${PREFIX}"
	rm -rf ${PREFIX}/tmm4py*
	rm -rf ${PREFIX}/tmm
	cd src/binding/tmm4py && ${PIP} install --no-build-isolation --target ${PREFIX} .
	@echo "Now add ${PREFIX} to the PYTHONPATH environment variable"
 endif
	
cleanall:
	rm -f $(OBJTMM)
	rm -f $(BUILDDIR)/*.mod $(BUILDDIR)/*.i $(BUILDDIR)/*.i90
