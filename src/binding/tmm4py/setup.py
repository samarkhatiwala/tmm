#!/usr/bin/env python3

"""
TMM4Py: Python Interface for Transport Matrix Method Library
============================================================

TMM4Py provides Python bindings for the Transport Matrix Method (TMM) library,
enabling efficient numerical methods for ocean circulation modeling and
biogeochemical simulations through a Python interface.
"""

import os
import sys
import numpy
import warnings
from pathlib import Path
from setuptools import setup, Extension
from Cython.Build import cythonize
import site
import logging
import sys
import subprocess

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

metadata = {
    "provides": ["tmm4py"],
    "zip_safe": False,
}


description = __doc__.split("\n")[1:-1]
del description[1:3]

classifiers = """
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: Science/Research
Operating System :: POSIX
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries
"""


def get_site_package_dir(package_name):
    """Find out where the final package will be installed in site-packages"""
    return os.path.join(site.getsitepackages()[0], package_name)


# I don't know what the best way is to ensure the correct numpy and Cython versions are used to build 
# tmm4py or make available when pip finishes. 
# If there is no pre-installed petsc4py, the petsc4py build will likely install the latest numpy/Cython. 
# However, these may not be the same as the ones optionally specified by the user with NUMPY_VER and CYTHON_VER. 
# resulting in different versions being used for petsc4py and tmm4py. Moreover, they may clash with the ones 
# already installed in the user's environment. Replacing them with the numpy/Cython used for building 
# petsc4py/tmm4py does not sound like a good idea. Ideally, we want the user to install petsc4py separately 
# and indicate it with PETSC_HAS_PETSC4PY=1. Alternatively, they can install the desired versions of numpy 
# and Cython first and install petsc4py/tmm4py with --no-build-isolation.
# For the time being I'm not adding either numpy or Cython to extra_install_requires.


# we should fail to build if the build info file is not found ... except if --no-build-isolation 
# is used
try:
    build_info_file = Path(__file__).parent / "_build_info.txt"
    build_info = build_info_file.read_text().strip()
    print("Reading build info from file: ", build_info)

    petsc_setup, petsc_ver, numpy_ver, cython_ver = [x.split(":")[1] for x in build_info.split("|")]
    if petsc_ver == "LATEST":
        required_petsc = "petsc"
    else:
        required_petsc = f"petsc=={petsc_ver}"
    if petsc_ver == "LATEST":
        required_petsc4py = "petsc4py"
    else:
        required_petsc4py = f"petsc4py=={petsc_ver}"
    if numpy_ver == "LATEST":
        required_numpy = "numpy"
    else:
        required_numpy = f"numpy=={numpy_ver}"
    if cython_ver == "LATEST":
        required_cython = "Cython"
    else:
        required_cython = f"Cython=={cython_ver}"
    # (We're not doing this. See above) We add Cython so that users can then build Python wrappers for their models 
    # without having to install it themselves
    extra_install_requires = []
except:
    # If we're here then (presumably) --no-build-isolation was used and we 
    # assume PETSc, petsc4py and all required packages are already installed
    petsc_setup = "USER_PETSC_WITH_PETSC4PY_NO_BUILD_ISOLATION"
    extra_install_requires = []

log.info(f"PETSc setup type: {petsc_setup}")

if petsc_setup == "PIP_PETSC_NO_PETSC4PY":
    log.info("Using pip PETSc installation")
    import petsc
    petsc_dir = petsc.get_petsc_dir()
    os.environ["PETSC_DIR"] = petsc_dir
    petsc_arch = ""
    extra_install_requires = extra_install_requires + [required_petsc] + [required_petsc4py]
elif petsc_setup == "USER_PETSC_WITH_PETSC4PY":
    log.info("Using user PETSc installation with petsc4py")
    petsc_dir = os.environ["PETSC_DIR"]
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    # Note: tmm4py.pyx looks for petsc4py/PETSc.pxd and the correct path needs to be 
    # appended to sys.path so that it is found by Cython during the build. If petsc4py 
    # was configured and installed from source along with PETSc, petsc4py/ will be in 
    # PETSC_DIR/PETSC_ARCH/lib/ and we first look for it there:
    p4pdirFound = False
    p4pdir = os.path.join(petsc_dir, petsc_arch, "lib", "petsc4py")
    if os.path.isdir(p4pdir):
        s = os.path.join(petsc_dir, petsc_arch, "lib")
        log.info(f"Existing petsc4py found at {p4pdir}; appending {s} to sys.path")
        sys.path.append(s)
        p4pdirFound = True
    else:
        # However, if petsc4py was pip installed, it will be in site-packages/ (but we need 
        # to determine where exactly)
        if sys.prefix != sys.base_prefix:
            # We're in a virtual environment
            log.info(f"DEBUG: virtual environment detected")
            s = site.getsitepackages()[0]
        else:
            log.info(f"DEBUG: assuming regular environment")
            s = site.getusersitepackages()
        p4pdir=os.path.join(s,"petsc4py")
        if os.path.isdir(p4pdir):
            log.info(f"Existing petsc4py found at {p4pdir}; appending {s} to sys.path")
            sys.path.append(s)
            p4pdirFound = True
    if not p4pdirFound: raise RuntimeError(f"Could not find installed petsc4py in expected locations!")
    # The user is responsible for ensuring the correct numpy version is installed and adding the petc4py to PYTHONPATH.
    # We will try and install latest numpy version but they may have to specify their own and use --no-build-isolation
    #extra_install_requires = extra_install_requires + [required_numpy]
elif petsc_setup == "USER_PETSC_NO_PETSC4PY":
    log.info("Using user PETSc installation without petsc4py")
    petsc_dir = os.environ["PETSC_DIR"]
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    extra_install_requires = extra_install_requires + [required_petsc4py]
elif petsc_setup == "USER_PETSC_WITH_PETSC4PY_NO_BUILD_ISOLATION":
    log.info("Using user PETSc installation with petsc4py and all other required packages")
    petsc_dir = os.environ["PETSC_DIR"]
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    # If petsc4py was configured with PETSc then we assume PYTHONPATH has been correctly set by the user
    # If petsc4py was installed with pip then nothing needs to be done
else:
    raise RuntimeError(f"Unknown PETSc setup type: {petsc_setup}")

log.info("Using PETSC_DIR: " + petsc_dir + " and PETSC_ARCH: " + petsc_arch)

if os.environ.get("TMM_DIR", "") != "":
    tmm_root = os.environ["TMM_DIR"]
    pip_tmm = False
else:
    import tmmlib
    pip_tmm = True
    tmm_root = tmmlib.get_tmm_dir()
    os.environ["TMM_DIR"] = tmm_root
    extra_install_requires = ["tmmlib"] + extra_install_requires

log.info(f"Using TMM_DIR: {tmm_root}")

log.info(f"Extra install requires: {extra_install_requires}")

tmm_lib_path = os.path.join(tmm_root, "lib")
tmm_include_path = os.path.join(tmm_root, "include")

# log the files in the lib directory of tmm_root
log.info("Files in TMM lib directory: " + str(os.listdir(tmm_lib_path)))

# If we are installing tmmlib with pip we have to link the tmm library to the site packages directory 
# where it will be ultimately installed because tmm_lib_path is currently pointed towards the tmmlib 
# version built in the temporary build enviroment. This is tricky to determine and quite fragile. 
# Ideally, we want the user to install TMM separately and set TMM_DIR. Or indicate the ultimate location 
# by setting the PETSCPYTHONSITE environment variable.
# For now we assume this will be in the user site packages directory. Unfortunately, this 
# location depends on whether pip was invoked from within a virtual environment or not. 
# If a virtual environment, then the install location will be site.getsitepackages()[0]. If 
# from a 'regular' environment then it will be site.getusersitepackages(). At least that's 
# what I've determined from experimenting.

if pip_tmm:
    if os.environ.get('PETSCPYTHONSITE',"")=="":
        if sys.prefix != sys.base_prefix:
            # We're in a virtual environment
            log.info(f"DEBUG: virtual environment detected")
            s = site.getsitepackages()[0]
        else:
            log.info(f"DEBUG: assuming regular environment")
            s = site.getusersitepackages()
        log.info(f"DEBUG: assuming tmmlib will be installed in {s}")
    else:
        s = os.environ["PETSCPYTHONSITE"]
        log.info(f"DEBUG: PETSCPYTHONSITE environment variable found; assuming tmmlib will be installed in {s}")
    tmm_extra_link_path = os.path.join(s,"tmmlib", "lib")
else:
    tmm_extra_link_path = tmm_lib_path

log.info("Getting TMM version")
verscript = os.path.join(tmm_root, "lib", "tmmversion")
tmm_ver = subprocess.run(
    [f"{verscript}"], capture_output=True, text=True
).stdout.strip()

petscvariables = os.path.join(
    petsc_dir, petsc_arch, "lib", "petsc", "conf", "petscvariables"
)
if os.path.isfile(petscvariables):
    fid = open(petscvariables, "r")
    for line in fid:
        if line.startswith("CC "):
            CC = line.strip().split("=")[1].strip("\n").strip()
    fid.close()
    if "${PETSC_DIR}" in CC:
        CC = CC.replace("${PETSC_DIR}", petsc_dir)
else:
    warnings.warn(
        f"File {petscvariables} not found. Trying to determine C compiler from CC environment variable"
    )
    if "CC" in os.environ:
        CC = os.environ["CC"]
    else:
        raise RuntimeError(
            f"Environment variable CC not set. Cannot determine C compiler"
        )

os.environ["CC"] = CC

# Deal with MacOS/clang
os.environ["ARCHFLAGS"] = ""

include_paths = [
    os.path.join(petsc_dir, "include"),
    os.path.join(petsc_dir, petsc_arch, "include"),
    tmm_include_path,
]

extensions = [
    Extension(
        name="tmm4py.tmm4py_core",
        sources=["tmm4py/tmm4py_core.pyx"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_9_API_VERSION")],
        include_dirs=[
            numpy.get_include(),
            os.path.join(petsc_dir, "include"),
            os.path.join(petsc_dir, petsc_arch, "include"),
            tmm_include_path,
            tmm_lib_path,
        ],
        library_dirs=[
            os.path.join(petsc_dir, petsc_arch, "lib"),
            tmm_lib_path,
        ],
        libraries=["petsc", "tmm"],
        extra_compile_args=["-fPIC"],
        extra_link_args=[
            "-Wl,-rpath," + os.path.join(get_site_package_dir("petsc"), "lib"),
            "-Wl,-rpath," + tmm_extra_link_path,
        ],
    )
]


setup(
    name="tmm4py",
    version=tmm_ver,
    description=description.pop(0),
    long_description="\n".join(description),
    long_description_content_type="text/x-rst",
    classifiers=classifiers.split("\n")[1:-1],
    keywords=["TMM", "Ocean", "Transport", "Matrix", "Python", "Bindings"],
    platforms=["POSIX"],
    license="",
    url="https://github.com/samarkhatiwala/tmm",
    download_url=None,
    author="Samar Khatiwala, Jamie Carr",
    author_email="samkat6@gmail.com",
    maintainer="Samar Khatiwala",
    maintainer_email="samkat6@gmail.com",
    ext_modules=cythonize(extensions, language_level="3", include_path=include_paths),
    packages=["tmm4py", "tmm"],
    include_package_data=True,
    package_data={"tmm4py": ["*.pxd", "*.pyx", "*.so"]},
    **metadata,
    install_requires=extra_install_requires,
    setup_requires=[],
)
