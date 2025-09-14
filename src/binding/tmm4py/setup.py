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


# we should fail to build if the build info file is not found
build_info_file = Path(__file__).parent / "_build_info.txt"
build_info = build_info_file.read_text().strip()
print("Reading build info from file: ", build_info)

petsc_setup, petsc_ver, numpy_ver = [x.split(":")[1] for x in build_info.split("|")]
if petsc_ver == "LATEST":
    required_petsc = "petsc4py"
else:
    required_petsc = f"petsc4py=={petsc_ver}"
if numpy_ver == "LATEST":
    required_numpy = "numpy"
else:
    required_numpy = f"numpy=={numpy_ver}"


log.info(f"PETSc setup type: {petsc_setup}")

if petsc_setup == "PIP_PETSC_NO_PETSC4PY":
    log.info("Using pip PETSc installation")
    import petsc

    petsc_dir = petsc.get_petsc_dir()
    os.environ["PETSC_DIR"] = petsc_dir
    petsc_arch = ""
    extra_install_requires = [required_petsc]
elif petsc_setup == "USER_PETSC_WITH_PETSC4PY":
    log.info("Using user PETSc installation with petsc4py")
    petsc_dir = os.environ["PETSC_DIR"]
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    # must be appended to sys.path so that the petsc4py library can be found by Cython during the build
    sys.path.append(os.path.join(petsc_dir, petsc_arch, "lib"))
    # Note: tmm4py.pyx looks for petsc4py/PETSc.pxd. For USER_PETSC_WITH_PETSC4PY and
    # petsc4py has been pip installed, this is in site-packages/petsc4py/ and NOT
    # PETSC_DIR/PETSC_ARCH/lib/petsc4py/. We append the path to site-packages so that
    # the file can be found by Cython
    sys.path.append(get_site_package_dir(""))
    # The user is respoinsible for ensuring the correct numpy version is installed and adding the petc4py to PYTHONPATH.
    # We will try and install latest numpy version but they may have to specify their own and use --no-build-isolation
    extra_install_requires = [required_numpy]
elif petsc_setup == "USER_PETSC_NO_PETSC4PY":
    log.info("Using user PETSc installation without petsc4py")
    petsc_dir = os.environ["PETSC_DIR"]
    petsc_arch = os.environ.get("PETSC_ARCH", "")
    extra_install_requires = [required_petsc]
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

log.info("Using TMM_DIR: " + tmm_root)

# We add Cython so that users can then build Python wrappers for their models 
# without having to install it themselves
extra_install_requires = ["Cython"] + extra_install_requires

log.info(f"Extra install requires: {extra_install_requires}")

tmm_lib_path = os.path.join(tmm_root, "lib")
tmm_include_path = os.path.join(tmm_root, "include")
# If we are installing from pip we have to link the tmm library to the site packages directory
# because tmm_lib_path is pointed towards the tmm version built in the temporary build enviroment.

# log the files in the lib directory of tmm_root
log.info("Files in TMM lib directory: " + str(os.listdir(tmm_lib_path)))


tmm_extra_link_path = (
    os.path.join(get_site_package_dir("tmmlib"), "lib") if pip_tmm else tmm_lib_path
)

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
    author="TMM Team",
    author_email="tmm-maint@example.com",
    maintainer="TMM Team",
    maintainer_email="tmm-maint@example.com",
    ext_modules=cythonize(extensions, language_level="3", include_path=include_paths),
    packages=["tmm4py", "tmm"],
    include_package_data=True,
    package_data={"tmm4py": ["*.pxd", "*.pyx", "*.so"]},
    **metadata,
    install_requires=extra_install_requires,
    setup_requires=[],
)
