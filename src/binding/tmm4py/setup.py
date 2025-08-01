#!/usr/bin/env python3

"""
TMM4Py: Python Interface for Transport Matrix Method Library
============================================================

TMM4Py provides Python bindings for the Transport Matrix Method (TMM) library,
enabling efficient numerical methods for ocean circulation modeling and
biogeochemical simulations through a Python interface.
"""

import re
import os
import sys
import numpy
import warnings
from setuptools import setup, Extension
from Cython.Build import cythonize
import site
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

metadata = {
    "provides": ["tmm4py"],
    "zip_safe": False,
}


description = __doc__.split("\n")[1:-1]
del description[1:3]

classifiers = """
Development Status :: 4 - Beta
Intended Audience :: Developers
Intended Audience :: Science/Research
Operating System :: POSIX
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries
"""


def version():
    # Read version from include/tmm/tmmversion.h
    version_re = {
        "major": re.compile(r"#define\s+TMM_VERSION_MAJOR\s+(\d+)"),
        "minor": re.compile(r"#define\s+TMM_VERSION_MINOR\s+(\d+)"),
        "subminor": re.compile(r"#define\s+TMM_VERSION_SUBMINOR\s+(\d+)"),
    }
    tmmversion_h = os.path.join("..", "..", "..", "include", "tmm", "tmmversion.h")
    if os.path.exists(tmmversion_h):
        try:
            with open(tmmversion_h, "r") as f:
                data = f.read()
            major = int(version_re["major"].search(data).groups()[0])
            minor = int(version_re["minor"].search(data).groups()[0])
            subminor = int(version_re["subminor"].search(data).groups()[0])
            return "%d.%d.%d" % (major, minor, subminor)
        except:
            pass

    # default version if header file doesn't exist or can't be read
    return "1.0.0"


def get_site_package_dir(package_name):
    """Find out where the final package will be installed in site-packages"""
    return os.path.join(site.getsitepackages()[0], package_name)


# Get PETSc and TMM paths from environment variables
if "PETSC_DIR" in os.environ and os.environ["PETSC_DIR"] != "":
    petsc_dir = os.environ["PETSC_DIR"]
else:
    import petsc

    petsc_dir = petsc.get_petsc_dir()
    os.environ["PETSC_DIR"] = petsc_dir

log.info("Using PETSC_DIR: " + petsc_dir)

if "PETSC_ARCH" in os.environ:
    petsc_arch = os.environ["PETSC_ARCH"]
else:
    log.error("PETSC_ARCH not set. Using default PETSC_ARCH")
    petsc_arch = ""

if "TMM_DIR" in os.environ and os.environ["TMM_DIR"] != "":
    tmm_root = os.environ["TMM_DIR"]
else:
    import tmmlib

    tmm_root = tmmlib.get_tmm_dir()
    os.environ["TMM_DIR"] = tmm_root

log.info("Using TMM_DIR: " + tmm_root)

tmm_lib_path = os.path.join(tmm_root, "lib")
tmm_include_path = os.path.join(tmm_root, "include")

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
        f"File {petscvariables} not found. Trying tp determing C compiler from CC environment variable"
    )
    if "CC" in os.environ:
        CC = os.environ["CC"]
    else:
        raise RuntimeError(
            f"Environment variable CC not set. Cannot determine C compiler"
        )

os.environ["CC"] = CC

include_paths = [
    os.path.join(petsc_dir, "include"),
    os.path.join(petsc_dir, petsc_arch, "include"),
    tmm_include_path,
]

extensions = [
    Extension(
        name="tmm4py",
        sources=["tmm4py.pyx"],
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
            "-Wl,-rpath," + os.path.join(get_site_package_dir("tmmlib"), "lib"),
        ],
    )
]


setup(
    name="tmm4py",
    version=version(),
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
    **metadata,
    install_requires=[
        "petsc4py",
        "tmmlib @ file:///home/jamie/repos/tmm_repo/tmm",
        "mpich",
    ],
    setup_requires=[],
)
