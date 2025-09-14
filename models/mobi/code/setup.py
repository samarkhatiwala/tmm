from setuptools import setup, Extension
from Cython.Build import cythonize
import os
import numpy
import sys
import warnings

model="mobi"

if 'MODEL_DIR' in os.environ:
  model_root = os.environ['MODEL_DIR']
else:
  model_root = os.getcwd()
model_lib_path = model_root
model_include_path = os.path.join(model_root,'src')

# Get PETSc and TMM paths from environment variables
if 'PETSC_DIR' in os.environ:
  petsc_dir = os.environ['PETSC_DIR']
else:
  raise RuntimeError("You must set the PETSC_DIR environment variable before building this library")

if 'PETSC_ARCH' in os.environ:
  petsc_arch = os.environ['PETSC_ARCH']
else:
  raise RuntimeError("You must set the PETSC_ARCH environment variable before building this library")

if 'TMM_DIR' in os.environ:
  tmm_root = os.environ['TMM_DIR']
  tmm_lib_path = os.path.join(tmm_root,'lib')
  tmm_include_path = os.path.join(tmm_root,'include')
else:
  raise RuntimeError("You must set the TMM_DIR environment variable before building this library")

petscvariables = os.path.join(petsc_dir,petsc_arch,'lib','petsc','conf','petscvariables')
if os.path.isfile(petscvariables):
  fid = open(petscvariables, 'r')
  for line in fid:
    if line.startswith('CC '):
      CC = line.strip().split('=')[1].strip('\n').strip()
  fid.close()
  if "${PETSC_DIR}" in CC:
    CC=CC.replace('${PETSC_DIR}',petsc_dir)
else:
  warnings.warn(f"File {petscvariables} not found. Trying tp determing C compiler from CC environment variable")
  if 'CC' in os.environ:
    CC=os.environ["CC"]
  else:
    raise RuntimeError(f"Environment variable CC not set. Cannot determine C compiler")

os.environ["CC"]=CC

# Deal with MacOS/clang
os.environ["ARCHFLAGS"] = ""

include_paths=[
  os.path.join(petsc_dir, "include"),
  os.path.join(petsc_dir, petsc_arch, "include"),
  model_include_path,
  tmm_include_path
]

# Define the extension module
extensions = [
    Extension(
        name=model,
        sources=["model_tmm_interface.pyx"],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_9_API_VERSION')],
        include_dirs=[
            numpy.get_include(),
            os.path.join(petsc_dir, "include"),
            os.path.join(petsc_dir, petsc_arch, "include"),
            tmm_include_path,
            tmm_lib_path,
            model_include_path
        ],
        library_dirs=[
            os.path.join(petsc_dir, petsc_arch, "lib"),
            tmm_lib_path,
            model_lib_path,
        ],
        libraries=["petsc", "tmm", model],
        extra_compile_args=["-fPIC"],
        extra_link_args=[
            "-Wl,-rpath," + os.path.join(petsc_dir, petsc_arch, "lib"),
            "-Wl,-rpath," + tmm_lib_path,
            "-Wl,-rpath," + model_lib_path,
        ],
    )
]

setup(name=model, ext_modules=cythonize(extensions,language_level="3",include_path=include_paths))
