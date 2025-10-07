import os
import logging
from pathlib import Path
import subprocess

# Import the original setuptools backend with an alias
from setuptools import build_meta as _build_meta

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

PY_TMMLIB_VER = "3.2.1"


def _determine_petsc_setup_type():
    """Determines if the the PETSc version is a user-defined setup or we are using the pip PETSc version"""

    setup = "UNDETERMINED"
    log.info("Determining PETSc setup type")

    # (1) If PETSc is installed by the user via pip it won't be visible to us because of
    #     pip's build isolation. If a user wants to use this pip-installed PETSc they
    #     should indicate it by setting PETSC_DIR/PETSC_ARCH. Otherwise, if PETSC_DIR is
    #     not set we assume PETSc is not installed and we need to install it by setting
    #     setup=PIP_PETSC_NO_PETSC4PY.
    # (2) If petsc4py was installed by the user during installation of PETSc from source
    #     (by specifying --with-petsc4py during configuration), that petsc4py won't be
    #     visible to us even if PYTHONPATH is set. If a user wants to use this petsc4py,
    #     they should indicate it by setting the environment variable PETSC_HAS_PETSC4PY=1

    if os.environ.get("PETSC_DIR", "") == "":
        setup = "PIP_PETSC_NO_PETSC4PY"
    else:
        if os.environ.get("PETSC_HAS_PETSC4PY", "") == "1":
            setup = "USER_PETSC_WITH_PETSC4PY"
        else:
            setup = "USER_PETSC_NO_PETSC4PY"

    return setup


def _determine_petsc_ver(petsc_setup):
    if petsc_setup == "USER_PETSC_NO_PETSC4PY":
        petsc_dir = os.environ["PETSC_DIR"]
        log.info(
            "PETSC_VER is not set. Attempting to determine version of installed PETSc ..."
        )
        try:
            verscript = os.path.join(petsc_dir, "lib", "petsc", "bin", "petscversion")
            petsc_ver = subprocess.run(
                [f"{verscript}"], capture_output=True, text=True
            ).stdout.strip()
            log.info(f"PETSc version {petsc_ver} detected")
        except Exception as e:
            log.error(f"Error determining PETSc version: {e}")
            petsc_ver = "LATEST"
            log.info(f"PETSc version could not be detected. Using latest version")
    else:  # PIP_PETSC_NO_PETSC4PY
        petsc_ver = "LATEST"
        log.info(
            f"PETSC_VER is not set and no installed PETSc found. Using latest version"
        )
    return petsc_ver


def _extra_build_requires():

    # If tmmlib is installed by the user via pip it won't be visible to us because of
    # pip's build isolation. If a user wants to use this pip-installed tmmlib they
    # should indicate it by setting TMM_DIR. Otherwise, if TMM_DIR is not set we
    # assume tmmlib is not installed and we need to install by making it a build
    # requirement.

    log.info("Getting extra build requires")

    build_requires = []
    if os.environ.get("TMM_DIR", "") == "":
        build_requires.append(f"tmmlib=={PY_TMMLIB_VER}")

    numpy_ver = os.environ.get("NUMPY_VER", "")
    if not numpy_ver:
        numpy_ver = "LATEST"
        build_requires.append("numpy")
    else:
        build_requires.append(f"numpy=={numpy_ver}")

    cython_ver = os.environ.get("Cython_VER", "")
    if not cython_ver:
        cython_ver = "LATEST"
        build_requires.append("Cython")
    else:
        build_requires.append(f"Cython=={cython_ver}")

    petsc_setup = _determine_petsc_setup_type()
    if petsc_setup == "USER_PETSC_WITH_PETSC4PY":
        log.info("petsc4py already installed")
        petsc_ver = ""
    else:
        log.info("Adding petsc4py to build requires")
        petsc_ver = os.environ.get("PETSC_VER", "")
        if not petsc_ver:
            petsc_ver = _determine_petsc_ver(petsc_setup)
        else:
            log.info(f"Using specified PETSc version {petsc_ver}")

        if petsc_ver == "LATEST":
            build_requires.append("petsc4py")
        else:
            build_requires.append(f"petsc4py=={petsc_ver}")
    # Write everything to a build_info.txt file
    # delete any previous build_info.txt file
    if (Path(__file__).parent / "_build_info.txt").exists():
        (Path(__file__).parent / "_build_info.txt").unlink()
    info_file = Path(__file__).parent / "_build_info.txt"
    info_file_str = (
        f"petsc_setup:{petsc_setup}|petsc_ver:{petsc_ver}|numpy_ver:{numpy_ver}|cython_ver:{cython_ver}"
    )
    log.info(f"Backend: Writing {info_file_str} to {info_file}")
    info_file.write_text(info_file_str)

    return build_requires


def get_requires_for_build_wheel(config_settings=None):
    return _extra_build_requires()


def get_requires_for_build_sdist(config_settings=None):
    return _extra_build_requires()


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _build_meta.build_wheel(wheel_directory, config_settings, metadata_directory)


def build_sdist(sdist_directory, config_settings=None):
    return _build_meta.build_sdist(sdist_directory, config_settings)


prepare_metadata_for_build_wheel = _build_meta.prepare_metadata_for_build_wheel
