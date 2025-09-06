from setuptools.build_meta import *
import os
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def _build_requires():

    log.info("Getting build requires for tmmlib")

    if "PETSC_DIR" in os.environ and os.environ["PETSC_DIR"] != "":
        log.info("Using existing PETSc installation to build tmmlib.")
        return []
    else:
        log.info("No PETSc in environment, installing from pip.")
        return ["petsc"]


def get_requires_for_build_wheel(config_settings=None):
    return _build_requires()


def get_requires_for_build_sdist(config_settings=None):
    return _build_requires()
