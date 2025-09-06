#!/usr/bin/env python3

"""
TMM: Transport Matrix Method Library
====================================

The Transport Matrix Method (TMM) library provides efficient numerical methods
for ocean circulation modeling and biogeochemical simulations.
"""

import re
import os
import sys
import shutil
from setuptools import setup
from setuptools.command.install import install as _install

try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ImportError:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

init_py = """\
def get_tmm_dir():
    import os
    return os.path.dirname(__file__)


def get_config():
    conf = {}
    conf['TMM_DIR'] = get_tmm_dir()
    return conf
"""

main_py = """\
if __name__ == "__main__":
    import sys
    if "--prefix" in sys.argv:
        from . import get_tmm_dir
        print(get_tmm_dir())
        del get_tmm_dir
    del sys
"""

metadata = {
    "provides": ["tmmlib"],
    "zip_safe": False,
}


def bootstrap():

    TMM_DIR = os.path.abspath(os.getcwd())
    os.environ["TMM_DIR"] = TMM_DIR
    sys.path.insert(0, os.path.join(TMM_DIR, "config"))

    pkgdir = os.path.join("config", "pypi")
    os.makedirs(pkgdir, exist_ok=True)
    for pyfile, contents in (
        ("__init__.py", init_py),
        ("__main__.py", main_py),
    ):
        with open(os.path.join(pkgdir, pyfile), "w") as fh:
            fh.write(contents)


def get_petsc_dir():
    if "PETSC_DIR" in os.environ and os.environ["PETSC_DIR"] != "":
        petsc_dir = os.environ["PETSC_DIR"]
        log.info(f"DEBUG: PETSC_DIR found in environment: {petsc_dir}")
    else:
        import petsc

        petsc_dir = petsc.get_petsc_dir()
        log.info(
            f"DEBUG: Using built-in PETSC_DIR found in pip petsc package: {petsc_dir}"
        )

    return petsc_dir


def build(dry_run=False):
    print("DEBUG: Starting TMM build")
    log.info("TMM: build")

    if dry_run:
        return
    make = shutil.which("make")
    command = [make, "all"]
    status = os.system(" ".join(command))
    if status != 0:
        raise RuntimeError(status)


def install(prefix, dry_run=False):
    log.info("TMM: install")

    if dry_run:
        return
    make = shutil.which("make")
    command = [make, "install", "PREFIX=" + prefix]
    status = os.system(" ".join(command))
    if status != 0:
        raise RuntimeError(status)


class context(object):
    def __init__(self):
        self.sys_argv = sys.argv[:]
        self.wdir = os.getcwd()

    def enter(self):
        del sys.argv[1:]
        tdir = os.environ["TMM_DIR"]
        os.environ["PETSC_DIR"] = get_petsc_dir()
        os.chdir(tdir)
        return self

    def exit(self):
        # explicitly unset PETSC_DIR incase we re-use build enviroment with --no-build-isolation
        os.environ.pop("PETSC_DIR", None)
        sys.argv[:] = self.sys_argv
        os.chdir(self.wdir)


class cmd_install(_install):

    def initialize_options(self):
        _install.initialize_options(self)

    def finalize_options(self):
        _install.finalize_options(self)
        self.install_lib = self.install_platlib
        self.install_libbase = self.install_lib
        self.old_and_unmanageable = True

    def run(self):
        root_dir = os.path.abspath(self.install_lib)
        prefix = os.path.join(root_dir, "tmmlib")

        ctx = context().enter()
        try:
            build(self.dry_run)
            install(prefix, self.dry_run)
        finally:
            ctx.exit()

        self.outputs = []
        if os.path.exists(prefix):
            for dirpath, dirnames, filenames in os.walk(prefix):
                for fn in filenames:
                    self.outputs.append(os.path.join(dirpath, fn))

        _install.run(self)

    def get_outputs(self):
        outputs = getattr(self, "outputs", [])
        outputs += _install.get_outputs(self)
        return outputs


class cmd_bdist_wheel(_bdist_wheel):

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False
        self.build_number = None

    def get_tag(self):
        plat_tag = super().get_tag()[-1]
        return (self.python_tag, "none", plat_tag)


def version():
    # Read version from include/tmm/tmmversion.h
    version_re = {
        "major": re.compile(r"#define\s+TMM_VERSION_MAJOR\s+(\d+)"),
        "minor": re.compile(r"#define\s+TMM_VERSION_MINOR\s+(\d+)"),
        "subminor": re.compile(r"#define\s+TMM_VERSION_SUBMINOR\s+(\d+)"),
    }
    tmmversion_h = os.path.join("include", "tmm", "tmmversion.h")
    if os.path.exists(tmmversion_h):
        try:
            with open(tmmversion_h, "r") as f:
                data = f.read()
            major = int(version_re["major"].search(data).groups()[0])
            minor = int(version_re["minor"].search(data).groups()[0])
            subminor = int(version_re["subminor"].search(data).groups()[0])
            return f"{major}.{minor}.{subminor}"
        except:
            pass

    # default version if header file doesn't exist or can't be read
    return "1.0.0"


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

bootstrap()
setup(
    name="tmmlib",
    version=version() + ".dev3",
    description=description.pop(0),
    long_description="\n".join(description),
    long_description_content_type="text/x-rst",
    classifiers=classifiers.split("\n")[1:-1],
    keywords=["TMM", "Ocean", "Transport", "Matrix"],
    platforms=["POSIX"],
    license="MIT",
    url="https://github.com/example/tmm",
    download_url=None,
    author="TMM Team",
    author_email="tmm-maint@example.com",
    maintainer="TMM Team",
    maintainer_email="tmm-maint@example.com",
    packages=["tmmlib"],
    package_dir={"tmmlib": "config/pypi"},
    cmdclass={
        "install": cmd_install,
        "bdist_wheel": cmd_bdist_wheel,
    },
    install_requires=[],
    setup_requires=[],
    **metadata,
)
