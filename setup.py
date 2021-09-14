from setuptools import find_packages
import pkgconfig
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
import sys

ED_libraries=pkgconfig.parse('mpi')['libraries']
ED_library_dirs=pkgconfig.parse('mpi')['library_dirs']
ED_include_dirs=pkgconfig.parse('mpi')['include_dirs']

ED_libraries.extend(pkgconfig.parse('edlat')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('edlat')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('edlat')['include_dirs'])

ED_libraries.extend(pkgconfig.parse('scifor')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('scifor')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('scifor')['include_dirs'])

ext1 = Extension(
    name='edlat2py',
    sources=['src/INPUT_VARS.f90','src/edlat2py/edlat2py.f90'],
    f2py_options=["--quiet"],
    libraries=ED_libraries,
    library_dirs=ED_library_dirs,
    include_dirs=ED_include_dirs,
    extra_f90_compile_args=["-O2", "-fPIC", "-ffree-line-length-none","-cpp","-D_MPI"])


setup(
    name = "edipy",
    version = "1.0.7",
    description = "EDlat python API",
    author = "Adriano Amaricci",
    author_email = "amaricci@sissa.it",
    url='https://github.com/QcmPlab/EDlattice',
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    ext_modules=[ext1])

