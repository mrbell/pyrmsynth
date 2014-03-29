from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
include_gsl_dir = "/usr/include/"
lib_gsl_dir = "/usr/lib/"

VERSION = '1.2.1'

# python setup.py build_ext --inplace

# Note that I need to include gslcblas otherwise I get import errors!!!
ext = Extension("grid_tools", ["grid_tools.pyx"], include_dirs=
                [numpy.get_include(), include_gsl_dir],
                library_dirs=[lib_gsl_dir],
                libraries=["gsl", "gslcblas"])

setup(name="rm_tools",
      description="rm_tools package from pyrmsynth",
      author="Michael Bell",
      author_email='mrbell@mpa-garching.mpg.de',
      packages=['rm_tools'],
      package_dir={'rm_tools': ''},
      version=VERSION,
      ext_modules=[ext], cmdclass={'build_ext': build_ext})
