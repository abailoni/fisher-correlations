from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys
import numpy

if sys.platform == "linux2":
    include_gsl_dir = "/usr/local/include/"
    lib_gsl_dir = "/usr/local/lib/"
elif sys.platform == "darwin":
    include_gsl_dir = "/opt/local/include/"
    lib_gsl_dir = "/opt/local/lib/"

ext = Extension('uFM', ['uFM.pyx', 'GSL_library.pxd'],
    include_dirs=[numpy.get_include(),
                  include_gsl_dir],
    library_dirs=[lib_gsl_dir],
    libraries=["gsl", "gslcblas"])

setup(
      name='Compute non-correlate Fisher Matrix',
      include_dirs = [numpy.get_include()],
      ext_modules=[ext],
      cmdclass={'build_ext': build_ext}
)
