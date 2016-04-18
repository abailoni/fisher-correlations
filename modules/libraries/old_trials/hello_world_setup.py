from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import cython_gsl


setup(
      name='Hello World Trial',
      include_dirs = [cython_gsl.get_include()],
      ext_modules=[Extension('hello_world', ['hello_world.pyx'],
                             libraries=cython_gsl.get_libraries(),
                             library_dirs=[cython_gsl.get_library_dir()],
                             include_dirs=[cython_gsl.get_cython_include_dir()])],
      cmdclass={'build_ext': build_ext}
)
