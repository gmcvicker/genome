from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("trackreader", ["trackreader.pyx"],
                         libraries=["genome", "z"]),
               Extension("seq", ["seq.pyx"],
                         libraries=["genome", "z"])]

setup(
  name = 'track file parser',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [numpy.get_include(), '.'],
  ext_modules = ext_modules
)
