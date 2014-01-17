from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("dist", ["dist.pyx"]),
               Extension("dseg", ["dseg.pyx"]),
               Extension("kmer", ["kmer.pyx"]),
               Extension("wig", ["wig.pyx"],
                         libraries=["genome"])]

setup(
  name = 'genome library',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [numpy.get_include(), '.'],
  ext_modules = ext_modules
)
