from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy


extensions = [
    Extension('cpp_lower_bounds', ['cpp_lower_bounds.pyx', 'lower_bounds.cpp', 'baldacci.cpp'],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++1y','-g0','-O3'],
              language='c++'
              ),
]

setup(
    ext_modules=cythonize(extensions),
    extra_compile_args=['-O3','-g0'],
    extra_link_args=["-O3",'-g0'],
)
