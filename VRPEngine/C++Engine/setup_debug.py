from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy


extensions = [
    Extension('cpp_lower_bounds', ['cpp_lower_bounds.pyx', 'lower_bounds.cpp', 'baldacci.cpp'],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++14','-mmacosx-version-min=10.9','-g'],
              language='c++'
              ),
]

setup(
    ext_modules=cythonize(extensions),
    gdb_debug=True,
    extra_compile_args=[ '-mmacosx-version-min=10.9','-g'],
    extra_link_args=['-g'],
)
