from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os.path

source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src')

include_dirs=[r'/Users/kmlee/.anaconda/include',
              r'/Users/kmlee/code/koretools/include',
              r'/Users/kmlee/.brew/include/eigen3',
              r'/Users/kmlee/.brew/include',
              numpy.get_include(),
              source_path,
]

library_dirs = [r'/Users/kmlee/.anaconda/lib',
                r'/Users/kmlee/.brew/lib',
                ]

extra_compile_args = ['-std=c++11', '-O2', '-DNDEBUG', '-DEIGEN_NO_DEBUG', ]

libraries = ['tbb',]

extensions = [
    Extension('paralinalg',
              ['cython/paralinalg.pyx', ],
              language='c++',
              extra_compile_args = extra_compile_args,
              libraries = libraries,
              include_dirs = include_dirs,
              library_dirs = library_dirs)
    ]
    

setup(
    name = 'KoreTools',
    version = '1.0',
    description='Kore Tools Module',
    author='Kyungmin Lee',
    author_email='kyungmin.lee.42@gmail.com',
    ext_modules = cythonize(extensions))
