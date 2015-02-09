from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os.path

source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src')
header_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'include')

include_dirs=[r'/data/kmlee/.mc3/envs/py27/include',
#              r'/data/kmlee/.local/pkg/tbb43_20140724oss/include',
              r'/data/kmlee/.local/pkg/eigen3/include/eigen3',
              numpy.get_include(),
              source_path,
              header_path,
]

library_dirs = [r'/data/kmlee/.mc3/envs/py27/lib',
                r'/data/kmlee/.local/pkg/tbb43_20140724oss/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.13.0_release'
                ]
extra_compile_args = ['-std=c++11', '-O2', ]
define_macros = [('NDEBUG', 1),
                 ('EIGEN_NO_DEBUG', 1),
                 ('USE_TR1', 1),
             ]
libraries = ['tbb',]

extensions = [
    Extension('paralinalg',
              ['cython/paralinalg.pyx', ],
              language='c++',
              define_macros = define_macros,
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
