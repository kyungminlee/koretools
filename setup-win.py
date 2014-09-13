from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import os.path

source_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src')

include_dirs=[r'C:\Program Files\Enthought\Canopy\App\appdata\canopy-1.4.1.1975.win-x86_64\include',
              r'C:\Users\kmlee\code\koretools\include',
              r'C:\Users\kmlee\code\dev\src\eigen',
              r'E:\dev\src\eigen',
              r'C:\Program Files (x86)\Intel\Composer XE\tbb\include',
              numpy.get_include(),
              source_path,
]

library_dirs = [r'C:\Program Files (x86)\Intel\Composer XE\tbb\lib\intel64\vc12',
                ]

extra_compile_args = ['/O2', '/DNDEBUG', '/DEIGEN_NO_DEBUG', '/EHsc',]

libraries = ['tbb',]

extensions = [
    Extension('koretools',
              ['cython/koretools.pyx', ],
              language='c++',
              extra_compile_args = extra_compile_args,
              libraries = libraries,
              include_dirs = include_dirs,
              library_dirs = library_dirs),
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
