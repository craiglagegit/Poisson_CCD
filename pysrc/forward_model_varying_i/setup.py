# setup.py
# Run with python setup.py build_ext --inplace

from distutils.core import setup, Extension
import os

import numpy
include_dirs = [numpy.get_include(), '.'] # This works but get deprecated warning

forward_module = Extension("forward", sources=["forward.cpp", "forward_convert.cpp"], include_dirs=include_dirs)

os.environ['CC'] = 'g++'

setup(name="forward",    
      ext_modules=[forward_module])
