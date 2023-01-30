# This file does merely include spscpp, which is installed in the conda/lib directory by cmake

# This is necessary since pip (to my knowledge) does not provide a way to install module.so files 
# (i.e. c/c++ code compiled as a python module) without compiling the module itself. 
import os
import sys

sys.path.append(os.environ["PREFIX"] + "/lib")
sys.path.append(os.environ["CONDA_PREFIX"] + "/lib")

from spscpp import *