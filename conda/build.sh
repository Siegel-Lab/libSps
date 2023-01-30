#!/bin/bash

# Compile C++ code

# create build directory
mkdir -p build
cd build

# compile code
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make VERBOSE=1 -j 8 sps stxxl stxxl_tool

# move compiled python library to output dir
mkdir -p $PREFIX/lib
cp sps*.so $PREFIX/lib

# move header library into correct directory
cp generated/inc/sps/* ../inc/sps
cmake --install . --prefix $PREFIX

# install the python package with pip @todo
#$PYTHON ../conda/setup.py install