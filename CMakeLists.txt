cmake_minimum_required (VERSION 3.8)

project("kdtree" VERSION 0.0.1)

add_subdirectory (pybind11)

pybind11_add_module( kdtree kdtree.cpp )
