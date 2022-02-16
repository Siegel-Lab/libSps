cmake_minimum_required (VERSION 3.8)

project("cstree" VERSION 0.0.1)

add_subdirectory (pybind11)

pybind11_add_module( cstree src/tree.cpp )
target_include_directories( cstree PUBLIC inc )
