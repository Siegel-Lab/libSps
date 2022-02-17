cmake_minimum_required (VERSION 3.8)

project("cstree" VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory (pybind11)

pybind11_add_module( cstree src/tree.cpp )
target_include_directories( cstree PUBLIC inc )
target_compile_definitions( cstree PRIVATE WITH_PYTHON )
target_link_libraries(cstree PUBLIC -lstdc++)

add_executable ( test src/main.cpp )
target_include_directories( test PUBLIC inc )
target_link_libraries(test PUBLIC -lstdc++)
