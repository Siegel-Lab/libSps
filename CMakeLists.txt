cmake_minimum_required (VERSION 3.9)

project("KdpsTree" VERSION 0.0.1)

# set building pythonj library as default
option(WITH_PYTHON "Build python library." ON)

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

# Add -O0 to remove optimizations when using gcc
IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

# Enable Link time optimization
## include(CheckIPOSupported)
## check_ipo_supported(RESULT supported OUTPUT error)
## if( supported )
##     message(STATUS "IPO / LTO enabled")
##     set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
## else()
##     message(STATUS "IPO / LTO not supported: <${error}>")
## endif()

# include subprojects
if(WITH_PYTHON)
    set(requested_python_version 3.5)
    add_subdirectory (contrib/pybind11 EXCLUDE_FROM_ALL)
endif()
add_subdirectory (contrib/stxxl EXCLUDE_FROM_ALL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# create python module if python has been found
if(WITH_PYTHON)
    pybind11_add_module( Sps src/tree.cpp )
    target_include_directories( Sps PUBLIC inc )
    target_include_directories( Sps SYSTEM PUBLIC ${STXXL_INCLUDE_DIRS} )
    target_compile_definitions( Sps PRIVATE WITH_PYTHON )
    target_link_libraries( Sps PUBLIC stxxl )
    #target_link_libraries( Sps PUBLIC -lpthread )
endif()

# build a cpp only test project (this wll show any linker errors that might be hidden with the python build)
add_executable ( test src/main.cpp )
target_include_directories( test PUBLIC inc )
target_include_directories( test SYSTEM PUBLIC ${STXXL_INCLUDE_DIRS} )
target_link_libraries(test PUBLIC stxxl )
#target_link_libraries(test PUBLIC -lpthread )

