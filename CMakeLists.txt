cmake_minimum_required (VERSION 3.9)

# only adjust major and minor version - patch is taken from git commit number
project("libSps" VERSION 0.2.0)

# set building python library as default
option(WITH_PYTHON "Build python library." ON)

set(DIMENSIONS_A "2" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_B "3" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_C "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_D "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_E "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")

set(ORTHOTOPE_A "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_B "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_C "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_D "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_E "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")

set(STORAGE_A "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_B "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_C "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_D "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_E "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")

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
include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)
if( supported )
    message(STATUS "IPO / LTO enabled")
    set(INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()

# include subprojects
if(WITH_PYTHON)
    set(requested_python_version 3.5)
    find_package(pybind11 REQUIRED)
endif()
add_subdirectory (contrib/stxxl EXCLUDE_FROM_ALL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# configure the version.h.in file
add_custom_target(
    sps_config_version ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_version.h.in.sh
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/version.h.in"
        "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        "${CMAKE_BINARY_DIR}/generated/inc/sps/version.h"
    VERBATIM)

add_custom_target(
    sps_config_build_time ALL
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/src/conf_build_time.h.in.sh
        "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/build_time.h.in"
        "${CMAKE_BINARY_DIR}/generated/inc/sps/build_time.h"
    VERBATIM)

# create python module if python has been found
if(WITH_PYTHON)
    pybind11_add_module( libSps src/tree.cpp )
    target_include_directories( libSps PUBLIC inc )
    target_include_directories( libSps SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )

    target_compile_definitions( libSps PRIVATE WITH_PYTHON )

    target_compile_definitions( libSps PRIVATE DIMENSIONS_A=${DIMENSIONS_A} )
    target_compile_definitions( libSps PRIVATE DIMENSIONS_B=${DIMENSIONS_B} )
    target_compile_definitions( libSps PRIVATE DIMENSIONS_C=${DIMENSIONS_C} )
    target_compile_definitions( libSps PRIVATE DIMENSIONS_D=${DIMENSIONS_D} )
    target_compile_definitions( libSps PRIVATE DIMENSIONS_E=${DIMENSIONS_E} )

    target_compile_definitions( libSps PRIVATE ORTHOTOPE_A=${ORTHOTOPE_A} )
    target_compile_definitions( libSps PRIVATE ORTHOTOPE_B=${ORTHOTOPE_B} )
    target_compile_definitions( libSps PRIVATE ORTHOTOPE_C=${ORTHOTOPE_C} )
    target_compile_definitions( libSps PRIVATE ORTHOTOPE_D=${ORTHOTOPE_D} )
    target_compile_definitions( libSps PRIVATE ORTHOTOPE_E=${ORTHOTOPE_E} )

    target_compile_definitions( libSps PRIVATE STORAGE_A=${STORAGE_A} )
    target_compile_definitions( libSps PRIVATE STORAGE_B=${STORAGE_B} )
    target_compile_definitions( libSps PRIVATE STORAGE_C=${STORAGE_C} )
    target_compile_definitions( libSps PRIVATE STORAGE_D=${STORAGE_D} )
    target_compile_definitions( libSps PRIVATE STORAGE_E=${STORAGE_E} )

    target_link_libraries( libSps PUBLIC stxxl )
    target_include_directories( libSps PUBLIC ${CMAKE_BINARY_DIR}/generated/inc)
    add_dependencies(libSps sps_config_version sps_config_build_time)
endif()


add_library( Sps INTERFACE )
target_include_directories( Sps SYSTEM INTERFACE inc )
target_include_directories( Sps SYSTEM INTERFACE ${STXXL_INCLUDE_DIRS} )
target_link_libraries( Sps INTERFACE stxxl )
target_include_directories( Sps SYSTEM INTERFACE ${CMAKE_BINARY_DIR}/generated/inc)
add_dependencies(Sps sps_config_version sps_config_build_time)

# Documentation setup
find_package(Doxygen)
find_package(Sphinx)
if (DOXYGEN_FOUND)

if(Sphinx_FOUND)

    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs/")

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs_conf/doxygen.config)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config.out)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( docs-doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )

    # configured documentation tools and intermediate build results
    set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_docs_build")
    
    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
    
    # HTML output directory
    set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs")

    # pdftex output
    set(SPHINX_TEX_DIR "${CMAKE_CURRENT_BINARY_DIR}/tex")

    add_custom_target(docs
        rm -rf ${SPHINX_CACHE_DIR} &&
        rm -rf ${BINARY_BUILD_DIR} &&
        cp -r ${CMAKE_CURRENT_SOURCE_DIR}/docs_conf ${BINARY_BUILD_DIR} &&
        ${SPHINX_EXECUTABLE}
            -E -a -q -b html
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.libSps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_HTML_DIR}" &&
        ${SPHINX_EXECUTABLE}
            -E -a -q -b latex
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.libSps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_TEX_DIR}" &&
        tectonic "${SPHINX_TEX_DIR}/*.tex"
        DEPENDS libSps docs-doxygen
        COMMENT "Building API documentation with Sphinx")

else(Sphinx_FOUND)
    message(WARNING "No Sphinx found. Documentation target not available.")
endif(Sphinx_FOUND)

else (DOXYGEN_FOUND)
    message("No Doxygen foud. Documentation target not available.")
endif (DOXYGEN_FOUND)


