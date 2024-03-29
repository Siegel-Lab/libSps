cmake_minimum_required (VERSION 3.17)

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")

# only adjust major and minor version - patch is taken from git commit number
project("libsps" VERSION ${LIB_SPS_VERSION})

##############
# Parameters #
##############

# set building python library as default
set(DIMENSIONS_A "2" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_B "3" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_C "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_D "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_E "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_F "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_G "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_H "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_I "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_J "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_K "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_L "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_M "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")
set(DIMENSIONS_N "0" CACHE STRING "Build a datastructure with this amount of dimensions (put zero to disable).")

set(ORTHOTOPE_A "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_B "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_C "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_D "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_E "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_F "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_G "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_H "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_I "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_J "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_K "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_L "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_M "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")
set(ORTHOTOPE_N "1" CACHE STRING "Datapoints in this index have this many orthotope dimensions.")

set(STORAGE_A "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_B "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_C "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_D "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_E "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_F "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_G "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_H "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_I "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_J "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_K "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_L "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_M "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")
set(STORAGE_N "0" CACHE STRING "The storage type of this index (one of Ram, Disk, Cached).")


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(WITH_STXXL_DEFAULT OFF)
    set(UNROLL_FOR_ALL_COMBINATIONS_DEFAULT ON)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(WITH_STXXL_DEFAULT ON)
    set(UNROLL_FOR_ALL_COMBINATIONS_DEFAULT ON)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(WITH_STXXL_DEFAULT OFF)
    set(UNROLL_FOR_ALL_COMBINATIONS_DEFAULT OFF)
else()
    set(WITH_STXXL_DEFAULT OFF)
    set(UNROLL_FOR_ALL_COMBINATIONS_DEFAULT OFF)
endif()

option(WITH_STXXL "Use stxxl to build the cached forms of the index" ${WITH_STXXL_DEFAULT})
option(UNROLL_FOR_ALL_COMBINATIONS "Use aggressive template unrolling. Will increase compiletime and -ram usage massiveley, but decrease runtime." ${UNROLL_FOR_ALL_COMBINATIONS_DEFAULT})

#################
# Configuration #
#################

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # enable warnings (always good)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")
    # throw error for undefined symbols
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-undefined")

    # Add -O0 to remove optimizations when using gcc
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # include(CheckCXXSourceCompiles)
    # check_cxx_source_compiles(
    # "#include <stxxl/bits/config.h>
    #  #include <cmath>

    #  inline double log2(double x){}"
    # HAVE_LOG2_DEFINITION)
    # if(NOT HAVE_LOG2_DEFINITION)
    #     add_compile_definitions(STXXL_MSVC_COMPATIBILITY_HEADER)
    # endif()
endif()

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
if(WITH_STXXL)
    add_subdirectory ( contrib/stxxl SYSTEM)
    # these things are not needed (we only want the header files)
    set_target_properties(stxxl test1 test2 PROPERTIES EXCLUDE_FROM_ALL 1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
endif(WITH_STXXL)

find_package (Python COMPONENTS Interpreter Development.Module)

if(Python_FOUND)
    set(conf_version_exec "${CMAKE_CURRENT_SOURCE_DIR}/src/conf_version.h.in.py")
    set(conf_build_time_exec "${CMAKE_CURRENT_SOURCE_DIR}/src/conf_build_time.h.in.py")
    # configure the version.h.in file
    add_custom_target(
        sps_config_version ALL
        COMMAND ${Python_EXECUTABLE} ${conf_version_exec}
            "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/version.h.in"
            "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
            "${CMAKE_BINARY_DIR}/generated/inc/sps/version.h"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        VERBATIM)

    add_custom_target(
        sps_config_build_time ALL
        COMMAND ${Python_EXECUTABLE} ${conf_build_time_exec}
            "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/build_time.h.in"
            "${CMAKE_BINARY_DIR}/generated/inc/sps/build_time.h"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        VERBATIM)
else(Python_FOUND)
    add_custom_target(sps_config_version)
    add_custom_command(TARGET sps_config_version PRE_BUILD
                     COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/generated/inc/sps/"
                     COMMAND ${CMAKE_COMMAND} -E copy 
                     "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/version.h.backup" 
                     "${CMAKE_BINARY_DIR}/generated/inc/sps/version.h")
    add_custom_target(sps_config_build_time)
    add_custom_command(TARGET sps_config_build_time PRE_BUILD
                     COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/generated/inc/sps/"
                     COMMAND ${CMAKE_COMMAND} -E copy 
                     "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/build_time.h.backup" 
                     "${CMAKE_BINARY_DIR}/generated/inc/sps/build_time.h")
endif(Python_FOUND)


########################
# Python Module Target #
########################

#set(requested_python_version 3.5)
find_package( pybind11 PATHS ${Python_SITEARCH} )
if (pybind11_FOUND)
    pybind11_add_module( sps src/tree.cpp )
    target_include_directories( sps PUBLIC inc )
    if(WITH_STXXL)
        target_include_directories( sps SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )
        target_link_libraries( sps PRIVATE stxxl )
        target_compile_definitions( sps PRIVATE "WITH_STXXL" )
    endif(WITH_STXXL)
    if(UNROLL_FOR_ALL_COMBINATIONS)
        target_compile_definitions( sps PRIVATE "UNROLL_FOR_ALL_COMBINATIONS" )
    endif(UNROLL_FOR_ALL_COMBINATIONS)

    target_compile_definitions( sps PRIVATE WITH_PYTHON )
    target_compile_definitions( sps PRIVATE CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}" )

    target_compile_definitions( sps PRIVATE DIMENSIONS_A=${DIMENSIONS_A} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_B=${DIMENSIONS_B} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_C=${DIMENSIONS_C} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_D=${DIMENSIONS_D} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_E=${DIMENSIONS_E} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_F=${DIMENSIONS_F} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_G=${DIMENSIONS_G} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_H=${DIMENSIONS_H} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_I=${DIMENSIONS_I} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_J=${DIMENSIONS_J} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_K=${DIMENSIONS_K} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_L=${DIMENSIONS_L} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_M=${DIMENSIONS_M} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_N=${DIMENSIONS_N} )

    target_compile_definitions( sps PRIVATE ORTHOTOPE_A=${ORTHOTOPE_A} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_B=${ORTHOTOPE_B} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_C=${ORTHOTOPE_C} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_D=${ORTHOTOPE_D} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_E=${ORTHOTOPE_E} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_F=${ORTHOTOPE_F} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_G=${ORTHOTOPE_G} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_H=${ORTHOTOPE_H} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_I=${ORTHOTOPE_I} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_J=${ORTHOTOPE_J} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_K=${ORTHOTOPE_K} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_L=${ORTHOTOPE_L} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_M=${ORTHOTOPE_M} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_N=${ORTHOTOPE_N} )

    target_compile_definitions( sps PRIVATE STORAGE_A=${STORAGE_A} )
    target_compile_definitions( sps PRIVATE STORAGE_B=${STORAGE_B} )
    target_compile_definitions( sps PRIVATE STORAGE_C=${STORAGE_C} )
    target_compile_definitions( sps PRIVATE STORAGE_D=${STORAGE_D} )
    target_compile_definitions( sps PRIVATE STORAGE_E=${STORAGE_E} )
    target_compile_definitions( sps PRIVATE STORAGE_F=${STORAGE_F} )
    target_compile_definitions( sps PRIVATE STORAGE_G=${STORAGE_G} )
    target_compile_definitions( sps PRIVATE STORAGE_H=${STORAGE_H} )
    target_compile_definitions( sps PRIVATE STORAGE_I=${STORAGE_I} )
    target_compile_definitions( sps PRIVATE STORAGE_J=${STORAGE_J} )
    target_compile_definitions( sps PRIVATE STORAGE_K=${STORAGE_K} )
    target_compile_definitions( sps PRIVATE STORAGE_L=${STORAGE_L} )
    target_compile_definitions( sps PRIVATE STORAGE_M=${STORAGE_M} )
    target_compile_definitions( sps PRIVATE STORAGE_N=${STORAGE_N} )

    target_include_directories( sps PUBLIC ${CMAKE_BINARY_DIR}/generated/inc )
    add_dependencies(sps sps_config_version sps_config_build_time)
else(pybind11_FOUND)
    message(WARNING "No pybind11 found. sps python module target will not be available.")
endif(pybind11_FOUND)

#################
# Documentation #
#################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Doxygen QUIET)
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
    set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs" CACHE STRING "The HTML output directory of the documentation.")
    message("SPHINX_HTML_DIR=${SPHINX_HTML_DIR}")

    # pdftex output
    set(SPHINX_TEX_DIR "${CMAKE_CURRENT_BINARY_DIR}/tex" CACHE STRING "The TEX output directory of the documentation.")

    add_custom_target(docs
        rm -rf ${SPHINX_CACHE_DIR} &&
        rm -rf ${BINARY_BUILD_DIR} &&
        cp -r ${CMAKE_CURRENT_SOURCE_DIR}/docs_conf ${BINARY_BUILD_DIR} &&
        ${SPHINX_EXECUTABLE}
            -E -a -q -b html
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.sps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_HTML_DIR}" &&
        ${SPHINX_EXECUTABLE}
            -E -a -q -b latex
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.sps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_TEX_DIR}" &&
        tectonic "${SPHINX_TEX_DIR}/*.tex"
        DEPENDS sps docs-doxygen
        COMMENT "Building API documentation with Sphinx")

else(Sphinx_FOUND)
    message(WARNING "No Sphinx found. Documentation target will not be available.")
endif(Sphinx_FOUND)

else (DOXYGEN_FOUND)
    message("No Doxygen foud. Documentation target will not be available.")
endif (DOXYGEN_FOUND)


################
# Installation #
################

include(CMakePackageConfigHelpers)

add_library( libsps INTERFACE )
target_include_directories( libsps INTERFACE 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:inc>
)
if(WITH_STXXL)
    target_link_libraries( libsps INTERFACE stxxl )
    target_compile_definitions( libsps INTERFACE "WITH_STXXL" )
endif(WITH_STXXL)
if(UNROLL_FOR_ALL_COMBINATIONS)
    target_compile_definitions( libsps INTERFACE "UNROLL_FOR_ALL_COMBINATIONS" )
endif(UNROLL_FOR_ALL_COMBINATIONS)
add_dependencies(libsps sps_config_version sps_config_build_time)

write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/libspsConfigVersion.cmake"
    VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
    ARCH_INDEPENDENT
)

install(TARGETS libsps
    EXPORT libspsTargets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION inc COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/libspsConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/libspsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/libsps
)

install(EXPORT libspsTargets DESTINATION lib/cmake/libsps)
install(FILES "${PROJECT_BINARY_DIR}/libspsConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/libspsConfig.cmake"
        DESTINATION lib/cmake/libsps)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/inc/ ${CMAKE_BINARY_DIR}/generated/inc/ DESTINATION inc)

file(GLOB format_files
     "inc/sps/*.h"
     "src/*.cpp"
)

##########
# Format #
##########

add_custom_target(format COMMAND "clang-format" "-i" "-style=file" ${format_files}
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM)
