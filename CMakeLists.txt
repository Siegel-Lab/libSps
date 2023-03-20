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


if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(WITH_STXXL_DEFAULT OFF)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(WITH_STXXL_DEFAULT ON)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(WITH_STXXL_DEFAULT OFF)
else()
    set(WITH_STXXL_DEFAULT OFF)
endif()

#option(WITH_STXXL "Use stxxl to build the cached forms of the index" ${WITH_STXXL_DEFAULT})
option(WITH_STXXL "Use stxxl to build the cached forms of the index" OFF)

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
    #add_subdirectory ( contrib/stxxl EXCLUDE_FROM_ALL )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")
endif(WITH_STXXL)

find_package (Python COMPONENTS Interpreter Development)

if(Python_FOUND)
    set(conf_version_exec "${CMAKE_CURRENT_SOURCE_DIR}/src/conf_version.h.in.py")
    set(conf_build_time_exec "${CMAKE_CURRENT_SOURCE_DIR}/src/conf_build_time.h.in.py")
    # configure the version.h.in file
    add_custom_target(
        sps_config_version ALL
        COMMAND "python" ${conf_version_exec}
            "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/version.h.in"
            "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
            "${CMAKE_BINARY_DIR}/generated/inc/sps/version.h"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        VERBATIM)

    add_custom_target(
        sps_config_build_time ALL
        COMMAND "python" ${conf_build_time_exec}
            "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/build_time.h.in"
            "${CMAKE_BINARY_DIR}/generated/inc/sps/build_time.h"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        VERBATIM)
else(Python_FOUND)
    add_custom_target(sps_config_version)
    add_custom_command(TARGET sps_config_version PRE_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy 
                     "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/version.h.backup" 
                     "${CMAKE_BINARY_DIR}/generated/inc/sps/version.h")
    add_custom_target(sps_config_build_time)
    add_custom_command(TARGET sps_config_build_time PRE_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy 
                     "${CMAKE_CURRENT_SOURCE_DIR}/inc/sps/build_time.h.backup" 
                     "${CMAKE_BINARY_DIR}/generated/inc/sps/build_time.h")
endif(Python_FOUND)


########################
# Python Module Target #
########################

#set(requested_python_version 3.5)
find_package( pybind11 QUIET PATHS ${Python_SITEARCH} )
if (pybind11_FOUND)
    pybind11_add_module( sps src/tree.cpp )
    target_include_directories( sps PUBLIC inc )
    if(WITH_STXXL)
        target_include_directories( sps SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )
        target_compile_definitions( sps PRIVATE "WITH_STXXL" )
    endif(WITH_STXXL)

    target_compile_definitions( sps PRIVATE WITH_PYTHON )
    target_compile_definitions( sps PRIVATE CXX_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}" )

    target_compile_definitions( sps PRIVATE DIMENSIONS_A=${DIMENSIONS_A} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_B=${DIMENSIONS_B} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_C=${DIMENSIONS_C} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_D=${DIMENSIONS_D} )
    target_compile_definitions( sps PRIVATE DIMENSIONS_E=${DIMENSIONS_E} )

    target_compile_definitions( sps PRIVATE ORTHOTOPE_A=${ORTHOTOPE_A} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_B=${ORTHOTOPE_B} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_C=${ORTHOTOPE_C} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_D=${ORTHOTOPE_D} )
    target_compile_definitions( sps PRIVATE ORTHOTOPE_E=${ORTHOTOPE_E} )

    target_compile_definitions( sps PRIVATE STORAGE_A=${STORAGE_A} )
    target_compile_definitions( sps PRIVATE STORAGE_B=${STORAGE_B} )
    target_compile_definitions( sps PRIVATE STORAGE_C=${STORAGE_C} )
    target_compile_definitions( sps PRIVATE STORAGE_D=${STORAGE_D} )
    target_compile_definitions( sps PRIVATE STORAGE_E=${STORAGE_E} )

    if(WITH_STXXL)
        target_link_libraries( sps PUBLIC stxxl )
    endif(WITH_STXXL)
    target_include_directories( sps PUBLIC ${CMAKE_BINARY_DIR}/generated/inc )
    add_dependencies(sps sps_config_version sps_config_build_time)
else(pybind11_FOUND)
    message(WARNING "No pybind11 found. sps python module target will not be available.")
endif(pybind11_FOUND)

#################
# Documentation #
#################

find_package(Doxygen QUIET)
find_package(Sphinx QUIET)
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
     "inc/cm/*.h"
     "src/*.cpp"
)
add_custom_target(format COMMAND "clang-format" "-i" "-style=file" ${format_files}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} VERBATIM)
