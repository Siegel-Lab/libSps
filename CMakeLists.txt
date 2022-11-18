cmake_minimum_required (VERSION 3.9)

# only adjust major and minor version - patch is taken from git commit number
project("libSps" VERSION 0.2.0)

# set building python library as default
option(WITH_PYTHON "Build python library." ON)

set(NUM_DIMENSIONS_A "2" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_B "3" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_C "0" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_D "0" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")

option(W_CUBES "Build datastructures with cubes as data-elements" OFF)
option(W_RECTANGLES "Build datastructures with rectangles as data-elements" OFF)
option(W_INTERVALS "Build datastructures with intervals as data-elements" OFF)
option(W_POINTS "Build datastructures with points as data-elements" ON)

option(DISK "Build disk based datastructures (IO on load & unload)." OFF)
option(CACHED "Build cached datastructures (IO during runtime)." ON)
option(RAM "Build ram datastructures (Non-persistent)." ON)

option(DO_PROFILE "Add profileing code to the llibrary. Uses gperftools." OFF)

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

    target_compile_definitions( libSps PRIVATE NUM_DIMENSIONS_A=${NUM_DIMENSIONS_A} )
    target_compile_definitions( libSps PRIVATE NUM_DIMENSIONS_B=${NUM_DIMENSIONS_B} )
    target_compile_definitions( libSps PRIVATE NUM_DIMENSIONS_C=${NUM_DIMENSIONS_C} )
    target_compile_definitions( libSps PRIVATE NUM_DIMENSIONS_D=${NUM_DIMENSIONS_D} )

    if(${DISK})
        target_compile_definitions( libSps PRIVATE DISK )
    endif()
    if(${CACHED})
        target_compile_definitions( libSps PRIVATE CACHED )
    endif()
    if(${RAM})
        target_compile_definitions( libSps PRIVATE RAM )
    endif()

    if(${W_CUBES})
        target_compile_definitions( libSps PRIVATE W_CUBES )
    endif()
    if(${W_RECTANGLES})
        target_compile_definitions( libSps PRIVATE W_RECTANGLES )
    endif()
    if(${W_INTERVALS})
        target_compile_definitions( libSps PRIVATE W_INTERVALS )
    endif()
    if(${W_POINTS})
        target_compile_definitions( libSps PRIVATE W_POINTS )
    endif()
    target_link_libraries( libSps PUBLIC stxxl )
    target_include_directories( libSps PUBLIC ${CMAKE_BINARY_DIR}/generated/inc)
    add_dependencies(libSps sps_config_version sps_config_build_time)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) 
if(DO_PROFILE)
    # build a cpp only test project (this wll show any linker errors that might be hidden with the python build)
    #find_package(Gperftools REQUIRED)
    add_executable ( test src/main.cpp )
    target_include_directories( test PUBLIC inc )
    target_include_directories( test SYSTEM PUBLIC ${STXXL_INCLUDE_DIRS} )
    target_link_libraries(test PUBLIC stxxl )
    target_include_directories( test PUBLIC ${CMAKE_BINARY_DIR}/generated/inc)
    add_dependencies(test sps_config_version sps_config_build_time)

    # build benchmakring
    add_executable ( benchmark src/benchmark.cpp )

    target_compile_definitions( benchmark PRIVATE NUM_DIMENSIONS_A=2 )
    target_compile_definitions( benchmark PRIVATE NUM_DIMENSIONS_B=0 )
    target_compile_definitions( benchmark PRIVATE NUM_DIMENSIONS_C=0 )
    target_compile_definitions( benchmark PRIVATE NUM_DIMENSIONS_D=0 )
    target_compile_definitions( benchmark PRIVATE RAM )
    target_compile_definitions( benchmark PRIVATE W_POINTS )
    target_compile_definitions( benchmark PRIVATE DO_PROFILE )

    target_include_directories( benchmark PUBLIC /usr/include )
    target_include_directories( benchmark PUBLIC inc )
    target_include_directories( benchmark SYSTEM PUBLIC ${STXXL_INCLUDE_DIRS} )

    target_link_libraries(benchmark PUBLIC stxxl )
    target_link_libraries( benchmark PUBLIC profiler )

    add_dependencies(benchmark sps_config_version sps_config_build_time)
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
        DEPENDS Sps docs-doxygen
        COMMENT "Building API documentation with Sphinx")

else(Sphinx_FOUND)
    message(WARNING "No Sphinx found. Documentation target not available.")
endif(Sphinx_FOUND)

else (DOXYGEN_FOUND)
    message("No Doxygen foud. Documentation target not available.")
endif (DOXYGEN_FOUND)


