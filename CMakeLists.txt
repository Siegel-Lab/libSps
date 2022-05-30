cmake_minimum_required (VERSION 3.9)

project("libSps" VERSION 0.0.1)

# set building python library as default
option(WITH_PYTHON "Build python library." ON)

set(NUM_DIMENSIONS_A "2" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_B "3" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_C "0" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")
set(NUM_DIMENSIONS_D "0" CACHE STRING "Build datastructures with this amount of dimensions (put zero to disable).")

option(W_DEPENDENT_DIM "Build datastructures with dependent dimension." OFF)
option(WO_DEPENDENT_DIM "Build datastructures without dependent dimension." ON)

option(W_CUBES "" OFF)
option(W_RECTANGLES "" OFF)
option(W_INTERVALS "" OFF)
option(W_POINTS "" ON)

option(DISK "Build disk based datastructures (IO on load & unload)." OFF)
option(CACHED "Build cached datastructures (IO during runtime)." ON)
option(RAM "Build ram datastructures (Non-persistent)." ON)

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
    add_subdirectory (contrib/pybind11 EXCLUDE_FROM_ALL)
endif()
add_subdirectory (contrib/stxxl EXCLUDE_FROM_ALL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# create python module if python has been found
if(WITH_PYTHON)
    pybind11_add_module( Sps src/tree.cpp )
    target_include_directories( Sps PUBLIC inc )
    target_include_directories( Sps SYSTEM PRIVATE ${STXXL_INCLUDE_DIRS} )

    target_compile_definitions( Sps PRIVATE WITH_PYTHON )

    target_compile_definitions( Sps PRIVATE NUM_DIMENSIONS_A=${NUM_DIMENSIONS_A} )
    target_compile_definitions( Sps PRIVATE NUM_DIMENSIONS_B=${NUM_DIMENSIONS_B} )
    target_compile_definitions( Sps PRIVATE NUM_DIMENSIONS_C=${NUM_DIMENSIONS_C} )
    target_compile_definitions( Sps PRIVATE NUM_DIMENSIONS_D=${NUM_DIMENSIONS_D} )

    if(${W_DEPENDENT_DIM})
        target_compile_definitions( Sps PRIVATE W_DEPENDANT_DIM )
    endif()
    if(${WO_DEPENDENT_DIM})
        target_compile_definitions( Sps PRIVATE WO_DEPENDANT_DIM )
    endif()
    
    if(${DISK})
        target_compile_definitions( Sps PRIVATE DISK )
    endif()
    if(${CACHED})
        target_compile_definitions( Sps PRIVATE CACHED )
    endif()
    if(${RAM})
        target_compile_definitions( Sps PRIVATE RAM )
    endif()

    if(${W_CUBES})
        target_compile_definitions( Sps PRIVATE W_CUBES )
    endif()
    if(${W_RECTANGLES})
        target_compile_definitions( Sps PRIVATE W_RECTANGLES )
    endif()
    if(${W_INTERVALS})
        target_compile_definitions( Sps PRIVATE W_INTERVALS )
    endif()
    if(${W_POINTS})
        target_compile_definitions( Sps PRIVATE W_POINTS )
    endif()
    target_link_libraries( Sps PUBLIC stxxl )
    #target_link_libraries( Sps PUBLIC -lpthread )
endif()

# build a cpp only test project (this wll show any linker errors that might be hidden with the python build)
add_executable ( test src/main.cpp )
target_include_directories( test PUBLIC inc )
target_include_directories( test SYSTEM PUBLIC ${STXXL_INCLUDE_DIRS} )
target_link_libraries(test PUBLIC stxxl )
#target_link_libraries(test PUBLIC -lpthread )

# Documentation setup


set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) 
find_package(Doxygen)
find_package(Sphinx)
if (DOXYGEN_FOUND)

if(Sphinx_FOUND)

    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs/")

    add_custom_target(copy_docs_conf ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory 
        "${CMAKE_CURRENT_SOURCE_DIR}/docs_conf" "${CMAKE_CURRENT_BINARY_DIR}/_docs_build")

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs_conf/doxygen.config)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config.out)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( docs-cpp ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS copy_docs_conf
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

    add_custom_target(docs-py ALL
        ${SPHINX_EXECUTABLE}
            -E -a -q -b html
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.libSps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_HTML_DIR}"
        DEPENDS Sps docs-cpp copy_docs_conf
        COMMENT "Building HTML documentation with Sphinx")

    add_custom_target(docs-py-latex ALL
        ${SPHINX_EXECUTABLE}
            -E -a -q -b latex
            -d "${SPHINX_CACHE_DIR}"
            -Dbreathe_projects.libSps="${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml"
            "${BINARY_BUILD_DIR}"
            "${SPHINX_TEX_DIR}"
        DEPENDS Sps docs-cpp copy_docs_conf
        COMMENT "Building tex documentation with Sphinx") # tectonic *.tex

else(Sphinx_FOUND)
    message(WARNING "No Sphinx found. Documentation target not available.")
endif(Sphinx_FOUND)

else (DOXYGEN_FOUND)
    message("No Doxygen foud. Documentation target not available.")
endif (DOXYGEN_FOUND)


