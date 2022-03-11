cmake_minimum_required (VERSION 3.8)

project("KdpsTree" VERSION 0.0.1)

# disallow in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory (contrib/pybind11 EXCLUDE_FROM_ALL)
add_subdirectory (contrib/stxxl EXCLUDE_FROM_ALL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

pybind11_add_module( KdpsTree src/tree.cpp )
target_include_directories( KdpsTree PUBLIC inc )
target_include_directories( KdpsTree PUBLIC SYSTEM  ${STXXL_INCLUDE_DIRS} )
target_compile_definitions( KdpsTree PRIVATE WITH_PYTHON )
target_link_libraries(KdpsTree PUBLIC -lstdc++ ${STXXL_LIBRARIES} )

#add_executable ( test src/main.cpp )
#target_include_directories( test PUBLIC inc )
#target_include_directories( test PUBLIC SYSTEM ${STXXL_INCLUDE_DIRS} )
#link_libraries(test PUBLIC -lstdc++ ${STXXL_LIBRARIES} )
