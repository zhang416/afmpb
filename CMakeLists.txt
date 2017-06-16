cmake_minimum_required (VERSION 3.4)

project (afmpb) 

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "-O3") 

find_package (PkgConfig)
pkg_search_module (HPX REQUIRED hpx)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/example)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(contrib)

