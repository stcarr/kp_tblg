# Specify the minimum version for CMake

cmake_minimum_required(VERSION 2.4.0)

# Project's name
PROJECT(tblg_kp)
SET(CMAKE_BUILD_TYPE Debug)

OPTION(PYBIND_MAKE "Uses pybind11 to make a .so that can be imported to python" ON)

# Define name and version of package
SET(PACKAGE_NAME "tblg_kp" CACHE STRING "package_name" FORCE)
SET(PACKAGE_VERSION "1.0" CACHE STRING "package_version" FORCE)

# set code language
ENABLE_LANGUAGE(CXX)
#SET(CMAKE_CXX_COMPILER /usr/bin/c++)

# setup packages

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/module/")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-variable -Wno-unused-parameter -fPIC")
find_package(Eigen3 REQUIRED)


#MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )

if(PYBIND_MAKE)
	cmake_minimum_required(VERSION 2.8.12)
	add_subdirectory(externals/pybind11)
	pybind11_add_module(tblg_kpy pysrc/main.cpp)
	target_link_libraries(tblg_kpy PRIVATE lib_kp_tblg)
endif()

# include required directories
include_directories(${EIGEN3_INCLUDE_DIR})

include(CMakeBackwardCompatibilityCXX)
include(CheckLibraryExists)

# link necessary subdirectory library paths

add_subdirectory(src)
add_subdirectory(app)
