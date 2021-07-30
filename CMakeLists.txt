cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(stencil LANGUAGES CXX VERSION 0.1.0.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# only static libraries are suported by the device linker?
# so I think we need STATIC with CUDA_SEPARABLE_COMPILATION
add_library(mm INTERFACE
include/mm/mm.hpp)

# require c++11
target_compile_features(mm INTERFACE cxx_std_11)

# "this command should be in the source directory root for CTest to find the test file"
enable_testing() 

add_subdirectory(test)
