cmake_minimum_required(VERSION 3.1)
include(CheckCXXCompilerFlag)

# project settings
project(PaSGAL)

##### General Compilation Settings
set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} --std=c++11")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

#Set default cmake build type to RelWithDebInfo during prototyping
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo"
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )

# Check if the user want to build google test applications
OPTION(BUILD_TESTS "Inform whether test applications should be built" ON)
# Set vtune support
OPTION(ENABLE_VTUNE "Turn on/off the vtune support with intel compiler" OFF)
# Check if user wants to print verbose debug messages
OPTION(DEFINE_DEBUG_MACRO "Print verbose debug messages" OFF) 

IF(DEFINE_DEBUG_MACRO)
  ADD_DEFINITIONS(-DDEBUG)
ENDIF(DEFINE_DEBUG_MACRO)

# Assume avx512 availability by default
SET(SIMD_SUPPORT avx512 CACHE STRING "Choose SIMD instrinsics, options are: avx512 avx2 none")
SET_PROPERTY(CACHE SIMD_SUPPORT PROPERTY STRINGS avx512 avx2 none)

#User expected to provide the SIMD option
if (SIMD_SUPPORT STREQUAL avx2 OR SIMD_SUPPORT STREQUAL avx512)
  # Check if compiler supports AVX
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORT_AVX)
    if(COMPILER_SUPPORT_AVX)
      message(STATUS "GNU C++ compiler supports AVX")
    else()
      message(FATAL_ERROR "Compiler version required: GCC (>= 4.9) or Intel (>= 15)")
    endif()
  endif()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    check_cxx_compiler_flag("-xCORE-AVX512" COMPILER_SUPPORT_AVX)
    if(COMPILER_SUPPORT_AVX)
      message(STATUS "Intel C++ compiler supports AVX")
    else()
      message(FATAL_ERROR "Compiler version required: GCC (>= 4.9) or  Intel (>= 15)")
    endif()
  endif()

  #compiler flag to enable AVX512 SIMD instructions
  if (SIMD_SUPPORT STREQUAL avx512)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -DPASGAL_ENABLE_AVX512")
    message(STATUS "Enabling AVX512...")
  endif()

  #compiler flag to enable AVX2 simd instructions
  if (SIMD_SUPPORT STREQUAL avx2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -DPASGAL_ENABLE_AVX2")
    message(STATUS "Enabling AVX2...")
  endif()

else()

  message(STATUS "AVX option not provided, disabling SIMD support...")

endif()

# --------------------------------------------------------------
# Indicate CMake 2.7 and above that we don't want to mix relative
#  and absolute paths in linker lib lists.
# Run "cmake --help-policy CMP0003" for more information.
# --------------------------------------------------------------
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)

  #for recognizing numbers and boolean constants in if()
  if(CMAKE_VERSION VERSION_GREATER 2.6.4)
    cmake_policy(SET CMP0012 NEW)
  endif()
endif()

#### OpenMP
find_package(OpenMP)
# FindOpenMP defines the OpenMP_C_FLAGS and OpenMP_CXX_FLAGS.
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message(STATUS "Found OpenMP")
endif()

#compile testing code
if (BUILD_TESTS)
  #Add googleTest support
  add_subdirectory("${PROJECT_SOURCE_DIR}/ext/gtest/googletest")
endif(BUILD_TESTS)

##
#User is expected to set the following variable PROTOBUF_DIR 
##
set(PROTOBUF_DIR CACHE STRING "absolute path to protobuf installation directory")

set(PROTOBUF_INCLUDE_DIR "${PROTOBUF_DIR}/include")
set(PROTOBUF_LIBRARY "${PROTOBUF_DIR}/lib/libprotobuf.so")
set(PROTOBUF_PROTOC_EXECUTABLE "${PROTOBUF_DIR}/bin/protoc")

if(EXISTS "${PROTOBUF_INCLUDE_DIR}" AND EXISTS "${PROTOBUF_LIBRARY}" AND EXISTS "${PROTOBUF_PROTOC_EXECUTABLE}")
  message(STATUS "Protobuf libraries, header files and executables located")
else()
  message(FATAL_ERROR "Protobuf path is invalid, cannot find required libraries, header files or executables. Set the cmake variable PROTOBUF_DIR properly.")
endif()

find_package(Protobuf REQUIRED)

#generate C++ proto files
include_directories(${PROTOBUF_INCLUDE_DIR})
PROTOBUF_GENERATE_CPP(PROTO_SRCS PROTO_HDRS "${PROJECT_SOURCE_DIR}/ext/vg.proto")
add_library(vg_protobuf "vg.pb.cc")

#include external header files
include_directories("${PROJECT_SOURCE_DIR}/ext")

#include our own directories
include_directories("${PROJECT_SOURCE_DIR}/src/include")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory("${PROJECT_SOURCE_DIR}/tests")

#### To compile with Intel vtune library if requested
if(ENABLE_VTUNE)
  #Configure path to vtune_amplifier library 
  include_directories("/opt/intel/vtune_amplifier_2018/include")
  set(VTUNE_LIBRARY "/opt/intel/vtune_amplifier_2018/lib64/libittnotify.a")
  add_definitions(-DVTUNE_SUPPORT)
endif(ENABLE_VTUNE)

add_executable(PaSGAL "${PROJECT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(PaSGAL ${PROTOBUF_LIBRARY} vg_protobuf ${VTUNE_LIBRARY} -lz -lpthread)
