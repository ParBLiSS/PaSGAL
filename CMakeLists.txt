cmake_minimum_required(VERSION 2.6)

# compile using intel compiler 
# users can remove these two lines to compile with gcc
set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_C_COMPILER "icc")

# project settings
project(PaSGAL)

##### General Compilation Settings
set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} --std=c++11 -fopenmp")

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

# Set vtune support
OPTION(ENABLE_VTUNE "Turn on/off the vtune support with intel compiler" OFF)
# Check if the user want to build google test applications
OPTION(BUILD_TESTS "Inform whether test applications should be built" OFF)
# Add AVX512 support
option(BUILD_AVX512 "Support AVX512" ON)

#compiler flag to enable AVX512 simd instructions
if (BUILD_AVX512)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512")
endif (BUILD_AVX512)

#compile testing code
if (BUILD_TESTS)
  #Add googleTest support
  add_subdirectory("${PROJECT_SOURCE_DIR}/ext/googletest")
endif(BUILD_TESTS)

message(STATUS "Warning: Using hard coded paths to locate protobuf library")

############
#Configure PROTOBUF installation directory 
set(PROTOBUF_INSTALL_DIRECTORY "/work/04042/cjain7/stampede2/Intel/utility/protobuf-3.6.0/local_install")
############


set(PROTOBUF_INCLUDE_DIR "${PROTOBUF_INSTALL_DIRECTORY}/include")
set(PROTOBUF_LIBRARY "${PROTOBUF_INSTALL_DIRECTORY}/lib/libprotobuf.so")
set(PROTOBUF_PROTOC_EXECUTABLE "${PROTOBUF_INSTALL_DIRECTORY}/bin/protoc")
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

#### For vtune support
if(ENABLE_VTUNE)
  #Configure path to vtune_amplifier library 
  include_directories("/opt/intel/vtune_amplifier_2018/include")
  set(VTUNE_LIBRARY "/opt/intel/vtune_amplifier_2018/lib64/libittnotify.a")
  add_definitions(-DVTUNE_SUPPORT)
endif(ENABLE_VTUNE)

add_executable(PaSGAL "${PROJECT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(PaSGAL ${PROTOBUF_LIBRARY} vg_protobuf ${VTUNE_LIBRARY} -lz -lpthread)
