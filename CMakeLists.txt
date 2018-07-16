cmake_minimum_required(VERSION 2.6)

# compile using intel compiler 
# users can remove these two lines to compile with gcc
set(CMAKE_CXX_COMPILER "icc")
set(CMAKE_C_COMPILER "icc")

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

# Set vtune support
OPTION(ENABLE_VTUNE "Turn on/off the vtune support with intel compiler" OFF)
# Check if the user want to build google test applications
OPTION(BUILD_TESTS "Inform whether test applications should be built" OFF)
# Add AVX512 support
option(BUILD_AVX512 "Support AVX512" ON)


if (BUILD_AVX512)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512")
endif (BUILD_AVX512)


#Add OpenMP support
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if (BUILD_TESTS)
  #Add googleTest support
  add_subdirectory("${PROJECT_SOURCE_DIR}/ext/googletest")
endif(BUILD_TESTS)


#Configure Protobuf (path needs to be configured by each user)
#TODO: automate generation of protobuf cpp file
set(PROTOBUF_INCLUDE_DIRS "/data/nfs_home/nasamo_data_backup/mpatwary_copy/current/ib01_bakcup/minconda_install_OLD/pkgs/libprotobuf-3.0.0-0/include")
set(PROTOBUF_LIBRARY "/data/nfs_home/nasamo_data_backup/mpatwary_copy/current/ib01_bakcup/minconda_install_OLD/pkgs/libprotobuf-3.0.0-0/lib/libprotobuf.so")
include_directories(${PROTOBUF_INCLUDE_DIRS})
add_library(vg_protobuf "${PROJECT_SOURCE_DIR}/ext/vg.pb.cc")

#include external header files
include_directories("${PROJECT_SOURCE_DIR}/ext")

#include our own directories
include_directories("${PROJECT_SOURCE_DIR}/src/include")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests")

#### For vtune support
if(ENABLE_VTUNE)
  include_directories("/opt/intel/vtune_amplifier_2018/include")
  set(VTUNE_LIBRARY "/opt/intel/vtune_amplifier_2018/lib64/libittnotify.a")
  add_definitions(-DVTUNE_SUPPORT)
endif(ENABLE_VTUNE)

add_executable(PaSGAL "${PROJECT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(PaSGAL ${PROTOBUF_LIBRARY} vg_protobuf ${VTUNE_LIBRARY} -lz -lpthread)
