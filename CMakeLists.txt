cmake_minimum_required(VERSION 2.6)

# project settings
project(PaSGAL)

##### General Compilation Settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")

#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")

#Set default cmake build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo"
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


#Configure Protobuf (path needs to be configured by each user)
set(PROTOBUF_INCLUDE_DIRS "/projects/PCL-GBB/chirag/software/protobuf/build/include")
set(PROTOBUF_LIBRARY "/projects/PCL-GBB/chirag/software/protobuf/build/lib/libprotobuf.so")

include_directories("${PROJECT_SOURCE_DIR}/ext")
include_directories("${PROJECT_SOURCE_DIR}/src/include")
include_directories(${PROTOBUF_INCLUDE_DIRS})


add_executable(PaSGAL "${PROJECT_SOURCE_DIR}/src/main.cpp" "${PROJECT_SOURCE_DIR}/src/vg.pb.cc")
target_link_libraries(PaSGAL ${PROTOBUF_LIBRARY} -lz -lpthread)

