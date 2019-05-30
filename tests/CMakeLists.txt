cmake_minimum_required(VERSION 2.6)

# project settings
project(PaSGAL-test)

##### add definitions to find path to test data files
add_definitions(-DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/data)

if (BUILD_TESTS)

  add_executable(test-graph_load test_graph_load.cpp)
  target_link_libraries(test-graph_load gtest_main ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY} ${VTUNE_LIBRARY} -lz -lpthread)

  add_executable(test-local_alignment test_local_alignment.cpp)
  target_link_libraries(test-local_alignment gtest_main ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY} ${VTUNE_LIBRARY} -lz -lpthread)

  add_executable(test-local_alignment_uniform_len test_local_alignment_uniform_len.cpp)
  target_link_libraries(test-local_alignment_uniform_len gtest_main ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY} ${VTUNE_LIBRARY} -lz -lpthread)

  add_executable(test-pasgal_all test_pasgal_all.cpp)
  target_link_libraries(test-pasgal_all gtest_main ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY} ${VTUNE_LIBRARY} -lz -lpthread)

endif(BUILD_TESTS)
