/**
 * @file    test_pasgal_all.cpp
 * @brief   include all tests into a single executable
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "test_graph_load.cpp" 
#include "test_local_alignment.cpp"
#include "test_local_alignment_uniform_len.cpp"

TEST(printEnv, print) 
{
  psgl::showExecutionEnv();
  ASSERT_EQ(1, 1); 
}
