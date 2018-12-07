/**
 * @file    test_local_alignment_uniform_len.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "graphLoad.hpp"
#include "align.hpp"
#include "base_types.hpp"
#include "parseCmdArgs.hpp"
#include "googletest/include/gtest/gtest.h"

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define FOLDER STR(PROJECT_TEST_DATA_DIR)

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          16 query sequences to it using single thread.
 *          - All query sequences have equal length
 *          - The count 16 ensures all vector lanes are used
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignmentUniformLen, multipleQuerySequentialUniformLength_vg) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.vg";
  auto qfile = dir + "/BRCA1_16_uniform_len.fastq";

  std::vector<char> QFILE(qfile.c_str(), qfile.c_str() + qfile.size() + 1u);
  std::vector<char> RFILE(rfile.c_str(), rfile.c_str() + rfile.size() + 1u);
  char *mode = "vg";
  char *threads = "1"; 

  char *argv[] = {"PaSGAL", "-m", mode, "-q", QFILE.data(), "-r", RFILE.data(),
                  "-t", threads, "-o", "/dev/null" , nullptr};
  int argc = 11;

  psgl::Parameters parameters;        
  psgl::parseandSave(argc, argv, parameters);

  std::vector< psgl::BestScoreInfo > bestScoreVector;
  psgl::alignToDAG (parameters, psgl::MODE::LOCAL, bestScoreVector);

  //NOTE: Ground truth calculated using unit scoring system

  ASSERT_EQ(bestScoreVector.size(), 16); 

  ASSERT_EQ(bestScoreVector[0].score, 89);       
  ASSERT_EQ(bestScoreVector[0].strand, '-');    

  ASSERT_EQ(bestScoreVector[1].score, 80);       
  ASSERT_EQ(bestScoreVector[1].strand, '+');    

  ASSERT_EQ(bestScoreVector[2].score, 93);       
  ASSERT_EQ(bestScoreVector[2].strand, '-');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '+');    

  ASSERT_EQ(bestScoreVector[4].score, 87);       
  ASSERT_EQ(bestScoreVector[4].strand, '-');    

  ASSERT_EQ(bestScoreVector[5].score, 77);       
  ASSERT_EQ(bestScoreVector[5].strand, '+');    

  ASSERT_EQ(bestScoreVector[6].score, 85);       
  ASSERT_EQ(bestScoreVector[6].strand, '-');    

  ASSERT_EQ(bestScoreVector[7].score, 86);       
  ASSERT_EQ(bestScoreVector[7].strand, '+');    

  ASSERT_EQ(bestScoreVector[8].score, 80);       
  ASSERT_EQ(bestScoreVector[8].strand, '-');    

  ASSERT_EQ(bestScoreVector[9].score, 82);       
  ASSERT_EQ(bestScoreVector[9].strand, '+');    

  ASSERT_EQ(bestScoreVector[10].score, 74);       
  ASSERT_EQ(bestScoreVector[10].strand, '-');    

  ASSERT_EQ(bestScoreVector[11].score, 75);       
  ASSERT_EQ(bestScoreVector[11].strand, '+');    

  ASSERT_EQ(bestScoreVector[12].score, 65);       
  ASSERT_EQ(bestScoreVector[12].strand, '-');    

  ASSERT_EQ(bestScoreVector[13].score, 81);       
  ASSERT_EQ(bestScoreVector[13].strand, '+');    

  ASSERT_EQ(bestScoreVector[14].score, 82);       
  ASSERT_EQ(bestScoreVector[14].strand, '-');    

  ASSERT_EQ(bestScoreVector[15].score, 81);       
  ASSERT_EQ(bestScoreVector[15].strand, '+');    
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          16 query sequences to it using single thread.
 *          - All query sequences have equal length
 *          - The count 16 ensures all vector lanes are used
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignmentUniformLen, multipleQuerySequentialUniformLength_txt) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_16_uniform_len.fastq";

  std::vector<char> QFILE(qfile.c_str(), qfile.c_str() + qfile.size() + 1u);
  std::vector<char> RFILE(rfile.c_str(), rfile.c_str() + rfile.size() + 1u);
  char *mode = "txt";
  char *threads = "1"; 

  char *argv[] = {"PaSGAL", "-m", mode, "-q", QFILE.data(), "-r", RFILE.data(),
                  "-t", threads, "-o", "/dev/null" , nullptr};
  int argc = 11;

  psgl::Parameters parameters;        
  psgl::parseandSave(argc, argv, parameters);

  std::vector< psgl::BestScoreInfo > bestScoreVector;
  psgl::alignToDAG (parameters, psgl::MODE::LOCAL, bestScoreVector);

  //NOTE: Ground truth calculated using unit scoring system

  ASSERT_EQ(bestScoreVector.size(), 16); 

  ASSERT_EQ(bestScoreVector[0].score, 89);       
  ASSERT_EQ(bestScoreVector[0].strand, '+');    

  ASSERT_EQ(bestScoreVector[1].score, 80);       
  ASSERT_EQ(bestScoreVector[1].strand, '-');    

  ASSERT_EQ(bestScoreVector[2].score, 93);       
  ASSERT_EQ(bestScoreVector[2].strand, '+');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '-');    

  ASSERT_EQ(bestScoreVector[4].score, 87);       
  ASSERT_EQ(bestScoreVector[4].strand, '+');    

  ASSERT_EQ(bestScoreVector[5].score, 77);       
  ASSERT_EQ(bestScoreVector[5].strand, '-');    

  ASSERT_EQ(bestScoreVector[6].score, 85);       
  ASSERT_EQ(bestScoreVector[6].strand, '+');    

  ASSERT_EQ(bestScoreVector[7].score, 86);       
  ASSERT_EQ(bestScoreVector[7].strand, '-');    

  ASSERT_EQ(bestScoreVector[8].score, 80);       
  ASSERT_EQ(bestScoreVector[8].strand, '+');    

  ASSERT_EQ(bestScoreVector[9].score, 82);       
  ASSERT_EQ(bestScoreVector[9].strand, '-');    

  ASSERT_EQ(bestScoreVector[10].score, 74);       
  ASSERT_EQ(bestScoreVector[10].strand, '+');    

  ASSERT_EQ(bestScoreVector[11].score, 75);       
  ASSERT_EQ(bestScoreVector[11].strand, '-');    

  ASSERT_EQ(bestScoreVector[12].score, 65);       
  ASSERT_EQ(bestScoreVector[12].strand, '+');    

  ASSERT_EQ(bestScoreVector[13].score, 81);       
  ASSERT_EQ(bestScoreVector[13].strand, '-');    

  ASSERT_EQ(bestScoreVector[14].score, 82);       
  ASSERT_EQ(bestScoreVector[14].strand, '+');    

  ASSERT_EQ(bestScoreVector[15].score, 81);       
  ASSERT_EQ(bestScoreVector[15].strand, '-');    
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          16 query sequences to it using multiple threads.
 *          - All query sequences have equal length
 *          - The count 16 ensures all vector lanes are used
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignmentUniformLen, multipleQueryParallelUniformLength_vg) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.vg";
  auto qfile = dir + "/BRCA1_16_uniform_len.fastq";

  std::vector<char> QFILE(qfile.c_str(), qfile.c_str() + qfile.size() + 1u);
  std::vector<char> RFILE(rfile.c_str(), rfile.c_str() + rfile.size() + 1u);
  char *mode = "vg";
  char *threads = "8"; 

  char *argv[] = {"PaSGAL", "-m", mode, "-q", QFILE.data(), "-r", RFILE.data(),
                  "-t", threads, "-o", "/dev/null" , nullptr};
  int argc = 11;

  psgl::Parameters parameters;        
  psgl::parseandSave(argc, argv, parameters);
  
  std::vector< psgl::BestScoreInfo > bestScoreVector;
  psgl::alignToDAG (parameters, psgl::MODE::LOCAL, bestScoreVector);

  //NOTE: Ground truth calculated using unit scoring system

  ASSERT_EQ(bestScoreVector.size(), 16); 

  ASSERT_EQ(bestScoreVector[0].score, 89);       
  ASSERT_EQ(bestScoreVector[0].strand, '-');    

  ASSERT_EQ(bestScoreVector[1].score, 80);       
  ASSERT_EQ(bestScoreVector[1].strand, '+');    

  ASSERT_EQ(bestScoreVector[2].score, 93);       
  ASSERT_EQ(bestScoreVector[2].strand, '-');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '+');    

  ASSERT_EQ(bestScoreVector[4].score, 87);       
  ASSERT_EQ(bestScoreVector[4].strand, '-');    

  ASSERT_EQ(bestScoreVector[5].score, 77);       
  ASSERT_EQ(bestScoreVector[5].strand, '+');    

  ASSERT_EQ(bestScoreVector[6].score, 85);       
  ASSERT_EQ(bestScoreVector[6].strand, '-');    

  ASSERT_EQ(bestScoreVector[7].score, 86);       
  ASSERT_EQ(bestScoreVector[7].strand, '+');    

  ASSERT_EQ(bestScoreVector[8].score, 80);       
  ASSERT_EQ(bestScoreVector[8].strand, '-');    

  ASSERT_EQ(bestScoreVector[9].score, 82);       
  ASSERT_EQ(bestScoreVector[9].strand, '+');    

  ASSERT_EQ(bestScoreVector[10].score, 74);       
  ASSERT_EQ(bestScoreVector[10].strand, '-');    

  ASSERT_EQ(bestScoreVector[11].score, 75);       
  ASSERT_EQ(bestScoreVector[11].strand, '+');    

  ASSERT_EQ(bestScoreVector[12].score, 65);       
  ASSERT_EQ(bestScoreVector[12].strand, '-');    

  ASSERT_EQ(bestScoreVector[13].score, 81);       
  ASSERT_EQ(bestScoreVector[13].strand, '+');    

  ASSERT_EQ(bestScoreVector[14].score, 82);       
  ASSERT_EQ(bestScoreVector[14].strand, '-');    

  ASSERT_EQ(bestScoreVector[15].score, 81);       
  ASSERT_EQ(bestScoreVector[15].strand, '+'); 
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          16 query sequences to it using multiple threads.
 *          - All query sequences have equal length
 *          - The count 16 ensures all vector lanes are used
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignmentUniformLen, multipleQueryParallelUniformLength_txt) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_16_uniform_len.fastq";

  std::vector<char> QFILE(qfile.c_str(), qfile.c_str() + qfile.size() + 1u);
  std::vector<char> RFILE(rfile.c_str(), rfile.c_str() + rfile.size() + 1u);
  char *mode = "txt";
  char *threads = "8"; 

  char *argv[] = {"PaSGAL", "-m", mode, "-q", QFILE.data(), "-r", RFILE.data(),
                  "-t", threads, "-o", "/dev/null" , nullptr};
  int argc = 11;

  psgl::Parameters parameters;        
  psgl::parseandSave(argc, argv, parameters);

  std::vector< psgl::BestScoreInfo > bestScoreVector;
  psgl::alignToDAG (parameters, psgl::MODE::LOCAL, bestScoreVector);

  //NOTE: Ground truth calculated using unit scoring system
  ASSERT_EQ(bestScoreVector.size(), 16); 

  ASSERT_EQ(bestScoreVector[0].score, 89);       
  ASSERT_EQ(bestScoreVector[0].strand, '+');    

  ASSERT_EQ(bestScoreVector[1].score, 80);       
  ASSERT_EQ(bestScoreVector[1].strand, '-');    

  ASSERT_EQ(bestScoreVector[2].score, 93);       
  ASSERT_EQ(bestScoreVector[2].strand, '+');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '-');    

  ASSERT_EQ(bestScoreVector[4].score, 87);       
  ASSERT_EQ(bestScoreVector[4].strand, '+');    

  ASSERT_EQ(bestScoreVector[5].score, 77);       
  ASSERT_EQ(bestScoreVector[5].strand, '-');    

  ASSERT_EQ(bestScoreVector[6].score, 85);       
  ASSERT_EQ(bestScoreVector[6].strand, '+');    

  ASSERT_EQ(bestScoreVector[7].score, 86);       
  ASSERT_EQ(bestScoreVector[7].strand, '-');    

  ASSERT_EQ(bestScoreVector[8].score, 80);       
  ASSERT_EQ(bestScoreVector[8].strand, '+');    

  ASSERT_EQ(bestScoreVector[9].score, 82);       
  ASSERT_EQ(bestScoreVector[9].strand, '-');    

  ASSERT_EQ(bestScoreVector[10].score, 74);       
  ASSERT_EQ(bestScoreVector[10].strand, '+');    

  ASSERT_EQ(bestScoreVector[11].score, 75);       
  ASSERT_EQ(bestScoreVector[11].strand, '-');    

  ASSERT_EQ(bestScoreVector[12].score, 65);       
  ASSERT_EQ(bestScoreVector[12].strand, '+');    

  ASSERT_EQ(bestScoreVector[13].score, 81);       
  ASSERT_EQ(bestScoreVector[13].strand, '-');    

  ASSERT_EQ(bestScoreVector[14].score, 82);       
  ASSERT_EQ(bestScoreVector[14].strand, '+');    

  ASSERT_EQ(bestScoreVector[15].score, 81);       
  ASSERT_EQ(bestScoreVector[15].strand, '-');    
}
