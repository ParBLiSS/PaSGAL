/**
 * @file    test_local_alignment.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "graphLoad.hpp"
#include "align.hpp"
#include "base_types.hpp"
#include "parseCmdArgs.hpp"
#include "gtest/googletest/include/gtest/gtest.h"

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define FOLDER STR(PROJECT_TEST_DATA_DIR)

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          single query sequence to it.
 *          This routine checks for alignment strand,
 *          score and cigar
 **/
TEST(localAlignment, singleQuerySequential_vg) 
{
  //get file name
  std::string dir = FOLDER;
  std::string rfile = dir + "/BRCA1_seq_graph.vg";
  std::string qfile = dir + "/BRCA1_1_read.fastq";

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

  ASSERT_EQ(bestScoreVector.size(), 1); 
  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '-');       
  ASSERT_EQ(psgl::seqUtils::cigarScore (bestScoreVector[0].cigar, parameters), 482);
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          single query sequence to it.
 *          This routine checks for alignment strand,
 *          score and cigar
 **/
TEST(localAlignment, singleQuerySequential_txt) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_1_read.fastq";

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

  ASSERT_EQ(bestScoreVector.size(), 1); 
  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '+');       
  ASSERT_EQ(psgl::seqUtils::cigarScore (bestScoreVector[0].cigar, parameters), 482);
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          5 query sequences to it in parallel.
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignment, multipleQueryParallelScore_vg) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.vg";
  auto qfile = dir + "/BRCA1_5_reads.fastq";

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

  ASSERT_EQ(bestScoreVector.size(), 5); 

  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '-');    

  ASSERT_EQ(bestScoreVector[1].score, 122);       
  ASSERT_EQ(bestScoreVector[1].strand, '+');    

  ASSERT_EQ(bestScoreVector[2].score, 441);       
  ASSERT_EQ(bestScoreVector[2].strand, '-');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '+');    

  ASSERT_EQ(bestScoreVector[4].score, 259);       
  ASSERT_EQ(bestScoreVector[4].strand, '-');    
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          5 query sequences to it in parallel.
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignment, multipleQueryParallelScore_txt) 
{
  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_5_reads.fastq";

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

  ASSERT_EQ(bestScoreVector.size(), 5); 

  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '+');    

  ASSERT_EQ(bestScoreVector[1].score, 122);       
  ASSERT_EQ(bestScoreVector[1].strand, '-');    

  ASSERT_EQ(bestScoreVector[2].score, 441);       
  ASSERT_EQ(bestScoreVector[2].strand, '+');    

  ASSERT_EQ(bestScoreVector[3].score, 90);       
  ASSERT_EQ(bestScoreVector[3].strand, '-');    

  ASSERT_EQ(bestScoreVector[4].score, 259);       
  ASSERT_EQ(bestScoreVector[4].strand, '+');    
}
