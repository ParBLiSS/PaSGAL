/**
 * @file    test_local_alignment.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "graphLoad.hpp"
#include "align.hpp"
#include "base_types.hpp"
#include "googletest/include/gtest/gtest.h"

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
  //Run sequentially
  omp_set_num_threads( 1 ); 

  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.vg";
  auto qfile = dir + "/BRCA1_1_read.fastq";

  using ScoreType = int32_t;

  //load graph

  psgl::graphLoader g;
  g.loadFromVG(rfile);

  std::vector< psgl::BestScoreInfo<ScoreType> > bestScoreVector;
  psgl::alignToDAG<ScoreType> (qfile, g.diCharGraph, bestScoreVector, psgl::MODE::LOCAL);  

  //NOTE: Ground truth calculated using unit scoring system

  ASSERT_EQ(bestScoreVector.size(), 1); 
  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '-');       
  ASSERT_EQ(psgl::seqUtils::cigarScore<ScoreType> (bestScoreVector[0].cigar), 482);
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          single query sequence to it.
 *          This routine checks for alignment strand,
 *          score and cigar
 **/
TEST(localAlignment, singleQuerySequential_txt) 
{
  //Run sequentially
  omp_set_num_threads( 1 ); 

  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_1_read.fastq";

  using ScoreType = int32_t;

  //load graph

  psgl::graphLoader g;
  g.loadFromTxt(rfile);

  std::vector< psgl::BestScoreInfo<ScoreType> > bestScoreVector;
  psgl::alignToDAG<ScoreType> (qfile, g.diCharGraph, bestScoreVector, psgl::MODE::LOCAL);  

  //NOTE: Ground truth calculated using unit scoring system

  ASSERT_EQ(bestScoreVector.size(), 1); 
  ASSERT_EQ(bestScoreVector[0].score, 482);       
  ASSERT_EQ(bestScoreVector[0].strand, '+');       
  ASSERT_EQ(psgl::seqUtils::cigarScore<ScoreType> (bestScoreVector[0].cigar), 482);
}

/**
 * @brief   builds a graph from BRCA1 sequence and aligns
 *          5 query sequences to it in parallel.
 *          This routine checks for alignment strands and 
 *          scores
 **/
TEST(localAlignment, multipleQueryParallelScore_vg) 
{
  //Run in parallel
  omp_set_num_threads( 8 ); 

  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.vg";
  auto qfile = dir + "/BRCA1_5_reads.fastq";

  using ScoreType = int32_t;

  //load graph

  psgl::graphLoader g;
  g.loadFromVG(rfile);

  std::vector< psgl::BestScoreInfo<ScoreType> > bestScoreVector;
  psgl::alignToDAG<ScoreType> (qfile, g.diCharGraph, bestScoreVector, psgl::MODE::LOCAL);  

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
  //Run in parallel
  omp_set_num_threads( 8 ); 

  //get file name
  std::string dir = FOLDER;
  auto rfile = dir + "/BRCA1_seq_graph.txt";
  auto qfile = dir + "/BRCA1_5_reads.fastq";

  using ScoreType = int32_t;

  //load graph

  psgl::graphLoader g;
  g.loadFromTxt(rfile);

  std::vector< psgl::BestScoreInfo<ScoreType> > bestScoreVector;
  psgl::alignToDAG<ScoreType> (qfile, g.diCharGraph, bestScoreVector, psgl::MODE::LOCAL);  

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
