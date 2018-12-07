/**
 * @file    base_types.hpp
 * @brief   basic type formats
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PSGL_BASETYPES_HPP
#define PSGL_BASETYPES_HPP

#include <immintrin.h>

#define psgl_max(a,b) (((a)>(b))?(a):(b))
#define ASSUMED_CPU_FREQ 2100000000
#define PSGL_STATUS_OK 0

namespace psgl
{

  /**
   * @brief     input parameters that are expected 
   *            as command line arguments
   **/
  struct Parameters
  {
    std::string rfile;        //reference graph file
    std::string mode;         //reference graph format
    std::string qfile;        //query sequence file
    std::string ofile;        //output file
    int threads;              //thread count

    int match;                //match score 
    int mismatch;             //mismatch penalty (abs. value) 
    int ins;                  //insertion penalty (abs. value) 
    int del;                  //deletion penalty (abs. value)
  };

  /**
   * @brief     alignment modes
   */
  enum MODE
  {
    GLOBAL,     //TODO
    LOCAL,
    SEMIGLOBAL  //TODO
  };  

  /**
   * @brief                   container to save info about best score
   */
  struct BestScoreInfo
  {
    //coordinates in complete DP matrix where optimal alignment begins and ends (both inclusive)
    //these are 0-based offsets
    int32_t refColumnStart;
    int32_t refColumnEnd;

    int32_t qryRowStart;
    int32_t qryRowEnd;

    //score value
    int32_t score;

    char strand; // '+' or '-'

    //TODO: Storing cigar may be expensive, consider removing later
    std::string cigar;

    /**
     * @brief   constructor
     */
    BestScoreInfo()
    {
      this->score = 0;
    }
  };
}

#endif
