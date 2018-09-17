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
   * @brief     scores
   * @details   assume edits wrt. reference sequence
   */
  struct SCORE
  {
    static const int8_t match = 1;

    //penalties
    static const int8_t mismatch = -1;
    static const int8_t ins = -1;
    static const int8_t del = -1;
  };

  /**
   * @brief     alignment modes
   */
  enum MODE
  {
    GLOBAL,     
    LOCAL,
    SEMIGLOBAL
  };  

  /**
   * @brief                   container to save info about best score
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   */
  template <typename ScoreType>
    struct BestScoreInfo
    {
      //coordinates in complete DP matrix where optimal alignment begins and ends (both inclusive)
      //these are 0-based offsets
      int32_t refColumnStart;
      int32_t refColumnEnd;

      int32_t qryRowStart;
      int32_t qryRowEnd;

      //score value
      ScoreType score;

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
