/**
 * @file    base_types.hpp
 * @brief   basic type formats
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PSGL_BASETYPES_HPP
#define PSGL_BASETYPES_HPP


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
    static const int8_t mismatch = 1;
    static const int8_t ins = 1;
    static const int8_t del = 1;
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
  template <typename ScoreType, typename VertexIdType>
    struct BestScoreInfo
    {
      //positioning in complete DP matrix where optimal alignment ends
      std::size_t refColumn;
      std::size_t qryRow;

      //score value
      ScoreType score;

      char strand; // '+' or '-'

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
