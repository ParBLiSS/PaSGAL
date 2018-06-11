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
   */
  struct SCORE
  {
    static const int match = 1;

    //penalties
    static const int mismatch = 1;
    static const int ins = 1;
    static const int del = 1;
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

}

#endif
