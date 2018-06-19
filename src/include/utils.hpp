/**
 * @file    utils.hpp
 * @brief   functions for common use cases
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PSGL_UTILS_HPP
#define PSGL_UTILS_HPP

#include <random>
#include <iterator>
#include <algorithm>

#include "base_types.hpp"

namespace psgl
{
  namespace random
  {

    /**
     * @brief       select a random element from C++ container 
     *              given a random generator
     */
    template<typename Iter, typename RandomGenerator>
      Iter select(Iter start, Iter end, RandomGenerator& g) 
      {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(g));
        return start;
      }

    /**
     * @brief       select a random element from C++ container
     * @details     borrowed from stackoverflow answer by Christopher Smith
     *              https://stackoverflow.com/questions/6942273
     */
    template<typename Iter>
      Iter select(Iter start, Iter end) 
      {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return select(start, end, gen);
      }
  }

  namespace seqUtils
  {
    /**
     * @brief                     string compaction, e.g. replace ..MMMMM.. with ..5M..
     * @param[in/out]     cigar   the cigar string
     */
    void cigarCompact (std::string &cigar)
    {
      std::string cigarCpy = cigar;

      cigar.clear();

      for (auto it = cigarCpy.begin(); it != cigarCpy.end(); )
      {
        char current = *it;

        auto it2 = std::find_if (it, cigarCpy.end(), [&](const char &c) { return c != current; });

        std::size_t numeric = std::distance (it , it2);
        cigar.append (std::to_string( numeric ));
        cigar.push_back ( *it );

        it = it2;
      }
    }

    /**
     * @brief                 compute alignment score from cigar string
     * @param[in]     cigar   the cigar string
     */
    template <typename ScoreType>
      ScoreType cigarScore (const std::string &cigar)
      {
        if (cigar.empty()) 
          return 0;

        ScoreType score = 0;
        int currentNumeric = 0;

        for(int i = 0; i < cigar.length(); i++)
        {
          char c = cigar.at(i);

          if ( isdigit(c) )
          {
            currentNumeric = (currentNumeric * 10) + (c - '0');
          }
          else 
          {
            assert (c == '=' || c == 'X' || c == 'I' || c == 'D');

            if ( c == '=' )
              score += SCORE::match * currentNumeric;
            else if ( c == 'X')
              score -= SCORE::mismatch * currentNumeric;
            else if (c == 'I')
              score -= SCORE::ins * currentNumeric;
            else  // c == 'D'
              score -= SCORE::del * currentNumeric;

            currentNumeric = 0;
          }
        }

        return score;
      }
  }
}

#endif
