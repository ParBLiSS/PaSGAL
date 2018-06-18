/**
 * @file    utils.hpp
 * @brief   functions for common use cases
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PSGL_UTILS_HPP
#define PSGL_UTILS_HPP

#include  <random>
#include  <iterator>

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
  }

}

#endif
