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
#include <thread>
#include <chrono>
#include <fstream>
#include <omp.h>

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
     * @brief   reverse complement of string
     * @note    assumes dest is pre-allocated
     */
    void reverseComplement(const std::string &src, std::string &dest) 
    {
      assert(src.length() == dest.length());

      for ( int i = 0; i < src.length(); i++ )
      {    
        char base = src.at(i);

        assert(std::isupper(base));

        switch ( base )
        {    
          case 'A': base = 'T'; break;
          case 'C': base = 'G'; break;
          case 'G': base = 'C'; break;
          case 'T': base = 'A'; break;
          default: break;
        }    

        dest[src.length() - i - 1] = base;
      }    
    }

    /**
     * @brief               convert DNA or AA alphabets to upper case
     * @param[in]   seq     pointer to input sequence
     * @param[in]   len     length of input sequence
     */
    void makeUpperCase(char* seq, std::size_t len)
    {
      for ( int i = 0; i < len; i++ )
      {
        if (seq[i] > 96 && seq[i] < 123)
        {
          seq[i] -= 32;
        }
      }
    }

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

  namespace timer
  {
    /**
     * @brief  get CPU cycle count
     */
    uint64_t rdtsc()
    {
      return __rdtsc();
    }

    /**
     * @brief   return count of CPU cycles per second
     */
    uint64_t cycles_per_sec()
    {
      uint64_t tick1 = rdtsc();

      //sleep for a second
      std::this_thread::sleep_for (std::chrono::seconds(1));

      uint64_t tick2 = rdtsc();

      return tick2 - tick1;
    }
  }

  /**
   * @brief   print thread count
   */
  void printThreadCount()
  {
#pragma omp parallel
    {
      int tid = omp_get_thread_num();

      if (tid == 0) 
        std::cout << "INFO, psgl::printThreadCount, Number of openmp threads available = " << omp_get_num_threads() << std::endl;
    }
  }

  /**
   * @brief     check if file is accessible
   */
  bool fileExists(const std::string &filename)
  {
    std::ifstream infile(filename);
    return infile.good();
  }

}

#endif
