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
#include <immintrin.h>

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
     * @brief   reverse string
     * @note    assumes dest is pre-allocated
     */
    void reverse(const std::string &src, std::string &dest) 
    {
      assert(src.length() == dest.length());

      for ( int i = 0; i < src.length(); i++ )
        dest[src.length() - i - 1] = src.at(i);
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
    int32_t cigarScore (const std::string &cigar)
    {
      if (cigar.empty()) 
        return 0;

      int32_t score = 0;
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
            score += SCORE::mismatch * currentNumeric;
          else if (c == 'I')
            score += SCORE::ins * currentNumeric;
          else  // c == 'D'
            score += SCORE::del * currentNumeric;

          currentNumeric = 0;
        }
      }

      return score;
    }
  }

  template<typename T> struct simdUtils {};

  template<>
    struct simdUtils<int32_t> 
    {
      static void print_avx_num(const __m512i &var)
      {
        int32_t *val = (int32_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3], val[4],  val[5], 
            val[6],  val[7],  val[8],  val[9], val[10], val[11],
            val[12], val[13], val[14], val[15]);
      }

      static void print_avx_num(const __m256i &var)
      {
        int32_t *val = (int32_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3], val[4],  val[5], 
            val[6],  val[7]);
      }
    };

  template<>
    struct simdUtils<int16_t> 
    { 
      static void print_avx_num(const __m512i &var)
      {
        int16_t *val = (int16_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3],  val[4],  val[5], 
            val[6],  val[7],  val[8],  val[9],  val[10], val[11],
            val[12], val[13], val[14], val[15], val[16], val[17],
            val[18], val[19], val[20], val[21], val[22], val[23],
            val[24], val[25], val[26], val[27], val[28], val[29],
            val[30], val[31]);
      }

      static void print_avx_num(const __m256i &var)
      {
        int16_t *val = (int16_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3],  val[4],  val[5], 
            val[6],  val[7],  val[8],  val[9],  val[10], val[11],
            val[12], val[13], val[14], val[15]);
      }
    };

  template<>
    struct simdUtils<int8_t> 
    { 
      static void print_avx_num(const __m512i &var)
      {
        int8_t *val = (int8_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3],  val[4],  val[5], 
            val[6],  val[7],  val[8],  val[9],  val[10], val[11],
            val[12], val[13], val[14], val[15], val[16], val[17],
            val[18], val[19], val[20], val[21], val[22], val[23],
            val[24], val[25], val[26], val[27], val[28], val[29],
            val[30], val[31], val[32], val[33], val[34], val[35],
            val[36], val[37], val[38], val[39], val[40], val[41],
            val[42], val[43], val[44], val[45], val[46], val[47],
            val[48], val[49], val[50], val[51], val[52], val[53],
            val[54], val[55], val[56], val[57], val[58], val[59],
            val[60], val[61], val[62], val[63]);
      }

      static void print_avx_num(const __m256i &var)
      {
        int8_t *val = (int8_t*) &var;

        printf("Numerical: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i \n", 
            val[0],  val[1],  val[2],  val[3],  val[4],  val[5], 
            val[6],  val[7],  val[8],  val[9],  val[10], val[11],
            val[12], val[13], val[14], val[15], val[16], val[17],
            val[18], val[19], val[20], val[21], val[22], val[23],
            val[24], val[25], val[26], val[27], val[28], val[29],
            val[30], val[31]);
      }
    };

  /**
   * @brief     check if file is accessible
   */
  bool fileExists(const std::string &filename)
  {
    std::ifstream infile(filename);
    return infile.good();
  }

  /**
   * @brief       print few important execution env. variables
   */
  void showExecutionEnv()
  {
    std::cout << "--------" << "\n";

    {
      std::cout << "Assert() checks"; 
#ifdef NDEBUG
      std::cout << "\t\t\tOFF \n";
#else
      std::cout << "\t\t\tON \n";
#endif
    }

    {
      std::cout << "AVX SIMD support";
#ifdef PASGAL_ENABLE_AVX512
      std::cout << "\t\tON (AVX512) \n";
#elif PASGAL_ENABLE_AVX2
      std::cout << "\t\tON (AVX2) \n";
#else
      std::cout << "\t\tOFF\n";
#endif
    }

    {
      std::cout << "VTUNE profiling";
#ifdef VTUNE_SUPPORT
      std::cout << "\t\t\tON \n";
#else
      std::cout << "\t\t\tOFF\n";
#endif
    }

#pragma omp parallel
    {
      int tid = omp_get_thread_num();

      if (tid == 0) 
      {
        std::cout << "Default OMP thread count\t" << omp_get_num_threads() << "\n";
      }
    }

    std::cout << "--------\n" << std::endl;
  }

  /**
   * @brief                     print minimum, maximum and mean of thread runtimes
   * @param[in] threadRunTime   vector containing individual thread execution time
   * @return                    string containing min, max and mean
   */
  std::string printStats (std::vector<double> &threadRunTime)
  {
    auto min = *std::min_element(threadRunTime.begin(), threadRunTime.end());
    auto max = *std::max_element(threadRunTime.begin(), threadRunTime.end());
    auto mean = std::accumulate( threadRunTime.begin(), threadRunTime.end(), 0.0)/threadRunTime.size(); 

    std::string stats;

    stats += "min = " + std::to_string(min);
    stats += ", max = " + std::to_string(max);
    stats += ", mean = " + std::to_string(mean);
    stats += ", count = " + std::to_string(threadRunTime.size());

    return stats;
  }
}

#endif
