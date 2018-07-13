/**
 * @file    align_vectorized.hpp
 * @brief   vectorized routines to perform alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ALIGN_VEC_HPP
#define GRAPH_ALIGN_VEC_HPP

#include <immintrin.h>

#include "graphLoad.hpp"
#include "csr_char.hpp"
#include "graph_iter.hpp"
#include "base_types.hpp"
#include "utils.hpp"

//External includes
#include "kseq.h"
#include "prettyprint.hpp"

//KSEQ_INIT(gzFile, gzread)

#define psgl_max(a,b) (((a)>(b))?(a):(b))

//define common SIMD operations
#define _ZERO       _mm512_setzero_epi32()
#define _ADD        _mm512_add_epi32
#define _MAX        _mm512_max_epi32 
#define _EQUAL      _mm512_cmpeq_epi32_mask
#define _SET        _mm512_set1_epi32
#define _SET_MASK   _mm512_mask_set1_epi32
#define _STORE      _mm512_store_epi32
#define _LOAD       _mm512_load_epi32
#define SIMD_WIDTH  16 

namespace psgl
{
  /**
   * @brief                         execute first phase of alignment i.e. compute DP and 
   *                                find locations of the best alignment of each read
   * @tparam[in]  ScoreType         a signed numeric type to store scores in DP matrix
   * @param[in]   readSet_SOA       read characters arranged as S.O.A. for vectorization
   * @param[in]   graph
   * @param[out]  bestScores        best DP scores of reads
   * @param[out]  bestCols          columns where best alignment ends (for traceback later)
   * @param[out]  bestRows          rows where best alignment ends
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal_Phase1_vectorized( size_t readCount, size_t readLength,
                                            const std::vector<char> &readSet_SOA,
                                            const CSR_char_container<VertexIdType, EdgeIdType> &graph,
                                            std::vector< __m512i > &bestScores,
                                            std::vector< __m512i > &bestCols,
                                            std::vector< __m512i > &bestRows)
    {
      //few checks
      static_assert(std::is_same<ScoreType, int32_t>::value, "ScoreType needs to be int32_t for now");
      assert (bestScores.size() == readCount / SIMD_WIDTH);
      assert (bestCols.size() == readCount / SIMD_WIDTH);
      assert (bestRows.size() == readCount / SIMD_WIDTH);
      assert (readSet_SOA.size() == readCount * readLength);
      assert (readCount % SIMD_WIDTH == 0);

#ifdef VTUNE_SUPPORT
      __itt_resume();
#endif

      //for time profiling within phase 1
      auto tick1 = __rdtsc();

      //init best score vector to zero bits
      std::fill (bestScores.begin(), bestScores.end(), _ZERO);
      std::fill (bestCols.begin(), bestCols.end(), _ZERO);
      std::fill (bestRows.begin(), bestRows.end(), _ZERO);

      //init score simd vectors
      __m512i match512    = _SET (SCORE::match);
      __m512i mismatch512 = _SET (SCORE::mismatch);
      __m512i del512      = _SET (SCORE::del);
      __m512i ins512      = _SET (SCORE::ins);

#pragma omp parallel
      {
        //initialize 2D vector called matrix of size 2 x graph size
        std::vector< std::vector< __m512i > > matrix(2, std::vector< __m512i >(graph.numVertices));

        //buffer to save read charactes
        std::vector<int32_t> readCharsInt (SIMD_WIDTH);

        //process SIMD_WIDTH reads in a single iteration
#pragma omp for
        for (size_t readno = 0; readno < readCount; readno += SIMD_WIDTH)
        {
          size_t readBatch = readno/SIMD_WIDTH;

          //reset buffer
          std::fill(matrix[1].begin(), matrix[1].end(), _ZERO);

          //iterate over read length
          for (std::size_t i = 0; i < readLength; i++)
          {
            //convert read character to int32_t
            for (int j = 0; j < SIMD_WIDTH; j++)
              readCharsInt [j] = readSet_SOA [readno*readLength + i*SIMD_WIDTH + j];

            //load read characters
            __m512i readChars = _LOAD (readCharsInt.data());

            //iterate over characters in reference graph
            for (std::size_t j = 0; j < graph.numVertices; j++)
            {
              //current reference character
              __m512i graphChar = _SET ((int32_t) graph.vertex_label[j] );

              //current best score, init to 0
              __m512i currentMax = _ZERO;

              //see if query and reference character match
              __mmask16 compareChar = _EQUAL (readChars, graphChar);
              __m512i sub512 = _mm512_mask_blend_epi32 (compareChar, mismatch512, match512);

              //match-mismatch edit
              currentMax = _MAX (currentMax, sub512); //local alignment can also start with a match at this char 

              for(auto k = graph.offsets_in[j]; k < graph.offsets_in[j+1]; k++)
              {
                //paths with match mismatch edit
                __m512i substEdit = _ADD ( matrix[(i-1) & 1][graph.adjcny_in[k]], sub512);
                currentMax = _MAX (currentMax, substEdit); 
                //'& 1' is same as doing modulo 2

                //paths with deletion edit
                __m512i delEdit = _ADD ( matrix[i & 1][graph.adjcny_in[k]], del512);
                currentMax = _MAX (currentMax, delEdit); 
              }

              //insertion edit
              __m512i insEdit = _ADD (matrix[(i-1) & 1][j], ins512);
              currentMax = _MAX (currentMax, insEdit);

              //update best score observed yet
              bestScores[readBatch] = _MAX (currentMax, bestScores[readBatch]);

              //on which lanes is the best score updated
              __mmask16 updated = _EQUAL (currentMax, bestScores[readBatch]);

              //update row and column values accordingly
              bestRows[readBatch] = _SET_MASK (bestRows[readBatch], updated, (int32_t) i); 
              bestCols[readBatch] = _SET_MASK (bestCols[readBatch], updated, (int32_t) j); 

              matrix[i & 1][j] = currentMax;
            } // end of row computation
          } // end of DP
        } // all reads done
      } //end of omp parallel

      auto tick2 = __rdtsc();

      std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_vectorized, CPU cycles spent in phase 1 (without wrapper) = " << tick2 - tick1
                << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";

#ifdef VTUNE_SUPPORT
        __itt_pause();
#endif
    }

  /**
   * @brief                                 wrapper function for phase 1 routine 
   * @tparam[in]  ScoreType                 a signed numeric type to store scores in DP matrix
   * @param[in]   readSet                   vector of input query sequences to align
   * @param[in]   graph
   * @param[out]  outputBestScoreVector     vector to keep value and location of best scores,
   *                                        vector size is same as count of the reads
   * @note                                  reverse complement of the read is not handled here
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal_Phase1_vectorized_wrapper( const std::vector<std::string> &readSet,
        const CSR_char_container<VertexIdType, EdgeIdType> &graph,
        std::vector< BestScoreInfo<ScoreType, VertexIdType> > &outputBestScoreVector)
    {
      /**
       * Assumptions for now:
       *  - count of reads is a multiple of SIMD_WIDTH
       *  - all read lengths are equal
       */
      assert (readSet.size() > 0);
      assert (readSet.size() % SIMD_WIDTH == 0);
      assert (outputBestScoreVector.size() == readSet.size());
      auto qLen = readSet[0].length();
      for (auto &read : readSet)
        assert (qLen == read.length());

      //Convert input read vector into SOA for vectorization
      std::vector<char>  readSet_SOA ( readSet.size() * qLen);

      //re-arrange read characters for vectorized processing; SIMD_WIDTH reads at a time
      for (size_t readno = 0; readno < readSet.size(); readno += SIMD_WIDTH)
        for (size_t i = 0; i < qLen; i++)
          for (size_t j = 0; j < SIMD_WIDTH; j++)
            readSet_SOA [readno * qLen + i * SIMD_WIDTH + j] = readSet[readno + j][i];

      //Assumption
      assert (sizeof(VertexIdType) <= 32);

      // modified containers to hold best-score info
      std::vector <__m512i> _bestScoreVector (readSet.size() / SIMD_WIDTH);
      std::vector <__m512i> _bestScoreColVector (readSet.size() / SIMD_WIDTH);
      std::vector <__m512i> _bestScoreRowVector (readSet.size() / SIMD_WIDTH);

      alignToDAGLocal_Phase1_vectorized<ScoreType> (readSet.size(), qLen,
                                                    readSet_SOA, graph, 
                                                    _bestScoreVector, _bestScoreColVector, _bestScoreRowVector); 

      //parse best scores from vector registers
      std::vector<int32_t> storeScores (SIMD_WIDTH);
      std::vector<int32_t> storeCols   (SIMD_WIDTH);
      std::vector<int32_t> storeRows   (SIMD_WIDTH);

      for (size_t i = 0; i < _bestScoreVector.size(); i++)
      {
        _STORE (storeScores.data(), _bestScoreVector[i]);
        _STORE (storeCols.data(), _bestScoreColVector[i]);
        _STORE (storeRows.data(), _bestScoreRowVector[i]);

        for (size_t j = 0; j < SIMD_WIDTH; j++)
        {
          outputBestScoreVector [i*SIMD_WIDTH + j].score = storeScores[j];
          outputBestScoreVector [i*SIMD_WIDTH + j].refColumn = storeCols[j];
          outputBestScoreVector [i*SIMD_WIDTH + j].qryRow = storeRows[j];

          //std::cout << "INFO, psgl::alignToDAGLocal_Phase1_vectorized_wrapper, aligning read #" << i*SIMD_WIDTH + j << ", score " << storeScores[j] << "\n";
        }
      }
    }

}

#undef _ZERO     
#undef _ADD      
#undef _MAX      
#undef _EQUAL    
#undef _SET      
#undef _SET_MASK 
#undef _STORE    
#undef _LOAD     
#undef SIMD_WIDTH

#endif
