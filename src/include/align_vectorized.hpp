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
#include "aligned_allocator.hpp"

//define common SIMD operations
#define _ZERO       _mm512_setzero_epi32()
#define _ADD        _mm512_add_epi32
#define _MAX        _mm512_max_epi32 
#define _EQUAL      _mm512_cmpeq_epi32_mask
#define _SET1       _mm512_set1_epi32
#define _SET1_MASK  _mm512_mask_set1_epi32
#define _STORE      _mm512_store_epi32
#define _LOAD       _mm512_load_epi32
#define _BLEND      _mm512_mask_blend_epi32
#define SIMD_WIDTH  16 

namespace psgl
{
  /**
   * @brief   Supports phase 1 DP in forward direction
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    class Phase1_Vectorized
    {
      private:

        //input reads
        const std::vector<std::string> &readSet;

        //reference graph
        const CSR_char_container<VertexIdType, EdgeIdType> &graph;

        //for converting input reads into SOA to enable vectorization
        std::vector<char>  readSet_SOA;

        // pre-compute which graph vertices are connected with hop longer
        // than 'smallBufferWidth'
        std::vector<bool> withLongHop;

      public:

        //small temporary storage buffer for DP scores
        //should be a power of 2
        static constexpr size_t smallBufferWidth = 8; 

        //process these many vertical cells in a go
        //should be a power of 2
        static constexpr size_t matrixHeight = 8;   

        /**
         * @brief                   public constructor
         * @param[in]   readSet     vector of input query sequences to align
         * @param[in]   g           input reference graph
         */
        Phase1_Vectorized(const std::vector<std::string> &readSet, 
            const CSR_char_container<VertexIdType, EdgeIdType> &g) :
          readSet (readSet), graph (g)
        {
          //Type checks
          static_assert(std::is_same<ScoreType, int32_t>::value, "ScoreType needs to be int32_t for now");
          assert (sizeof(VertexIdType) <= 32);

          this->convertToSOA();
          this->computeLongHops();
        };

        /**
         * @brief                                 wrapper function for phase 1 routine 
         * @param[out]  outputBestScoreVector     vector to keep value and location of best scores,
         *                                        vector size is same as count of the reads
         * @note                                  reverse complement of reads are not handled
         *                                        inside this class
         */
        template <typename Vec>
          void alignToDAGLocal_Phase1_vectorized_wrapper(Vec &outputBestScoreVector) const
          {
            assert (outputBestScoreVector.size() == readSet.size());

            // modified containers to hold best-score info
            // vector elements aligned to 64 byte boundaries
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreVector (readSet.size() / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreColVector (readSet.size() / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreRowVector (readSet.size() / SIMD_WIDTH);

            this->alignToDAGLocal_Phase1_vectorized (_bestScoreVector, _bestScoreColVector, _bestScoreRowVector); 

            //parse best scores from vector registers
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeScores (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeCols   (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeRows   (SIMD_WIDTH);

            for (size_t i = 0; i < _bestScoreVector.size(); i++)
            {
              _STORE (storeScores.data(), _bestScoreVector[i]);
              _STORE (storeCols.data(), _bestScoreColVector[i]);
              _STORE (storeRows.data(), _bestScoreRowVector[i]);

              for (size_t j = 0; j < SIMD_WIDTH; j++)
              {
                outputBestScoreVector [i*SIMD_WIDTH + j].score = storeScores[j];
                outputBestScoreVector [i*SIMD_WIDTH + j].refColumnEnd = storeCols[j];
                outputBestScoreVector [i*SIMD_WIDTH + j].qryRowEnd = storeRows[j];

#ifdef DEBUG
                std::cout << "INFO, psgl::Phase1_Vectorized::alignToDAGLocal_Phase1_vectorized_wrapper, read # " << i * SIMD_WIDTH + j << ",  score = " << storeScores[j] << "\n";
#endif
              }
            }
          }

      private:

        /**
         * @brief       convert input reads characters into SOA to enable vectorization
         */
        void convertToSOA()
        {
          assert (readSet.size() > 0);
          assert (readSet_SOA.size() == 0);

          auto qLen = readSet[0].length();

          //assuming all read lengths are equal
          for (auto &read : readSet)
            assert (qLen == read.length());

          //assuming read lengths are multiple of vertical blocking length
          assert (qLen % matrixHeight == 0);

          //assuming count of reads is a multiple of SIMD_WIDTH
          assert (readSet.size() % SIMD_WIDTH == 0);

          //make space to save read charaacters
          readSet_SOA.resize ( readSet.size() * qLen);

          //re-arrange read characters for vectorized processing; SIMD_WIDTH reads at a time
          for (size_t readno = 0; readno < readSet.size(); readno += SIMD_WIDTH)
            for (size_t i = 0; i < qLen; i++)
              for (size_t j = 0; j < SIMD_WIDTH; j++)
                readSet_SOA [readno * qLen + i * SIMD_WIDTH + j] = readSet[readno + j][i];
        }

        /**
         * @brief     precompute which vertices are connected with long edge hops
         *            i.e. longer than 'smallBufferWidth'
         */
        void computeLongHops()
        {
          assert (withLongHop.size() == 0);

          this->withLongHop.resize(graph.numVertices, false);

          for(VertexIdType i = 0; i < graph.numVertices; i++)
          {
            for(auto j = graph.offsets_in[i]; j < graph.offsets_in[i+1]; j++)
            {
              auto from_pos = graph.adjcny_in[j];
              auto to_pos = i;

              //compare hop distance to 'smallBufferWidth'
              if (to_pos - from_pos >= this->smallBufferWidth)
                this->withLongHop[from_pos] = true;
            }
          }

#ifdef DEBUG
          auto trueCount = std::count(withLongHop.begin(), withLongHop.end(), true);
          std::cout << "INFO, psgl::Phase1_Vectorized::computeLongHops, fraction of hops that are long: " << trueCount * 1.0 / graph.numVertices << "\n";
#endif
        }

        /**
         * @brief                         execute first phase of alignment i.e. compute DP and 
         *                                find locations of the best alignment of each read
         * @param[out]  bestScores        best DP scores of reads
         * @param[out]  bestCols          columns where best alignment ends (for traceback later)
         * @param[out]  bestRows          rows where best alignment ends
         */
        template <typename Vec>
          void alignToDAGLocal_Phase1_vectorized (Vec &bestScores, Vec &bestCols, Vec &bestRows) const
          {
            size_t readCount = readSet.size();
            size_t readLength = readSet[0].length();

            //few checks
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
            __m512i match512    = _SET1 ((int32_t) SCORE::match);
            __m512i mismatch512 = _SET1 ((int32_t) SCORE::mismatch);
            __m512i del512      = _SET1 ((int32_t) SCORE::del);
            __m512i ins512      = _SET1 ((int32_t) SCORE::ins);

#pragma omp parallel
            {
              //initialize 2D vector called matrix of size 2 x graph size
              using AlignedVecType = std::vector <__m512i, aligned_allocator<__m512i, 64> >;

              //buffer to save selected columns (associated with long hops) of DP matrix
              std::vector< AlignedVecType > fartherColumns(graph.numVertices);

              //only allocate memory for selected vertices
              for(VertexIdType i = 0; i < graph.numVertices; i++)
              {
                if ( this->withLongHop[i] )
                  fartherColumns[i].resize(this->matrixHeight);
              }

              //buffer to save neighboring column scores
              std::vector< AlignedVecType > nearbyColumns (this->smallBufferWidth, AlignedVecType(this->matrixHeight));

              //buffer to save scores of last row in each iteration
              //one row for writing and one for reading
              std::vector< AlignedVecType > lastBatchRow (2, AlignedVecType (graph.numVertices));

              //buffer to save read charactes
              std::vector<int32_t, aligned_allocator<int32_t, 64> > readCharsInt (SIMD_WIDTH * this->matrixHeight);

              //process SIMD_WIDTH reads in a single iteration
#pragma omp for
              for (size_t readno = 0; readno < readCount; readno += SIMD_WIDTH)
              {
                size_t readBatch = readno/SIMD_WIDTH;

                __m512i bestScores512 = _ZERO;
                __m512i bestRows512   = _ZERO;
                __m512i bestCols512   = _ZERO;

                //reset DP 'lastBatchRow' buffer
                std::fill (lastBatchRow[1].begin(), lastBatchRow[1].end(), _ZERO);

                //iterate over read length (process more than 1 characters in batch)
                for (int32_t i = 0; i < readLength; i += this->matrixHeight)
                {
                  //loop counter 
                  size_t loopI = i / (this->matrixHeight);

                  //convert read character to int32_t
                  for (int32_t j = 0; j < SIMD_WIDTH * this->matrixHeight ; j++)
                  {
                    readCharsInt [j] = readSet_SOA [readno*readLength + i*SIMD_WIDTH + j];
                  }

                  //iterate over characters in reference graph
                  for (int32_t j = 0; j < graph.numVertices; j++)
                  {
                    //current reference character
                    __m512i graphChar = _SET1 ((int32_t) graph.vertex_label[j] );

                    //current best score, init to 0
                    __m512i currentMax512;

                    //iterate over read characters
                    for (size_t k = 0; k < this->matrixHeight; k++)
                    {
                      //load read characters
                      __m512i readChars = _LOAD ( &readCharsInt[k * SIMD_WIDTH] );

                      //current best score, init to 0
                      currentMax512 = _ZERO;

                      //see if query and reference character match
                      __mmask16 compareChar = _EQUAL (readChars, graphChar);
                      __m512i sub512 = _BLEND (compareChar, mismatch512, match512);

                      //match-mismatch edit
                      currentMax512 = _MAX (currentMax512, sub512); //local alignment can also start with a match at this char 

                      //iterate over graph neighbors
                      if (k == 0)
                      {
                        for(size_t l = graph.offsets_in[j]; l < graph.offsets_in[j+1]; l++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit = _ADD ( lastBatchRow[(loopI - 1) & 1][ graph.adjcny_in[l] ], sub512);
                          currentMax512 = _MAX (currentMax512, substEdit); 

                          //paths with deletion edit
                          __m512i delEdit;

                          if (j - graph.adjcny_in[l] < smallBufferWidth)
                            delEdit = _ADD ( nearbyColumns[graph.adjcny_in[l] & (smallBufferWidth-1) ][k], del512);
                          else
                            delEdit = _ADD ( fartherColumns[graph.adjcny_in[l]][k], del512);

                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (lastBatchRow[(loopI - 1) & 1][j], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }
                      else
                      {
                        for(size_t l = graph.offsets_in[j]; l < graph.offsets_in[j+1]; l++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit;

                          //paths with deletion edit
                          __m512i delEdit;

                          if (j - graph.adjcny_in[l] < smallBufferWidth)
                          {
                            substEdit = _ADD ( nearbyColumns[graph.adjcny_in[l] & (smallBufferWidth-1)][k-1], sub512);
                            delEdit = _ADD ( nearbyColumns[graph.adjcny_in[l] & (smallBufferWidth-1)][k], del512);
                          }
                          else
                          {
                            substEdit = _ADD ( fartherColumns[graph.adjcny_in[l]][k-1], sub512);
                            delEdit = _ADD ( fartherColumns[graph.adjcny_in[l]][k], del512);
                          }

                          currentMax512 = _MAX (currentMax512, substEdit); 
                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (nearbyColumns[j & (smallBufferWidth-1)][k-1], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }

                      //update best score observed yet
                      bestScores512 = _MAX (currentMax512, bestScores512);

                      //on which lanes is the best score updated
                      __mmask16 updated = _EQUAL (currentMax512, bestScores512);

                      //update row and column values accordingly
                      bestRows512 = _SET1_MASK (bestRows512, updated, (int32_t) (i + k));
                      bestCols512 = _SET1_MASK (bestCols512, updated, (int32_t) j);

                      //save current score in small buffer
                      nearbyColumns[j & (smallBufferWidth-1)][k] = currentMax512;

                      //save current score in large buffer if connected thru long hop
                      if ( this->withLongHop[j] )
                        fartherColumns[j][k] = currentMax512;
                    }

                    //save last score for next row-wise iteration
                    lastBatchRow[loopI & 1][j] = currentMax512; 

                  } // end of row computation
                } // end of DP

                bestScores[readBatch] = bestScores512;
                bestRows[readBatch]   = bestRows512; 
                bestCols[readBatch]   = bestCols512; 

              } // all reads done
            } //end of omp parallel

            auto tick2 = __rdtsc();

            std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_vectorized, CPU cycles spent in phase 1 (without wrapper) = " << tick2 - tick1
              << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";

#ifdef VTUNE_SUPPORT
            __itt_pause();
#endif
          }

    };

  /**
   * @brief   Supports phase 1 DP in reverse direction
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    class Phase1_Rev_Vectorized
    {
      private:

        //input reads
        const std::vector<std::string> &readSet;

        //reference graph
        const CSR_char_container<VertexIdType, EdgeIdType> &graph;

        //small temporary storage buffer for DP scores
        //should be a power of 2
        static constexpr size_t smallBufferWidth = Phase1_Vectorized<ScoreType, VertexIdType, EdgeIdType>::smallBufferWidth; 

        //process these many vertical cells in a go
        //should be a power of 2
        static constexpr size_t matrixHeight = Phase1_Vectorized<ScoreType, VertexIdType, EdgeIdType>::matrixHeight;   

        //for converting input reads into SOA to enable vectorization
        std::vector<char>  readSet_SOA;

        // pre-compute which graph vertices are connected with hop longer
        // than 'smallBufferWidth'
        std::vector<bool> withLongHop;

      public:

        /**
         * @brief                   public constructor
         * @param[in]   readSet     vector of input query sequences to align
         * @param[in]   g           input reference graph
         */
        Phase1_Rev_Vectorized(const std::vector<std::string> &readSet, 
            const CSR_char_container<VertexIdType, EdgeIdType> &g) :
          readSet (readSet), graph (g)
        {
          //Type checks
          static_assert(std::is_same<ScoreType, int32_t>::value, "ScoreType needs to be int32_t for now");
          assert (sizeof(VertexIdType) <= 32);

          this->convertToSOA();
          this->computeLongHops();
        };

        /**
         * @brief                                 wrapper function for phase 1 reverse DP routine 
         * @param[out]  outputBestScoreVector     vector to keep value and begin location of best scores,
         *                                        vector size is same as count of the reads
         * @note                                  reverse complement of reads are not handled
         *                                        inside this class
         */
        template <typename Vec>
          void alignToDAGLocal_Phase1_rev_vectorized_wrapper(Vec &outputBestScoreVector) const
          {
            assert (outputBestScoreVector.size() == readSet.size());

            // modified containers to hold best-score info
            // vector elements aligned to 64 byte boundaries
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreVector (readSet.size() / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreColVector (readSet.size() / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreRowVector (readSet.size() / SIMD_WIDTH);

            this->alignToDAGLocal_Phase1_rev_vectorized (outputBestScoreVector, _bestScoreVector, _bestScoreColVector, _bestScoreRowVector); 

            //parse best scores from vector registers
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeScores (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeCols   (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeRows   (SIMD_WIDTH);

            for (size_t i = 0; i < _bestScoreVector.size(); i++)
            {
              _STORE (storeScores.data(), _bestScoreVector[i]);
              _STORE (storeCols.data(), _bestScoreColVector[i]);
              _STORE (storeRows.data(), _bestScoreRowVector[i]);

              for (size_t j = 0; j < SIMD_WIDTH; j++)
              {
                assert (outputBestScoreVector [i*SIMD_WIDTH + j].score == storeScores[j] - 1);  //offset by 1
                outputBestScoreVector [i*SIMD_WIDTH + j].refColumnStart = storeCols[j];
                outputBestScoreVector [i*SIMD_WIDTH + j].qryRowStart = readSet[i*SIMD_WIDTH + j].length() - 1 - storeRows[j];   //convert rev. direction to fwd.

#ifdef DEBUG
                std::cout << "INFO, psgl::Phase1_Rev_Vectorized::alignToDAGLocal_Phase1_rev_vectorized_wrapper, read # " << i * SIMD_WIDTH + j << ",  score = " << storeScores[j] << "\n";
#endif
              }
            }
          }

      private:

        /**
         * @brief       convert input reads characters into SOA to enable vectorization
         */
        void convertToSOA()
        {
          assert (readSet.size() > 0);
          assert (readSet_SOA.size() == 0);

          auto qLen = readSet[0].length();

          //assuming all read lengths are equal
          for (auto &read : readSet)
            assert (qLen == read.length());

          //assuming read lengths are multiple of vertical blocking length
          assert (qLen % matrixHeight == 0);

          //assuming count of reads is a multiple of SIMD_WIDTH
          assert (readSet.size() % SIMD_WIDTH == 0);

          //make space to save read charaacters
          readSet_SOA.resize ( readSet.size() * qLen);

          //re-arrange read characters for vectorized processing; SIMD_WIDTH reads at a time
          for (size_t readno = 0; readno < readSet.size(); readno += SIMD_WIDTH)
            for (size_t i = 0; i < qLen; i++)
              for (size_t j = 0; j < SIMD_WIDTH; j++)
                readSet_SOA [readno * qLen + i * SIMD_WIDTH + j] = readSet[readno + j][i];
        }

        /**
         * @brief     precompute which vertices are connected with long edge hops
         *            i.e. longer than 'smallBufferWidth'
         */
        void computeLongHops()
        {
          assert (withLongHop.size() == 0);

          this->withLongHop.resize(graph.numVertices, false);

          for(VertexIdType i = 0; i < graph.numVertices; i++)
          {
            for(auto j = graph.offsets_in[i]; j < graph.offsets_in[i+1]; j++)
            {
              auto from_pos = graph.adjcny_in[j];
              auto to_pos = i;

              //compare hop distance to 'smallBufferWidth'
              if (to_pos - from_pos >= this->smallBufferWidth)
                this->withLongHop[to_pos] = true;
            }
          }

#ifdef DEBUG
          auto trueCount = std::count(withLongHop.begin(), withLongHop.end(), true);
          std::cout << "INFO, psgl::Phase1_Rev_Vectorized::computeLongHops, fraction of hops that are long: " << trueCount * 1.0 / graph.numVertices << "\n";
#endif
        }

        /**
         * @brief                               execute reverse variant of first phase of alignment i.e. 
         *                                      compute reverse DP and find begin locations of the best 
         *                                      alignment of each read
         * @param[in]   outputBestScoreVector   best scores and end locations computed during forward DP
         * @param[out]  bestScores              best DP scores of reads
         * @param[out]  bestCols                columns where best alignment starts
         * @param[out]  bestRows                rows where best alignment starts
         */
        template <typename Vec1, typename Vec2>
          void alignToDAGLocal_Phase1_rev_vectorized (const Vec1 &outputBestScoreVector,
                                                      Vec2 &bestScores, Vec2 &bestCols, Vec2 &bestRows) const
          {
            size_t readCount = readSet.size();
            size_t readLength = readSet[0].length();

            //few checks
            assert (bestScores.size() == readCount / SIMD_WIDTH);
            assert (bestCols.size() == readCount / SIMD_WIDTH);
            assert (bestRows.size() == readCount / SIMD_WIDTH);
            assert (readSet_SOA.size() == readCount * readLength);
            assert (readCount % SIMD_WIDTH == 0);

//#ifdef VTUNE_SUPPORT
            //__itt_resume();
//#endif

            //for time profiling within phase 1
            auto tick1 = __rdtsc();

            //init best score vector to zero bits
            std::fill (bestScores.begin(), bestScores.end(), _ZERO);
            std::fill (bestCols.begin(), bestCols.end(), _ZERO);
            std::fill (bestRows.begin(), bestRows.end(), _ZERO);

            //init score simd vectors
            __m512i match512    = _SET1 ((int32_t) SCORE::match);
            __m512i mismatch512 = _SET1 ((int32_t) SCORE::mismatch);
            __m512i del512      = _SET1 ((int32_t) SCORE::del);
            __m512i ins512      = _SET1 ((int32_t) SCORE::ins);

#pragma omp parallel
            {
              //initialize 2D vector called matrix of size 2 x graph size
              using AlignedVecType = std::vector <__m512i, aligned_allocator<__m512i, 64> >;

              //buffer to save selected columns (associated with long hops) of DP matrix
              std::vector< AlignedVecType > fartherColumns(graph.numVertices);

              //only allocate memory for selected vertices
              for(VertexIdType i = 0; i < graph.numVertices; i++)
              {
                if ( this->withLongHop[i] )
                  fartherColumns[i].resize(this->matrixHeight);
              }

              //buffer to save neighboring column scores
              std::vector< AlignedVecType > nearbyColumns (this->smallBufferWidth, AlignedVecType(this->matrixHeight));

              //buffer to save scores of last row in each iteration
              //one row for writing and one for reading
              std::vector< AlignedVecType > lastBatchRow (2, AlignedVecType (graph.numVertices));

              //buffer to save read charactes
              std::vector<int32_t, aligned_allocator<int32_t, 64> > readCharsInt (SIMD_WIDTH * this->matrixHeight);

              //process SIMD_WIDTH reads in a single iteration
#pragma omp for
              for (size_t readno = 0; readno < readCount; readno += SIMD_WIDTH)
              {
                size_t readBatch = readno/SIMD_WIDTH;

                __m512i fwdBestCols512;
                __m512i fwdBestRows512;

                //parse alignment locations of forward DP
                {
                  std::vector<int32_t, aligned_allocator<int32_t, 64> > fwdBestCols   (SIMD_WIDTH);
                  std::vector<int32_t, aligned_allocator<int32_t, 64> > fwdBestRows   (SIMD_WIDTH);

                  for (size_t j = 0; j < SIMD_WIDTH; j++)
                  {
                    fwdBestCols[j] = outputBestScoreVector[readno + j].refColumnEnd;
                    fwdBestRows[j] = readLength - 1 - outputBestScoreVector[readno + j].qryRowEnd;
                  }

                  fwdBestCols512 = _LOAD ( fwdBestCols.data() );
                  fwdBestRows512 = _LOAD ( fwdBestRows.data() );
                }

                __m512i bestScores512 = _ZERO;
                __m512i bestRows512   = _ZERO;
                __m512i bestCols512   = _ZERO;

                //reset DP 'lastBatchRow' buffer
                std::fill (lastBatchRow[1].begin(), lastBatchRow[1].end(), _ZERO);

                //iterate over read length (process more than 1 characters in batch)
                for (int32_t i = 0; i < readLength; i += this->matrixHeight)
                {
                  //loop counter 
                  size_t loopI = i / (this->matrixHeight);

                  //convert read character to int32_t
                  for (int32_t j = 0; j < SIMD_WIDTH * this->matrixHeight ; j++)
                  {
                    readCharsInt [j] = readSet_SOA [readno*readLength + i*SIMD_WIDTH + j];
                  }

                  //iterate over characters in reference graph
                  for (int32_t j = graph.numVertices - 1; j >= 0; j--)
                  {
                    //current reference character
                    __m512i graphChar = _SET1 ((int32_t) graph.vertex_label[j] );

                    //current best score
                    __m512i currentMax512;

                    //iterate over read characters
                    for (size_t k = 0; k < this->matrixHeight; k++)
                    {
                      //load read characters
                      __m512i readChars = _LOAD ( &readCharsInt[k * SIMD_WIDTH] );

                      //current best score, init to 0
                      currentMax512 = _ZERO;

                      //see if query and reference character match
                      __mmask16 compareChar = _EQUAL (readChars, graphChar);
                      __m512i sub512 = _BLEND (compareChar, mismatch512, match512);

                      //match-mismatch edit
                      currentMax512 = _MAX (currentMax512, sub512); //local alignment can also start with a match at this char 

                      //iterate over graph neighbors
                      if (k == 0)
                      {
                        for(size_t l = graph.offsets_out[j]; l < graph.offsets_out[j+1]; l++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit = _ADD ( lastBatchRow[(loopI - 1) & 1][ graph.adjcny_out[l] ], sub512);
                          currentMax512 = _MAX (currentMax512, substEdit); 

                          //paths with deletion edit
                          __m512i delEdit;

                          if (graph.adjcny_out[l] - j < smallBufferWidth)
                            delEdit = _ADD ( nearbyColumns[graph.adjcny_out[l] & (smallBufferWidth-1) ][k], del512);
                          else
                            delEdit = _ADD ( fartherColumns[graph.adjcny_out[l]][k], del512);

                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (lastBatchRow[(loopI - 1) & 1][j], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }
                      else
                      {
                        for(size_t l = graph.offsets_out[j]; l < graph.offsets_out[j+1]; l++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit;

                          //paths with deletion edit
                          __m512i delEdit;

                          if (graph.adjcny_out[l] - j < smallBufferWidth)
                          {
                            substEdit = _ADD ( nearbyColumns[graph.adjcny_out[l] & (smallBufferWidth-1)][k-1], sub512);
                            delEdit = _ADD ( nearbyColumns[graph.adjcny_out[l] & (smallBufferWidth-1)][k], del512);
                          }
                          else
                          {
                            substEdit = _ADD ( fartherColumns[graph.adjcny_out[l]][k-1], sub512);
                            delEdit = _ADD ( fartherColumns[graph.adjcny_out[l]][k], del512);
                          }

                          currentMax512 = _MAX (currentMax512, substEdit); 
                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (nearbyColumns[j & (smallBufferWidth-1)][k-1], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }

                      //update best score observed yet
                      bestScores512 = _MAX (currentMax512, bestScores512);

                      //on which lanes is the best score updated
                      __mmask16 updated = _EQUAL (currentMax512, bestScores512);

                      //update row and column values accordingly
                      bestCols512 = _SET1_MASK (bestCols512, updated, (int32_t) j);
                      bestRows512 = _SET1_MASK (bestRows512, updated, (int32_t) (i + k));

                      //detect and update score of the optimal alignment we want
                      {
                        __m512i currentCol = _SET1 ( (int32_t) j ); 
                        __m512i currentRow = _SET1 ( (int32_t) (i + k) );
                        
                        __mmask16 compareCell = _EQUAL (fwdBestCols512, currentCol) & _EQUAL (fwdBestRows512, currentRow);

                        currentMax512 = _SET1_MASK (currentMax512, compareCell, (int32_t) (SCORE::match + 1)); 
                      }

                      //save current score in small buffer
                      nearbyColumns[j & (smallBufferWidth-1)][k] = currentMax512;

                      //save current score in large buffer if connected thru long hop
                      if ( this->withLongHop[j] )
                        fartherColumns[j][k] = currentMax512;
                    }

                    //save last score for next row-wise iteration
                    lastBatchRow[loopI & 1][j] = currentMax512; 

                  } // end of row computation
                } // end of DP

                bestScores[readBatch] = bestScores512;
                bestRows[readBatch]   = bestRows512; 
                bestCols[readBatch]   = bestCols512; 

              } // all reads done
            } //end of omp parallel

            auto tick2 = __rdtsc();

            std::cout << "TIMER, psgl::Phase1_Rev_Vectorized::alignToDAGLocal_Phase1_rev_vectorized , CPU cycles spent in phase 1-R (without wrapper) = " << tick2 - tick1
              << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";

//#ifdef VTUNE_SUPPORT
            //__itt_pause();
//#endif
          }

    };
}

#undef _ZERO     
#undef _ADD      
#undef _MAX      
#undef _EQUAL    
#undef _SET1      
#undef _SET1_MASK
#undef _STORE    
#undef _LOAD     
#undef _BLEND
#undef SIMD_WIDTH

#endif
