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
#define DUMMY       'B'

namespace psgl
{
  /**
   * @brief   Supports phase 1 DP in forward direction
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    class Phase1_Vectorized
    {
      private:

        //reference graph
        const CSR_char_container<VertexIdType, EdgeIdType> &graph;

        // pre-compute which graph vertices are connected with hop longer
        // than 'blockWidth'
        std::vector<bool> withLongHop;

        //for converting input reads into SOA to enable vectorization
        std::vector<char> readSetSOA;

        //cumulative read batch sizes
        std::vector<size_t>  readSetSOAPrefixSum;

        //input reads
        const std::vector<std::string> &readSet;

        //read lengths in their 'sorted order'
        std::vector<size_t> sortedReadLengths;

        //sorted permutation order of input reads
        std::vector<size_t> sortedReadOrder;

      public:

        //small temporary storage buffer for DP scores
        //should be a power of 2
        static constexpr size_t blockWidth = 8; 

        //process these many vertical cells in a go
        //should be a power of 2
        static constexpr size_t blockHeight = 16;   

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

          this->sortReadsForLoadBalance();
          this->convertToSOA();
          this->computeLongHops();
        };

        /**
         * @brief                                 wrapper function for phase 1 routine 
         * @param[out]  outputBestScoreVector     vector to keep value and location of best scores,
         *                                        vector size is same as count of the reads
         * @note                                  this class won't care about rev. complement of sequences
         */
        template <typename Vec>
          void alignToDAGLocal_Phase1_vectorized_wrapper(Vec &outputBestScoreVector) const
          {
            assert (outputBestScoreVector.size() == readSet.size());

            // modified containers to hold best-score info
            // vector elements aligned to 64 byte boundaries
            std::size_t countReadBatches = std::ceil (readSet.size() * 1.0 / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreVector (countReadBatches);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreColVector (countReadBatches);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreRowVector (countReadBatches);

            // execute the alignment routine
            this->alignToDAGLocal_Phase1_vectorized (_bestScoreVector, _bestScoreColVector, _bestScoreRowVector); 

            //parse best scores from vector registers
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeScores (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeCols   (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeRows   (SIMD_WIDTH);

            for (size_t i = 0; i < countReadBatches; i++)
            {
              _STORE (storeScores.data(), _bestScoreVector[i]);
              _STORE (storeCols.data(), _bestScoreColVector[i]);
              _STORE (storeRows.data(), _bestScoreRowVector[i]);

              for (size_t j = 0; j < SIMD_WIDTH; j++)
              {
                if (i*SIMD_WIDTH + j < readSet.size())
                {
                  auto originalReadId = sortedReadOrder[i*SIMD_WIDTH + j];

                  outputBestScoreVector[originalReadId].score         = storeScores[j];
                  outputBestScoreVector[originalReadId].refColumnEnd  = storeCols[j];
                  outputBestScoreVector[originalReadId].qryRowEnd     = storeRows[j];

#ifdef DEBUG
                  std::cout << "INFO, psgl::Phase1_Vectorized::alignToDAGLocal_Phase1_vectorized_wrapper, read # " << originalReadId << ",  score = " << storeScores[j] << "\n";
#endif
                }
              }
            }
          }

      private:

        /**
         * @brief       compute the sorted order of sequences for load balancing
         * @details     sorting is done in decreasing length order
         */
        void sortReadsForLoadBalance()
        {
          typedef std::pair<size_t, size_t> pair_t;

          //vector of tuples of length and index of reads
          std::vector<pair_t> lengthTuples;
           
          for(size_t i = 0; i < this->readSet.size(); i++)
            lengthTuples.emplace_back (readSet[i].length(), i);

          //sort in descending order, longer reads first
          std::sort (lengthTuples.begin(), lengthTuples.end(), std::greater<pair_t>() );

          for(auto &e : lengthTuples)
          {
            this->sortedReadLengths.push_back (e.first);
            this->sortedReadOrder.push_back (e.second);
          }
        }

        /**
         * @brief       convert input reads characters into SOA to enable vectorization
         */
        void convertToSOA()
        {
          assert (readSet.size() > 0);
          assert (sortedReadOrder.size() == readSet.size());
          assert (readSetSOA.size() == 0);

          /**
           * Requirements from the padding process
           * - each read length be a multiple of 'blockHeight'
           * - count of reads be a multiple of SIMD_WIDTH
           * - each batch of SIMD_WIDTH reads should be of equal length
           */

          //re-arrange read characters for vectorized processing; 
          //SIMD_WIDTH reads at a time, sorted in decreasing order by length
          
          auto readCount = readSet.size();

          readSetSOAPrefixSum.push_back(0);

          for (size_t i = 0; i < readCount; i += SIMD_WIDTH)
          {
            auto batchLength = readSet[sortedReadOrder[i]].length();  //longest read in this batch
            batchLength += blockHeight - 1 - (batchLength - 1) % blockHeight; //round-up

            for (size_t j = 0; j < batchLength; j++)
              for (size_t k = 0; k < SIMD_WIDTH; k++)
                if ( i + k < readCount && j < readSet[sortedReadOrder[i+k]].length() )
                  readSetSOA.push_back ( readSet[sortedReadOrder[i + k]][j] );
                else
                  readSetSOA.push_back (DUMMY);

            readSetSOAPrefixSum.push_back (readSetSOA.size());
          }

          std::size_t countReadBatches = std::ceil (readSet.size() * 1.0 / SIMD_WIDTH);
          assert (readSetSOAPrefixSum.size() == countReadBatches + 1);
        }

        /**
         * @brief     precompute which vertices are connected with long edge hops
         *            i.e. longer than 'blockWidth'
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

              //compare hop distance to 'blockWidth'
              if (to_pos - from_pos >= this->blockWidth)
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
            std::size_t readCount = readSet.size();
            std::size_t countReadBatches = std::ceil (readCount * 1.0 / SIMD_WIDTH);

            //few checks
            assert (bestScores.size() == countReadBatches);
            assert (bestCols.size() == countReadBatches);
            assert (bestRows.size() == countReadBatches);

#ifdef VTUNE_SUPPORT
            __itt_resume();
#endif

            //init best score vector to zero bits
            std::fill (bestScores.begin(), bestScores.end(), _ZERO);
            std::fill (bestCols.begin(), bestCols.end(), _ZERO);
            std::fill (bestRows.begin(), bestRows.end(), _ZERO);

            //init score simd vectors
            __m512i match512    = _SET1 ((int32_t) SCORE::match);
            __m512i mismatch512 = _SET1 ((int32_t) SCORE::mismatch);
            __m512i del512      = _SET1 ((int32_t) SCORE::del);
            __m512i ins512      = _SET1 ((int32_t) SCORE::ins);

            std::vector<double> threadTimings (omp_get_max_threads(), 0);

#pragma omp parallel
            {
#pragma omp barrier
              threadTimings[omp_get_thread_num()] = omp_get_wtime();

              //create local copy of graph for faster access
              const CSR_char_container<VertexIdType, EdgeIdType> graphLocal = this->graph;
              const std::vector<bool> withLongHopLocal = withLongHop;

              //type def. for memory-aligned vector allocation for SIMD instructions
              using AlignedVecType = std::vector <__m512i, aligned_allocator<__m512i, 64> >;

              //buffer to save selected columns (associated with long hops) of DP matrix
              std::size_t countLongHops = std::count (withLongHopLocal.begin(), withLongHopLocal.end(), true);
              AlignedVecType fartherColumnsBuffer (countLongHops * this->blockHeight);

              //for convenient access to 2D buffer
              std::vector<__m512i*> fartherColumns(graphLocal.numVertices);
              {
                size_t j = 0;

                for(VertexIdType i = 0; i < graphLocal.numVertices; i++)
                  if ( withLongHopLocal[i] )
                    fartherColumns[i] = &fartherColumnsBuffer[ (j++) * this->blockHeight ];
              }

              //buffer to save neighboring column scores
              AlignedVecType nearbyColumnsBuffer (this->blockWidth * this->blockHeight);

              //for convenient access to 2D buffer
              std::vector<__m512i*> nearbyColumns (this->blockWidth);
              {
                for (std::size_t i = 0; i < this->blockWidth; i++)
                  nearbyColumns[i] = &nearbyColumnsBuffer[i * this->blockHeight];
              }

              //buffer to save scores of last row in each iteration
              //one row for writing and one for reading
              AlignedVecType lastBatchRowBuffer (2 * graphLocal.numVertices);

              //for convenient access to 2D buffer
              std::vector<__m512i*> lastBatchRow (2);
              {
                lastBatchRow[0] = &lastBatchRowBuffer[0];
                lastBatchRow[1] = &lastBatchRowBuffer[graphLocal.numVertices];
              }

              //buffer to save read charactes for innermost loop
              std::vector<int32_t, aligned_allocator<int32_t, 64> > readCharsInt (SIMD_WIDTH * this->blockHeight);

              //process SIMD_WIDTH reads in a single iteration
#pragma omp for schedule(dynamic) nowait
              for (size_t i = 0; i < countReadBatches; i++)
              {
                __m512i bestScores512 = _ZERO;
                __m512i bestRows512   = _ZERO;
                __m512i bestCols512   = _ZERO;

                //reset DP 'lastBatchRow' buffer
                std::fill (lastBatchRowBuffer.begin(), lastBatchRowBuffer.end(), _ZERO);

                int32_t qryBatchLength = readSet[sortedReadOrder[i*SIMD_WIDTH]].length();  //longest read in this batch
                qryBatchLength += this->blockHeight - 1 - (qryBatchLength - 1) % this->blockHeight; //round-up

                //iterate over read length (process more than 1 characters in batch)
                for (int32_t j = 0; j < qryBatchLength; j += this->blockHeight)
                {
                  //loop counter 
                  size_t loopJ = j / (this->blockHeight);

                  //convert read character to int32_t
                  for (int32_t k = 0; k < SIMD_WIDTH * this->blockHeight ; k++)
                  {
                    readCharsInt [k] = readSetSOA [readSetSOAPrefixSum[i] + j*SIMD_WIDTH + k];
                  }

                  //iterate over characters in reference graph
                  for (int32_t k = 0; k < graphLocal.numVertices; k++)
                  {
                    //current reference character
                    __m512i graphChar = _SET1 ((int32_t) graphLocal.vertex_label[k] );

                    //current best score, init to 0
                    __m512i currentMax512;

                    //iterate over 'blockHeight' read characters
                    for (size_t l = 0; l < this->blockHeight; l++)
                    {
                      //load read characters
                      __m512i readChars = _LOAD ( &readCharsInt[l * SIMD_WIDTH] );

                      //current best score, init to 0
                      currentMax512 = _ZERO;

                      //see if query and reference character match
                      __mmask16 compareChar = _EQUAL (readChars, graphChar);
                      __m512i sub512 = _BLEND (compareChar, mismatch512, match512);

                      //match-mismatch edit
                      currentMax512 = _MAX (currentMax512, sub512); //local alignment can also start with a match at this char 

                      //iterate over graph neighbors
                      //which buffers to access depends on the value of 'l'
                      if (l == 0)
                      {
                        for(size_t m = graphLocal.offsets_in[k]; m < graphLocal.offsets_in[k+1]; m++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit = _ADD ( lastBatchRow[(loopJ - 1) & 1][ graphLocal.adjcny_in[m] ], sub512);
                          currentMax512 = _MAX (currentMax512, substEdit); 

                          //paths with deletion edit
                          __m512i delEdit;

                          if (k - graphLocal.adjcny_in[m] < this->blockWidth)
                            delEdit = _ADD ( nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][l], del512);
                          else
                            delEdit = _ADD ( fartherColumns[graphLocal.adjcny_in[m]][l], del512);

                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (lastBatchRow[(loopJ - 1) & 1][k], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }
                      else
                      {
                        for(size_t m = graphLocal.offsets_in[k]; m < graphLocal.offsets_in[k+1]; m++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit;

                          //paths with deletion edit
                          __m512i delEdit;

                          if (k - graphLocal.adjcny_in[m] < this->blockWidth)
                          {
                            substEdit = _ADD ( nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][l-1], sub512);
                            delEdit = _ADD ( nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][l], del512);
                          }
                          else
                          {
                            substEdit = _ADD ( fartherColumns[graphLocal.adjcny_in[m]][l-1], sub512);
                            delEdit = _ADD ( fartherColumns[graphLocal.adjcny_in[m]][l], del512);
                          }

                          currentMax512 = _MAX (currentMax512, substEdit); 
                          currentMax512 = _MAX (currentMax512, delEdit); 

                        }

                        //insertion edit
                        __m512i insEdit = _ADD (nearbyColumns[k & (blockWidth-1)][l-1], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }

                      //update best score observed yet
                      bestScores512 = _MAX (currentMax512, bestScores512);

                      //on which lanes is the best score updated
                      __mmask16 updated = _EQUAL (currentMax512, bestScores512);

                      //update row and column values accordingly
                      bestRows512 = _SET1_MASK (bestRows512, updated, (int32_t) (j + l));
                      bestCols512 = _SET1_MASK (bestCols512, updated, (int32_t) k);

                      //save current score in small buffer
                      nearbyColumns[k & (blockWidth-1)][l] = currentMax512;

                      //save current score in large buffer if connected thru long hop
                      if ( withLongHopLocal[k] )
                        fartherColumns[k][l] = currentMax512;
                    }

                    //save last score for next row-wise iteration
                    lastBatchRow[loopJ & 1][k] = currentMax512; 

                  } // end of row computation
                } // end of DP

                bestScores[i] = bestScores512;
                bestRows[i]   = bestRows512; 
                bestCols[i]   = bestCols512; 

              } // all reads done

              threadTimings[omp_get_thread_num()] = omp_get_wtime() - threadTimings[omp_get_thread_num()];

            } //end of omp parallel

            std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_vectorized, individual thread timings (s) : " << printStats(threadTimings) << "\n"; 

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

        //reference graph
        const CSR_char_container<VertexIdType, EdgeIdType> &graph;

        // pre-compute which graph vertices are connected with hop longer
        // than 'blockWidth'
        std::vector<bool> withLongHop;

        //for converting input reads into SOA to enable vectorization
        std::vector<char> readSetSOA;

        //cumulative read batch sizes
        std::vector<size_t>  readSetSOAPrefixSum;

        //input reads
        const std::vector<std::string> &readSet;

        //read lengths in their 'sorted order'
        std::vector<size_t> sortedReadLengths;

        //sorted permutation order of input reads
        std::vector<size_t> sortedReadOrder;

        //small temporary storage buffer for DP scores
        //should be a power of 2
        static constexpr size_t blockWidth = Phase1_Vectorized<ScoreType, VertexIdType, EdgeIdType>::blockWidth; 

        //process these many vertical cells in a go
        //should be a power of 2
        static constexpr size_t blockHeight = Phase1_Vectorized<ScoreType, VertexIdType, EdgeIdType>::blockHeight;   

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

          this->sortReadsForLoadBalance();
          this->convertToSOA();
          this->computeLongHops();
        };

        /**
         * @brief                                 wrapper function for phase 1 reverse DP routine 
         * @param[out]  outputBestScoreVector     vector to keep value and begin location of best scores,
         *                                        vector size is same as count of the reads
         * @note                                  this class won't care about rev. complement of sequences
         */
        template <typename Vec>
          void alignToDAGLocal_Phase1_rev_vectorized_wrapper(Vec &outputBestScoreVector) const
          {
            assert (outputBestScoreVector.size() == readSet.size());

            // modified containers to hold best-score info
            // vector elements aligned to 64 byte boundaries
            std::size_t countReadBatches = std::ceil (readSet.size() * 1.0 / SIMD_WIDTH);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreVector (countReadBatches);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreColVector (countReadBatches);
            std::vector <__m512i, aligned_allocator<__m512i, 64> > _bestScoreRowVector (countReadBatches);

            this->alignToDAGLocal_Phase1_rev_vectorized (outputBestScoreVector, _bestScoreVector, _bestScoreColVector, _bestScoreRowVector); 

            //parse best scores from vector registers
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeScores (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeCols   (SIMD_WIDTH);
            std::vector<int32_t, aligned_allocator<int32_t, 64> > storeRows   (SIMD_WIDTH);

            for (size_t i = 0; i < countReadBatches; i++)
            {
              _STORE (storeScores.data(), _bestScoreVector[i]);
              _STORE (storeCols.data(), _bestScoreColVector[i]);
              _STORE (storeRows.data(), _bestScoreRowVector[i]);

              for (size_t j = 0; j < SIMD_WIDTH; j++)
              {
                if (i*SIMD_WIDTH + j < readSet.size())
                {
                  auto originalReadId = sortedReadOrder[i*SIMD_WIDTH + j];

                  assert (outputBestScoreVector[originalReadId].score == storeScores[j] - 1);  //offset by 1

                  outputBestScoreVector[originalReadId].refColumnStart  = storeCols[j];
                  outputBestScoreVector[originalReadId].qryRowStart     = readSet[originalReadId].length() - 1 - storeRows[j];;

#ifdef DEBUG
                  std::cout << "INFO, psgl::Phase1_Vectorized::alignToDAGLocal_Phase1_rev_vectorized_wrapper, read # " << originalReadId << ",  score = " << storeScores[j] << "\n";
#endif
                }
              }
            }
          }

      private:

        /**
         * @brief       compute the sorted order of sequences for load balancing
         * @details     sorting is done in decreasing length order
         */
        void sortReadsForLoadBalance()
        {
          typedef std::pair<size_t, size_t> pair_t;

          //vector of tuples of length and index of reads
          std::vector<pair_t> lengthTuples;
           
          for(size_t i = 0; i < this->readSet.size(); i++)
            lengthTuples.emplace_back (readSet[i].length(), i);

          //sort in descending order, longer reads first
          std::sort (lengthTuples.begin(), lengthTuples.end(), std::greater<pair_t>() );

          for(auto &e : lengthTuples)
          {
            this->sortedReadLengths.push_back (e.first);
            this->sortedReadOrder.push_back (e.second);
          }
        }

        /**
         * @brief       convert input reads characters into SOA to enable vectorization
         */
        void convertToSOA()
        {
          assert (readSet.size() > 0);
          assert (sortedReadOrder.size() == readSet.size());
          assert (readSetSOA.size() == 0);

          /**
           * Requirements from the padding process
           * - each read length be a multiple of 'blockHeight'
           * - count of reads be a multiple of SIMD_WIDTH
           * - each batch of SIMD_WIDTH reads should be of equal length
           */

          //re-arrange read characters for vectorized processing; 
          //SIMD_WIDTH reads at a time, sorted in decreasing order by length
          
          auto readCount = readSet.size();

          readSetSOAPrefixSum.push_back(0);

          for (size_t i = 0; i < readCount; i += SIMD_WIDTH)
          {
            auto batchLength = readSet[sortedReadOrder[i]].length();  //longest read in this batch
            batchLength += blockHeight - 1 - (batchLength - 1) % blockHeight; //round-up

            for (size_t j = 0; j < batchLength; j++)
              for (size_t k = 0; k < SIMD_WIDTH; k++)
                if ( i + k < readCount && j < readSet[sortedReadOrder[i+k]].length() )
                  readSetSOA.push_back ( readSet[sortedReadOrder[i + k]][j] );
                else
                  readSetSOA.push_back (DUMMY);

            readSetSOAPrefixSum.push_back (readSetSOA.size());
          }

          std::size_t countReadBatches = std::ceil (readSet.size() * 1.0 / SIMD_WIDTH);
          assert (readSetSOAPrefixSum.size() == countReadBatches + 1);
        }

        /**
         * @brief     precompute which vertices are connected with long edge hops
         *            i.e. longer than 'blockWidth'
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

              //compare hop distance to 'blockWidth'
              if (to_pos - from_pos >= this->blockWidth)
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
            std::size_t readCount = readSet.size();
            std::size_t countReadBatches = std::ceil (readCount * 1.0 / SIMD_WIDTH);

            //few checks
            assert (bestScores.size() == countReadBatches);
            assert (bestCols.size() == countReadBatches);
            assert (bestRows.size() == countReadBatches);

//#ifdef VTUNE_SUPPORT
            //__itt_resume();
//#endif

            //init best score vector to zero bits
            std::fill (bestScores.begin(), bestScores.end(), _ZERO);
            std::fill (bestCols.begin(), bestCols.end(), _ZERO);
            std::fill (bestRows.begin(), bestRows.end(), _ZERO);

            //init score simd vectors
            __m512i match512    = _SET1 ((int32_t) SCORE::match);
            __m512i mismatch512 = _SET1 ((int32_t) SCORE::mismatch);
            __m512i del512      = _SET1 ((int32_t) SCORE::del);
            __m512i ins512      = _SET1 ((int32_t) SCORE::ins);

            std::vector<double> threadTimings (omp_get_max_threads(), 0);

#pragma omp parallel
            {
#pragma omp barrier
              threadTimings[omp_get_thread_num()] = omp_get_wtime();

              //create local copy of graph for faster access
              const CSR_char_container<VertexIdType, EdgeIdType> graphLocal = this->graph;
              const std::vector<bool> withLongHopLocal = withLongHop;

              //type def. for memory-aligned vector allocation for SIMD instructions
              using AlignedVecType = std::vector <__m512i, aligned_allocator<__m512i, 64> >;

              //buffer to save selected columns (associated with long hops) of DP matrix
              std::size_t countLongHops = std::count (withLongHopLocal.begin(), withLongHopLocal.end(), true);
              AlignedVecType fartherColumnsBuffer (countLongHops * this->blockHeight);

              //for convenient access to 2D buffer
              std::vector<__m512i*> fartherColumns(graphLocal.numVertices);
              {
                size_t j = 0;

                for(VertexIdType i = 0; i < graphLocal.numVertices; i++)
                  if ( this->withLongHop[i] )
                    fartherColumns[i] = &fartherColumnsBuffer[ (j++) * this->blockHeight ];
              }

              //buffer to save neighboring column scores
              AlignedVecType nearbyColumnsBuffer (this->blockWidth * this->blockHeight);

              //for convenient access to 2D buffer
              std::vector<__m512i*> nearbyColumns (this->blockWidth);
              {
                for (std::size_t i = 0; i < this->blockWidth; i++)
                  nearbyColumns[i] = &nearbyColumnsBuffer[i * this->blockHeight];
              }

              //buffer to save scores of last row in each iteration
              //one row for writing and one for reading
              AlignedVecType lastBatchRowBuffer (2 * graphLocal.numVertices);

              //for convenient access to 2D buffer
              std::vector<__m512i*> lastBatchRow (2);
              {
                lastBatchRow[0] = &lastBatchRowBuffer[0];
                lastBatchRow[1] = &lastBatchRowBuffer[graphLocal.numVertices];
              }

              //buffer to save read charactes for innermost loop
              std::vector<int32_t, aligned_allocator<int32_t, 64> > readCharsInt (SIMD_WIDTH * this->blockHeight);

              //process SIMD_WIDTH reads in a single iteration
#pragma omp for schedule(dynamic) nowait
              for (size_t i = 0; i < countReadBatches; i++)
              {
                __m512i fwdBestCols512;
                __m512i fwdBestRows512;

                //parse alignment locations of forward DP
                {
                  std::vector<int32_t, aligned_allocator<int32_t, 64> > fwdBestCols   (SIMD_WIDTH);
                  std::vector<int32_t, aligned_allocator<int32_t, 64> > fwdBestRows   (SIMD_WIDTH);

                  for (size_t j = 0; j < SIMD_WIDTH; j++)
                  {
                    if (i*SIMD_WIDTH + j < readSet.size())
                    {
                      auto originalReadId = sortedReadOrder[i*SIMD_WIDTH + j];
                      fwdBestCols[j] = outputBestScoreVector[originalReadId].refColumnEnd;
                      fwdBestRows[j] = readSet[originalReadId].length() - 1 - outputBestScoreVector[originalReadId].qryRowEnd;
                    }
                  }

                  fwdBestCols512 = _LOAD ( fwdBestCols.data() );
                  fwdBestRows512 = _LOAD ( fwdBestRows.data() );
                }

                __m512i bestScores512 = _ZERO;
                __m512i bestRows512   = _ZERO;
                __m512i bestCols512   = _ZERO;

                //reset DP 'lastBatchRow' buffer
                std::fill (lastBatchRowBuffer.begin(), lastBatchRowBuffer.end(), _ZERO);

                int32_t qryBatchLength = readSet[sortedReadOrder[i*SIMD_WIDTH]].length();  //longest read in this batch
                qryBatchLength += this->blockHeight - 1 - (qryBatchLength - 1) % this->blockHeight; //round-up

                //iterate over read length (process more than 1 characters in batch)
                for (int32_t j = 0; j < qryBatchLength; j += this->blockHeight)
                {
                  //loop counter 
                  size_t loopJ = j / (this->blockHeight);

                  //convert read character to int32_t
                  for (int32_t k = 0; k < SIMD_WIDTH * this->blockHeight ; k++)
                  {
                    readCharsInt [k] = readSetSOA [readSetSOAPrefixSum[i] + j*SIMD_WIDTH + k];
                  }

                  //iterate over characters in reference graph
                  for (int32_t k = graphLocal.numVertices - 1; k >= 0; k--)
                  {
                    //current reference character
                    __m512i graphChar = _SET1 ((int32_t) graphLocal.vertex_label[k] );

                    //current best score, init to 0
                    __m512i currentMax512;

                    //iterate over 'blockHeight' read characters
                    for (size_t l = 0; l < this->blockHeight; l++)
                    {
                      //load read characters
                      __m512i readChars = _LOAD ( &readCharsInt[l * SIMD_WIDTH] );

                      //current best score, init to 0
                      currentMax512 = _ZERO;

                      //see if query and reference character match
                      __mmask16 compareChar = _EQUAL (readChars, graphChar);
                      __m512i sub512 = _BLEND (compareChar, mismatch512, match512);

                      //match-mismatch edit
                      currentMax512 = _MAX (currentMax512, sub512); //local alignment can also start with a match at this char 

                      //iterate over graph neighbors
                      //which buffers to access depends on the value of 'l'
                      if (l == 0)
                      {
                        for(size_t m = graphLocal.offsets_out[k]; m < graphLocal.offsets_out[k+1]; m++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit = _ADD ( lastBatchRow[(loopJ - 1) & 1][ graphLocal.adjcny_out[m] ], sub512);
                          currentMax512 = _MAX (currentMax512, substEdit); 

                          //paths with deletion edit
                          __m512i delEdit;

                          if (graphLocal.adjcny_out[m] - k < this->blockWidth)
                            delEdit = _ADD ( nearbyColumns[graphLocal.adjcny_out[m] & (blockWidth-1)][l], del512);
                          else
                            delEdit = _ADD ( fartherColumns[graphLocal.adjcny_out[m]][l], del512);

                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (lastBatchRow[(loopJ - 1) & 1][k], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }
                      else
                      {
                        for(size_t m = graphLocal.offsets_out[k]; m < graphLocal.offsets_out[k+1]; m++)
                        {
                          //paths with match mismatch edit
                          __m512i substEdit;

                          //paths with deletion edit
                          __m512i delEdit;

                          if (graphLocal.adjcny_out[m] - k < this->blockWidth)
                          {
                            substEdit = _ADD ( nearbyColumns[graphLocal.adjcny_out[m] & (blockWidth-1)][l-1], sub512);
                            delEdit = _ADD ( nearbyColumns[graphLocal.adjcny_out[m] & (blockWidth-1)][l], del512);
                          }
                          else
                          {
                            substEdit = _ADD ( fartherColumns[graphLocal.adjcny_out[m]][l-1], sub512);
                            delEdit = _ADD ( fartherColumns[graphLocal.adjcny_out[m]][l], del512);
                          }

                          currentMax512 = _MAX (currentMax512, substEdit); 
                          currentMax512 = _MAX (currentMax512, delEdit); 
                        }

                        //insertion edit
                        __m512i insEdit = _ADD (nearbyColumns[k & (blockWidth-1)][l-1], ins512);
                        currentMax512 = _MAX (currentMax512, insEdit);
                      }

                      //update best score observed yet
                      bestScores512 = _MAX (currentMax512, bestScores512);

                      //on which lanes is the best score updated
                      __mmask16 updated = _EQUAL (currentMax512, bestScores512);

                      //update row and column values accordingly
                      bestRows512 = _SET1_MASK (bestRows512, updated, (int32_t) (j + l));
                      bestCols512 = _SET1_MASK (bestCols512, updated, (int32_t) k);

                      //detect and update score of the optimal alignment we want
                      {
                        __m512i currentRow = _SET1 ( (int32_t) (j + l));
                        __m512i currentCol = _SET1 ( (int32_t) k ); 
                        
                        __mmask16 compareCell = _EQUAL (fwdBestRows512, currentRow) & _EQUAL (fwdBestCols512, currentCol);

                        currentMax512 = _SET1_MASK (currentMax512, compareCell, (int32_t) (SCORE::match + 1)); 
                      }

                      //save current score in small buffer
                      nearbyColumns[k & (blockWidth-1)][l] = currentMax512;

                      //save current score in large buffer if connected thru long hop
                      if ( withLongHopLocal[k] )
                        fartherColumns[k][l] = currentMax512;
                    }

                    //save last score for next row-wise iteration
                    lastBatchRow[loopJ & 1][k] = currentMax512; 

                  } // end of row computation
                } // end of DP

                bestScores[i] = bestScores512;
                bestRows[i]   = bestRows512; 
                bestCols[i]   = bestCols512; 

              } // all reads done

              threadTimings[omp_get_thread_num()] = omp_get_wtime() - threadTimings[omp_get_thread_num()];

            } //end of omp parallel

            std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_Rev_Vectorized, individual thread timings (s) : " << printStats(threadTimings) << "\n"; 

//#ifdef VTUNE_SUPPORT
            //__itt_pause();
//#endif
          }
    };
}

#endif
