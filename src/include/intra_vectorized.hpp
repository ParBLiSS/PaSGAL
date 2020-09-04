/**
 * @file    align_vectorized.hpp
 * @file    intra_vectorized.hpp    
 * @brief   vectorized routines (intra-task) to perform alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ALIGN_VEC_HPP
#define GRAPH_ALIGN_VEC_HPP

#include <immintrin.h>
#include <x86intrin.h>

#include "graphLoad.hpp"
#include "csr_char.hpp"
#include "graph_iter.hpp"
#include "base_types.hpp"
#include "utils.hpp"

//External includes

#define DUMMY                     'B'

#if defined(PASGAL_ENABLE_AVX512)
  #define SIMD_REG_SIZE           512 
  typedef __m512i __mxxxi;          
#elif defined(PASGAL_ENABLE_AVX2)
  #define SIMD_REG_SIZE           256 
  typedef __m256i __mxxxi;   
#endif

namespace psgl
{
  /**
   * Parameters and SIMD instructions specific for different score precisions required
   */

  template<typename T> struct SimdInst {};

  template<>
    struct SimdInst<int32_t> 
    {
      typedef int32_t type;

#if defined(PASGAL_ENABLE_AVX512)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));

#elif defined(PASGAL_ENABLE_AVX2)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));
#endif

    };

  template<>
    struct SimdInst<int16_t> 
    {
      typedef int16_t type;

#if defined(PASGAL_ENABLE_AVX512)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));

#elif defined(PASGAL_ENABLE_AVX2)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));
#endif

    };

  template<>
    struct SimdInst<int8_t> 
    {
      typedef int8_t type;

#if defined(PASGAL_ENABLE_AVX512)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));

#elif defined(PASGAL_ENABLE_AVX2)
      static constexpr int p = SIMD_REG_SIZE / (8 * sizeof(type));
#endif

    };


  static const int8_t kBaseTranslation[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    // A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //           T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    // a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    //           t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
  };

  /**
   * @brief   Supports phase 1 DP in forward direction
   */
  template <typename SIMD>
    class Phase1_Vectorized
    {
      private:

        //reference graph
        const CSR_char_container &graph;

        // pre-compute which graph vertices are connected with hop longer
        // than 'blockWidth'
        std::vector<bool> withLongHop;

        //input reads
        const std::vector<std::string> &readSet;

        //input parameters (e.g., scoring scheme)
        const Parameters &parameters;

      public:

        //small temporary storage buffer for DP scores
        //should be a power of 2
        static constexpr size_t blockWidth = 8 ; 

        //process these many vertical cells in a go
        //should be a power of 2
        static constexpr size_t blockHeight = 16 * SIMD::p;   

        /**
         * @brief                   public constructor
         * @param[in]   readSet     vector of input query sequences to align
         * @param[in]   g           input reference graph
         */
        Phase1_Vectorized(const std::vector<std::string> &readSet, 
            const CSR_char_container &g,
            const Parameters &p) :
          readSet (readSet), graph (g), parameters (p)
        {
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

            for (size_t i = 0; i < readCount; i += 1)
            {
              int _bestScore = -1, _bestScoreCol = -1, _bestScoreRow = -1;

              // execute the alignment routine
              this->alignToDAGLocal_Phase1_vectorized (&readSet[i], _bestScore, _bestScoreCol, _bestScoreRow); 

              assert (_bestScore > 0  && _bestScoreCol > 0 && _bestScoreRow > 0);

              outputBestScoreVector[i].score         = _bestScore;
              outputBestScoreVector[i].refColumnEnd  = _bestScoreCol;
              outputBestScoreVector[i].qryRowEnd     = _bestScoreRow;
            }
          }

      private:

        /**
         * @brief     precompute which vertices are connected with long edge hops
         *            i.e. longer than 'blockWidth'
         */
        void computeLongHops()
        {
          assert (withLongHop.size() == 0);

          this->withLongHop.resize(graph.numVertices, false);

          for(int32_t i = 0; i < graph.numVertices; i++)
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
         *  @brief    return string whose length is a multiple of blockHeight
         */
        void getPaddedRead(const std::string &in, std::string &out)
        {
          std::size_t originalLength = in.length();
          assert (originalLength > 0);

          std::size_t newLength = (originalLength / this->blockHeight) * this->blockHeight;

          //copy
          out = in;

          //append DUMMY characters
          if (newLength > originalLength)
            out.append(newLength - originalLength, DUMMY);

          assert (out.length() == newLength);
          assert (out.length () % this->blockHeight == 0);
        }

        template <typename Vec>
          void buildQueryProfile (Vec &queryProfileBuffer, const std::string &read, std::size_t beg, std::size_t end)
          {
            assert (end - beg == this->blockHeight);

            constexpr int t = this->blockHeight / SIMD::p;

            //first fill scalar buffer, we will load these values into SIMD vector afterwards
            // initialize the buffer with mismatch score
            std::vector<typename SIMD::type, aligned_alloc<typename SIMD::type, 64> > queryProfileScalarBuffer (4 * this->blockHeight, (typename SIMD::type) -1 * parameters.mismatch);

            //for convenient access to 2D buffer
            std::vector<SIMD::type*> queryProfileScalar (4);
            for (int i = 0; i < 4; i++)
              queryProfileScalar = &queryProfileScalarBuffer[i * this->blockHeight];

            //set appropriate values to match parameter
            for (int i = 0; i < t; i++)
              for (int j = 0; j < SIMD::p; j++)
              {
                if (read[beg + t*j + i] != DUMMY)
                  queryProfileScalar[kBaseTranslation(read[beg + t*j + i])][i*SIMD::p + j] = parameters.match;
              }

            //load queryProfileBuffer from queryProfileScalarBuffer
            assert (queryProfileScalarBuffer.size() == queryProfileBuffer.size() * SIMD::p);
            for (std::size_t i = 0; i < queryProfileBuffer.size(); i += 1)
              queryProfileBuffer[i] = SIMD::load ((const __mxxxi*) &queryProfileScalarBuffer [i * SIMD::p] ); 
          }

        void alignToDAGLocal_Phase1_vectorized (const std::string &read, int &_bestScore, int &_bestScoreCol, int &_bestScoreRow) const
        {
          std::string _read;
          
          getPaddedRead (read, _read);
          assert (_read.length() % blockHeight == 0);
          assert (blockHeight % SIMD::p == 0);

          //type def. for memory-aligned vector allocation for SIMD instructions
          using AlignedVecType = std::vector <__mxxxi, aligned_alloc<__mxxxi, 64> >;

          __mxxxi del512      = SIMD::set1 ((typename SIMD::type) -1 * parameters.del);
          __mxxxi ins512      = SIMD::set1 ((typename SIMD::type) -1 * parameters.ins);

          AlignedVecType queryProfileBuffer (4 * this->blockHeight / SIMD::p);
          //alphabet size = 4
          //for convenient access to above 2D buffer
          std::vector<__mxxxi*> queryProfile (4);
          {
            for (std::size_t i = 0; i < 4; i++)
              queryProfile[i] = &queryProfileBuffer[i * this->blockHeight / SIMD::p];
          } 

          //buffer to save scores of last row in each iteration
          //one row for writing and one for reading
          std::vector<typename SIMD::type, aligned_alloc<typename SIMD::type, 64> > lastBatchRowBuffer (2 * graphLocal.numVertices, 0);

          //temporary buffer to save all scalar values from any SIMD register
          std::vector<typename SIMD::type, aligned_alloc<typename SIMD::type, 64> > scalarVals (SIMD::p);

          //buffer to save neighboring column scores
          AlignedVecType nearbyColumnsBuffer (this->blockWidth * this->blockHeight / SIMD::p);
          //for convenient access to 2D buffer
          std::vector<__mxxxi*> nearbyColumns (this->blockWidth);
          {
            for (std::size_t i = 0; i < this->blockWidth; i++)
              nearbyColumns[i] = &nearbyColumnsBuffer[i * this->blockHeight / SIMD::p];
          }

          //2D buffer to save selected columns (associated with long hops) of DP matrix
          std::size_t countLongHops = std::count (withLongHop.begin(), withLongHop.end(), true);
          AlignedVecType fartherColumnsBuffer (countLongHops * this->blockHeight / SIMD::p);
          //pointer array for convenient data access in 2D buffer
          std::vector<__mxxxi*> fartherColumns(graphLocal.numVertices);
          {
            size_t j = 0;
            for(int32_t i = 0; i < graphLocal.numVertices; i++)
              if ( this->withLongHop[i] )
                fartherColumns[i] = &fartherColumnsBuffer[ (j++) * this->blockHeight / SIMD::p];
          }


          //for-loop for blockHeight-sized groups of rows
          for (std::size_t i = 0; i < _read.length()/this->blockHeight; i++)
          {
            //buffer to save query profile
            int read_pos_begin = i * this->blockHeight; //inclusive
            int read_pos_end = read_pos_begin + this->blockHeight; //exclusive
            buildQueryProfile (queryProfile, _read, read_pos_begin, read_pos_end);

            int8_t firstReadCharacter = kBaseTranslation (_read[read_pos_begin]);

            //iterate over characters in reference graph
            for (int32_t j = 0; j < graphLocal.numVertices; j++)
            {
              int8_t graphChar = kBaseTranslation (graphLocal.vertex_label[j]);  
              
              //iterate over SIMD::p -sized groups of rows
              for (std::size_t k = 0; k < this->blockHeight / SIMD::p; k++)
              {
                //current best score, init to 0
                __mxxxi currentMax512 = SIMD::zero();

                __mxxxi substEdit, correctedSubstEdit, rightShifted, delEdit, insEdit;

                if (k == 0)
                {

                  //load information from previous row
                  SIMD::type scoreFirstCell = std::max(0, lastBatchRow[(i-1) & 1][j] - parameters.ins);

                  //check match of read character at read_pos_begin
                  //use last row buffer to compute subst. score
                  SIMD::type matchReward = (graphChar == firstReadCharacter) ? parameters.match : -1 * parameters.mismatch;
                  if (j)
                    scoreFirstCell = std::max(scoreFirstCell, lastBatchRow[(i-1) & 1][j - 1] + matchReward);
                  else
                    scoreFirstCell = std::max(scoreFirstCell, matchReward);

                  //iterate over vertex neighbors in the graph
                  for(size_t m = graphLocal.offsets_in[j]; m < graphLocal.offsets_in[j+1]; m++)
                  {
                    if (j - graphLocal.adjcny_in[m] < this->blockWidth)
                    {
                      rightShifted = SIMD::rightshift (nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][this->blockHeight/SIMD::p  -1] );
                      substEdit = SIMD::add (rightShifted, queryProfile[graphChar][k]); 
                      delEdit = SIMD::add (nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][k], del512);
                    }
                    else
                    {
                      rightShifted = SIMD::rightshift (fartherColumns[graphLocal.adjcny_in[m]][this->blockHeight/SIMD::p  -1] );
                      substEdit = SIMD::add (rightShifted, queryProfile[graphChar][k]); 
                      delEdit = SIMD::add (fartherColumns[graphLocal.adjcny_in[m]][k], del512);
                    }
                    
                    correctedSubstEdit = SIMD::setFirst (substEdit, scoreFirstCell);
                    currentMax512 = SIMD::max (currentMax512, correctedSubstEdit); 
                    currentMax512 = SIMD::max (currentMax512, delEdit); 
                  }
                  //assume no insertion edit propagation here
                }
                else
                {
                  //iterate over vertex neighbors in the graph
                  for(size_t m = graphLocal.offsets_in[j]; m < graphLocal.offsets_in[j+1]; m++)
                  {
                    if (j - graphLocal.adjcny_in[m] < this->blockWidth)
                    {
                      substEdit = SIMD::add (nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][k-1], queryProfile[graphChar][k]); 
                      delEdit = SIMD::add (nearbyColumns[graphLocal.adjcny_in[m] & (blockWidth-1)][k], del512);
                    }
                    else
                    {
                      substEdit = SIMD::add (fartherColumns[graphLocal.adjcny_in[m]][k-1], queryProfile[graphChar][k]); 
                      delEdit = SIMD::add (fartherColumns[graphLocal.adjcny_in[m]][k], del512);
                    }

                    currentMax512 = SIMD::max (currentMax512, substEdit); 
                    currentMax512 = SIMD::max (currentMax512, delEdit); 
                  }

                  insEdit = SIMD::add (nearbyColumns[j & (blockWidth-1)][k-1], ins512); 
                  currentMax512 = SIMD::max (currentMax512, insEdit);
                }

                nearbyColumns[j & (blockWidth-1)][k] = currentMax512;
              }

              //check and run lazy f-loop
              rightShifted = SIMD::rightshift (nearbyColumns[j & (blockWidth-1)][this->blockHeight/SIMD::p  -1] );
              insEdit = SIMD::add (rightShifted, ins512);
              auto needToRecompute = SIMD::cmpgt (insEdit, nearbyColumns[j & (blockWidth-1)][0]);
              while (needToRecompute > 0)
              {
                nearbyColumns[j & (blockWidth-1)][0] = SIMD::max (insEdit, nearbyColumns[j & (blockWidth-1)][0]);
                
                for (std::size_t k = 1; k < this->blockHeight / SIMD::p; k++)
                {
                  insEdit = SIMD::add (nearbyColumns[j & (blockWidth-1)][k-1], ins512); 
                  nearbyColumns[j & (blockWidth-1)][k] = SIMD::max (insEdit, nearbyColumns[j & (blockWidth-1)][k-1]);
                }

                rightShifted = SIMD::rightshift (nearbyColumns[j & (blockWidth-1)][this->blockHeight/SIMD::p  -1] );
                insEdit = SIMD::add (rightShifted, ins512);
                needToRecompute = SIMD::cmpgt (insEdit, nearbyColumns[j & (blockWidth-1)][0]);
              }

              //save last-row scores into row-buffer
              SIMD::store ((__mxxxi*) scalarVals.data(), nearbyColumns[j & (blockWidth-1)][this->blockHeight / SIMD::p  -1] );
              lastBatchRow[i & 1][j] = scalarVals[SIMD::p  -1];

              
              //save currentMax into bestScore
              for (std::size_t k = 0; k < this->blockHeight / SIMD::p; k++)
              {
                __mxxxi score = nearbyColumns[j & (blockWidth-1)][k];  
                
                //search for maximum scalar value inside 'score' register
                int newBest = std::max (_bestScore, SIMD::reduceMax (score));

                //update
                if (newBest > _bestScore)
                {
                  //score
                  _bestScore = newBest;

                  //column
                  _bestScoreCol = j;

                  //row
                  __mxxxi bestScoreVec = SIMD::set1 ( (typename SIMD::type) (_bestScore));
                  auto cmpeq_mask = SIMD::cmpeq (score, bestScoreVec);
                  int rowIndex =  __builtin_ctz(cmpeq_mask);
                  _bestScoreRow = i * (_read.length()/this->blockHeight) + (this>blockWidth / SIMD::p) * (rowIndex % SIMD::p)   + (rowIndex / SIMD::p);
                }
              }

            }
          }

        }

    };
}

#endif
