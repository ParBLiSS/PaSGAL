/**
 * @file    align.hpp
 * @brief   routines to perform alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ALIGN_HPP
#define GRAPH_ALIGN_HPP

#include "graphLoad.hpp"
#include "csr.hpp"
#include "graph_iter.hpp"
#include "base_types.hpp"

namespace psgl
{

  /**
   * @brief                   container to save info about best score
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   */
  template <typename ScoreType, typename VertexIdType>
    struct BestScoreInfo
    {
      //score value
      ScoreType score;

      //vertex id and sequence offset within where optimal alignment ends
      VertexIdType vid;
      std::size_t vertexSeqOffset;

      //positioning in complete DP matrix where optimal alignment ends
      std::size_t refColumn;
      std::size_t qryRow;

      /**
       * @brief   constructor
       */
      BestScoreInfo()
      {
        this->score = 0;
      }
    };

  /**
   * @brief                   local alignment routine
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   * @param[in]   reads
   * @param[in]   graph
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal( const std::vector<std::string> &reads,
        const CSR_container<VertexIdType, EdgeIdType> &graph)
    {
      //width of score matrix that we need in memory
      auto width = graph.totalRefLength();

      //iterate over reads
      for (auto &read : reads)
      {
        std::size_t height = read.length();

        /**
         * PHASE 1 : COMPUTE COMPLETE DP MATRIX
         */

        BestScoreInfo<ScoreType, VertexIdType> best;

        {
          //initialize matrix of size 2 x width, init with zero
          //we will keep re-using rows to keep memory-usage low
          std::vector<std::vector<ScoreType>> matrix(2, std::vector<ScoreType>(width, 0));

          for (std::size_t i = 0; i < height; i++)
          {
            //iterate over characters in reference graph
            for (graphIter<VertexIdType, EdgeIdType> g(graph); !g.end(); g.next())
            {
              //current reference character
              char curChar = g.curChar();

              //current column number in DP matrix
              std::size_t j = g.getGlobalOffset();

              //deletion
              ScoreType fromDeletion = matrix[(i-1) % 2][j] - SCORE::del;

              //get preceeding dependency offsets from graph
              std::vector<std::size_t> preceedingOffsets;
              g.getNeighborOffsets(preceedingOffsets);

              //match-mismatch
              ScoreType matchScore = curChar == read[i] ? SCORE::match : -1 * SCORE::mismatch;

              ScoreType fromMatch = matchScore;   //handles the case when in-degree is zero 
              for(auto k : preceedingOffsets)
              {
                fromMatch = std::max (fromMatch, matrix[(i-1) % 2][k] + matchScore);
              }

              //insertion
              ScoreType fromInsertion = 0; 
              for(auto k : preceedingOffsets)
              {
                fromInsertion = std::max (fromInsertion, matrix[i % 2][k] - SCORE::ins);
              }

              //Evaluate current score
              matrix[i % 2][j] = std::max( std::max(fromDeletion, fromMatch) , std::max(fromInsertion, 0) );

              //Update best score observed till now
              if (best.score < matrix[i % 2][j])
              {
                best.score = matrix[i % 2][j];
                best.vid = g.getCurrentVertexId();
                best.vertexSeqOffset = g.getCurrentSeqOffset();
                best.refColumn = j;
                best.qryRow = i;
              }
            }
          }

          std::cout << "INFO, psgl::alignToDAGLocal, best score = " << best.score << ", with alignment ending at vertex id = " << best.vid << std::endl;
        }   

        /**
         * PHASE 2 : COMPUTE FARTHEST REACHABLE VERTEX 
         */

        VertexIdType leftMostReachable;

        {
          //assume that there is atmost 10% insertion error rate
          std::size_t maxDistance = read.length() + std::ceil( read.length() * 1.0/10 );
          leftMostReachable = graph.computeLeftMostReachableVertex(best.vid, maxDistance);  
          std::cout << "INFO, psgl::alignToDAGLocal, left most reachable vertex id = " << leftMostReachable << std::endl;
        }

        /**
         * PHASE 3 : RECOMPUTE DP MATRIX WITH TRACEBACK INFORMATION
         */

        //width of score matrix that we need in memory
        auto reducedWidth = graph.totalRefLength(leftMostReachable, best.vid);

        auto reducedHeight = best.qryRow + 1;   //i.e. row id ass. with best score -- plus one 

        //scores in the last row
        std::vector<ScoreType> finalRow(reducedWidth, 0);

        //complete matrix of size height x width
        //Note: to optimize storge, we only store vertical difference; absolute values of which is bounded by gap penalty
        std::vector<std::vector<int8_t>> completeMatrixLog(reducedHeight, std::vector<int8_t>(reducedWidth, 0));

        {
          //scoring matrix of size 2 x width, init with zero
          std::vector<std::vector<ScoreType>> matrix(2, std::vector<ScoreType>(reducedWidth, 0));

          for (std::size_t i = 0; i < reducedHeight; i++)
          {
            graphIter<VertexIdType, EdgeIdType> g(graph, leftMostReachable);
            std::size_t j0 = g.getGlobalOffset();   //beginning column

            //iterate over characters in reference graph
            for (std::size_t j = 0; j < reducedWidth; j++)
            {
              //current reference character
              char curChar = g.curChar();

              //deletion
              ScoreType fromDeletion = matrix[(i-1) % 2][j] - SCORE::del;

              //get preceeding dependency offsets from graph
              std::vector<std::size_t> preceedingOffsets;
              g.getNeighborOffsets(preceedingOffsets);

              //match-mismatch
              ScoreType matchScore = curChar == read[i] ? SCORE::match : -1 * SCORE::mismatch;

              ScoreType fromMatch = matchScore;   //handles the case when in-degree is zero 
              for(auto k : preceedingOffsets)
              {
                fromMatch = std::max (fromMatch, matrix[(i-1) % 2][k-j0] + matchScore);
              }

              //insertion
              ScoreType fromInsertion = 0; 
              for(auto k : preceedingOffsets)
              {
                fromInsertion = std::max (fromInsertion, matrix[i % 2][k-j0] - SCORE::ins);
              }

              //evaluate current score
              matrix[i % 2][j] = std::max( std::max(fromDeletion, fromMatch) , std::max(fromInsertion, 0) );

              //save vertical difference of scores
              completeMatrixLog[i][j] = matrix[i % 2][j] - matrix[(i-1) % 2][j];

              //advance graph iterator
              g.next();
            }

            //Save last row
            if (i == reducedHeight - 1) 
              finalRow = matrix[i % 2];
          }
        }

        /**
         * PHASE 3.1 : VERIFY CORRECTNESS OF RE-COMPUTE
         */

        {
          graphIter<VertexIdType, EdgeIdType> g(graph, leftMostReachable);
          std::size_t j0 = g.getGlobalOffset();   //beginning column

          ScoreType bestScoreReComputed = finalRow[ best.refColumn - j0 ];

          assert(bestScoreReComputed == best.score);
        }

        /**
         * PHASE 4 : COMPUTE CIGAR
         */

        {
          for (std::size_t i = 0; i < reducedHeight; i++)
          {
            //iterate over characters in reference graph
            for (std::size_t j = 0; j < reducedWidth; j++)
            {
               
            }
          }
        }

      }
    }

  /**
   * @brief                   global alignment routine
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   * @param[in]   reads
   * @param[in]   graph
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGGlobal( const std::vector<std::string> &reads,
                          const CSR_container<VertexIdType, EdgeIdType> &graph)
    {
    }

  /**
   * @brief                   semi-global alignment routine
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   * @param[in]   reads
   * @param[in]   graph
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGSemiGlobal( const std::vector<std::string> &reads,
        const CSR_container<VertexIdType, EdgeIdType> &graph)
    {
    }

  /**
   * @brief                   alignment routine
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   * @param[in]   reads
   * @param[in]   graph
   * @param[in]   mode
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAG(const std::vector<std::string> &reads, 
        const CSR_container<VertexIdType, EdgeIdType> &graph,
        const MODE mode)  
    {
      static_assert(std::is_signed<ScoreType>::value, 
          "ERROR, psgl::alignToDAG, ScoreType must be a signed type");

      switch(mode)
      {
        case GLOBAL : alignToDAGGlobal<ScoreType> (reads, graph); break;
        case LOCAL : alignToDAGLocal<ScoreType> (reads, graph); break;
        case SEMIGLOBAL: alignToDAGSemiGlobal<ScoreType> (reads, graph); break;
        default: std::cerr << "ERROR, psgl::alignToDAG, Invalid alignment mode"; exit(1);
      }
    }

}

#endif
