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
   * @brief                   local alignment routine
   * @tparam[in]  ScoreType   type to store scores in DP matrix
   * @param[in]   reads
   * @param[in]   graph
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal( const std::vector<std::string> &reads,
                          const CSR_container<VertexIdType, EdgeIdType> &graph)
    {
      //total chars in reference
      auto referenceLength = graph.totalRefLength();

      //dimension of score matrix that we need in memory
      auto breadth = graph.directedBandwidth() + 1;

      //iterate over reads
      for (auto &read : reads)
      {
        auto height = read.length() + 1;

        ScoreType best_score = 0;

        //initialize matrix of size height x breadth, init with zero
        std::vector<std::vector<ScoreType> > matrix(breadth, std::vector<ScoreType>(height, 0));

        //iterate over reference graph
        for (graphIter<VertexIdType, EdgeIdType> g(graph); !g.end(); g.next())
        {
          //current reference character
          char curChar = g.curChar();

          //current column in DP matrix
          std::size_t i = g.getGlobalOffset() % breadth;

          //get preceeding dependency offsets from graph
          std::vector<std::size_t> preceedingOffsets;
          g.getNeighborOffsets(preceedingOffsets);

          for (int j = 1; j < height; j++)
          {
            //deletion
            ScoreType fromDeletion = matrix[i][j - 1] - SCORE::del;

            //match-mismatch
            auto matchScore = curChar == read.at(j-1) ? SCORE::match : -1 * SCORE::mismatch;

            ScoreType fromMatch = matchScore;   //handles the case when in-degree is zero 
            for(auto k : preceedingOffsets)
            {
              fromMatch = std::max (fromMatch, matrix[(k) % breadth][j - 1] + matchScore);
            }
             

            ScoreType fromInsertion = 0; 
            for(auto k : preceedingOffsets)
            {
              fromInsertion = std::max (fromInsertion, matrix[(k) % breadth][j] - SCORE::ins);
            }

            matrix[i][j] = std::max( std::max(fromDeletion, fromMatch) , std::max(fromInsertion, 0) );
            best_score = std::max(best_score, matrix[i][j]);
          }
        }

        std::cout << "INFO, psgl::alignToDAGLocal, best score = " << best_score << std::endl;
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
