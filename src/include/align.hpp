/**
 * @file    align.hpp
 * @brief   routines to perform alignment
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


#ifndef GRAPH_ALIGN_HPP
#define GRAPH_ALIGN_HPP

#include <immintrin.h>

#include "graphLoad.hpp"
#include "csr_char.hpp"
#include "graph_iter.hpp"
#include "base_types.hpp"
#include "utils.hpp"
#include "align_vectorized.hpp"

//External includes
#include "kseq.h"
#include "prettyprint.hpp"

KSEQ_INIT(gzFile, gzread)

#define psgl_max(a,b) (((a)>(b))?(a):(b))

namespace psgl
{
  /**
   * @brief                         execute first phase of alignment i.e. compute DP and 
   *                                find locations of the best alignment of each read
   * @tparam[in]  ScoreType         type to store scores in DP matrix
   * @param[in]   readSet           vector of input query sequences to align
   * @param[in]   graph
   * @param[out]  bestScoreVector   vector to keep value and location of best scores,
   *                                vector size is same as count of the reads
   * @note                          reverse complement of the read is not handled here
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal_Phase1_scalar(  const std::vector<std::string> &readSet,
                                  const CSR_char_container<VertexIdType, EdgeIdType> &graph,
                                  std::vector< BestScoreInfo<ScoreType, VertexIdType> > &bestScoreVector)
    {
      assert (bestScoreVector.size() == readSet.size());

#ifdef VTUNE_SUPPORT
      __itt_resume();
#endif

      auto tick1 = __rdtsc();

#pragma omp parallel
      {

        //initialize matrix of size 2 x width, init with zero
        //we will keep re-using rows to keep memory-usage low
        std::vector< std::vector<ScoreType> > matrix(2, std::vector<ScoreType>(graph.numVertices, 0));

#pragma omp for
        for (size_t readno = 0; readno < readSet.size(); readno++)
        {
          //reset buffer
          std::fill(matrix[1].begin(), matrix[1].end(), 0);

          auto readLength = readSet[readno].length();

          //iterate over characters in read
          for (std::size_t i = 0; i < readLength; i++)
          {
            //iterate over characters in reference graph
            for (std::size_t j = 0; j < graph.numVertices; j++)
            {
              //current reference character
              char curChar = graph.vertex_label[j];

              ScoreType currentMax = 0;

              //see if query and ref. character match
              ScoreType matchScore = curChar == readSet[readno][i] ? SCORE::match : SCORE::mismatch;

              //match-mismatch edit
              currentMax = psgl_max (currentMax, matchScore);   //local alignment can also start with a match at this char

              for(auto k = graph.offsets_in[j]; k < graph.offsets_in[j+1]; k++)
              {
                //paths with match mismatch edit
                currentMax = psgl_max (currentMax, matrix[(i-1) & 1][ graph.adjcny_in[k] ] + matchScore);
                //'& 1' is same as doing modulo 2

                //paths with deletion edit
                currentMax = psgl_max (currentMax, matrix[i & 1][ graph.adjcny_in[k] ] + SCORE::del);
              }

              //insertion edit
              currentMax = psgl_max( currentMax, matrix[(i-1) & 1][j] + SCORE::ins );

              matrix[i & 1][j] = currentMax;

              //Update best score observed till now
              if (bestScoreVector[readno].score < matrix[i & 1][j])
              {
                bestScoreVector[readno].score = matrix[i & 1][j];
                bestScoreVector[readno].refColumn = j;
                bestScoreVector[readno].qryRow = i;
              }
            } // end of row computation
          } // end of DP

        } // all reads done
      } //end of omp parallel

      auto tick2 = __rdtsc();

      std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_scalar, CPU cycles spent in phase 1 = " << tick2 - tick1 
                << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";

#ifdef VTUNE_SUPPORT
        __itt_pause();
#endif

    }

  /**
   * @brief                         execute second phase of alignment i.e. compute cigar
   * @tparam[in]  ScoreType         type to store scores in DP matrix
   * @param[in]   readSet
   * @param[in]   graph
   * @param[in]   bestScoreVector   best score and alignment location for each read
   * @note                          we assume that query sequences are oriented properly
   *                                after executing the alignment phase 1
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal_Phase2(  const std::vector<std::string> &readSet,
                                  const CSR_char_container<VertexIdType, EdgeIdType> &graph,
                                  std::vector< BestScoreInfo<ScoreType, VertexIdType> > &bestScoreVector)
    {
      assert (bestScoreVector.size() == readSet.size());

      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        //for time profiling within phase 2
        uint64_t time_p2_1, time_p2_2, time_p2_3;

        //read length
        auto readLength = readSet[readno].length();

        //
        // PHASE 2.1 : COMPUTE FARTHEST REACHABLE VERTEX 
        //

        VertexIdType leftMostReachable;

        {
          auto tick1 = __rdtsc();

          std::size_t maxDistance = readLength + std::ceil( readLength * 1.0 * SCORE::match/SCORE::del );
          leftMostReachable = graph.computeLeftMostReachableVertex(bestScoreVector[readno].refColumn, maxDistance);  

#ifdef DEBUG
          std::cout << "INFO, psgl::alignToDAGLocal, left most reachable vertex id = " << leftMostReachable << std::endl;
#endif

          auto tick2 = __rdtsc();
          time_p2_1 = tick2 - tick1;
        }

        //
        // PHASE 2.2 : RECOMPUTE DP MATRIX WITH TRACEBACK INFORMATION
        //

        //width of score matrix that we need in memory
        std::size_t reducedWidth = bestScoreVector[readno].refColumn - leftMostReachable + 1;

        //for new beginning column
        std::size_t j0 = leftMostReachable; 

        //height of scoring matrix for re-computation
        std::size_t reducedHeight = bestScoreVector[readno].qryRow + 1;   //i.e. row id ass. with best score -- plus one 

        //scores in the last row
        std::vector<ScoreType> finalRow(reducedWidth, 0);

        //complete score matrix of size height x width to allow traceback
        //Note: to optimize storge, we only store vertical difference; absolute values of 
        //      which is bounded by gap penalty
        std::vector< std::vector<int8_t> > completeMatrixLog(reducedHeight, std::vector<int8_t>(reducedWidth, 0));

        {
          auto tick1 = __rdtsc();

          //scoring matrix of size 2 x width, init with zero
          std::vector<std::vector<ScoreType>> matrix(2, std::vector<ScoreType>(reducedWidth, 0));

          //iterate over characters in read
          for (std::size_t i = 0; i < reducedHeight; i++)
          {
            //iterate over characters in reference graph
            for (std::size_t j = 0; j < reducedWidth; j++)
            {
              //current reference character
              char curChar = graph.vertex_label[j + j0];

              //insertion edit
              ScoreType fromInsertion = matrix[(i-1) & 1][j] + SCORE::ins;
              //'& 1' is same as doing modulo 2

              //match-mismatch edit
              ScoreType matchScore = curChar == readSet[readno][i] ? SCORE::match : SCORE::mismatch;
              ScoreType fromMatch = matchScore;   //also handles the case when in-degree is zero 

              //deletion edit
              ScoreType fromDeletion  = -1; 

              for(auto k = graph.offsets_in[j + j0]; k < graph.offsets_in[j + j0 + 1]; k++)
              {
                //ignore edges outside the range 
                if ( graph.adjcny_in[k] >= j0)
                {
                  fromMatch = psgl_max (fromMatch, matrix[(i-1) & 1][ graph.adjcny_in[k] - j0] + matchScore);
                  fromDeletion = psgl_max (fromDeletion, matrix[i & 1][ graph.adjcny_in[k] - j0] + SCORE::del);
                }
              }

              //evaluate current score
              matrix[i & 1][j] = psgl_max ( psgl_max(fromInsertion, fromMatch) , psgl_max(fromDeletion, 0) );

              //save vertical difference of scores, used later for backtracking
              completeMatrixLog[i][j] = matrix[i & 1][j] - matrix[(i-1) & 1][j];
            }

            //Save last row
            if (i == reducedHeight - 1) 
              finalRow = matrix[i & 1];
          }

          ScoreType bestScoreReComputed = *std::max_element(finalRow.begin(), finalRow.end());

          //the recomputed score and its location should match our original calculation
          assert( bestScoreReComputed == bestScoreVector[readno].score );
          assert( bestScoreReComputed == finalRow[ bestScoreVector[readno].refColumn - j0 ] );

          auto tick2 = __rdtsc();
          time_p2_2 = tick2 - tick1;
        }

        //
        // PHASE 2.3 : COMPUTE CIGAR
        //
        
        std::string cigar;

        {
          auto tick1 = __rdtsc();

          std::vector<ScoreType> currentRowScores = finalRow; 
          std::vector<ScoreType> aboveRowScores (reducedWidth);

          int col = reducedWidth - 1;
          int row = bestScoreVector[readno].qryRow;

          while (col >= 0 && row >= 0)
          {
            if (currentRowScores[col] <= 0)
              break;

            //retrieve score values from vertical score differences
            for(std::size_t i = 0; i < reducedWidth; i++)
              aboveRowScores[i] = currentRowScores[i] - completeMatrixLog[row][i]; 

            //current reference character
            char curChar = graph.vertex_label[col + j0];

            //insertion edit
            ScoreType fromInsertion = aboveRowScores[col] + SCORE::ins;

            //match-mismatch edit
            ScoreType matchScore = curChar == readSet[readno][row] ? SCORE::match : SCORE::mismatch;

            ScoreType fromMatch = matchScore;   //also handles the case when in-degree is zero 
            std::size_t fromMatchPos = col;

            //deletion edit
            ScoreType fromDeletion = -1; 
            std::size_t fromDeletionPos;

            for(auto k = graph.offsets_in[col + j0]; k < graph.offsets_in[col + j0 + 1]; k++)
            {
              if ( graph.adjcny_in[k] >= j0)
              {
                auto fromCol = graph.adjcny_in[k] - j0;

                if (fromMatch < aboveRowScores[fromCol] + matchScore)
                {
                  fromMatch = aboveRowScores[fromCol] + matchScore;
                  fromMatchPos = fromCol;
                }

                if (fromDeletion < currentRowScores[fromCol] + SCORE::del)
                {
                  fromDeletion = currentRowScores[fromCol] + SCORE::del;
                  fromDeletionPos = fromCol;
                }
              }
            }

            //evaluate recurrence
            {
              if (currentRowScores[col] == fromMatch)
              {
                if (matchScore == SCORE::match)
                  cigar.push_back('=');
                else
                  cigar.push_back('X');

                //if alignment starts from this column, stop
                if (fromMatchPos == col)
                  break;

                //shift to preceeding column
                col = fromMatchPos;

                //shift to above row
                row--; currentRowScores = aboveRowScores;
              }
              else if (currentRowScores[col] == fromDeletion)
              {
                cigar.push_back('D');

                //shift to preceeding column
                col = fromDeletionPos;
              }
              else 
              {
                assert(currentRowScores[col] == fromInsertion);

                cigar.push_back('I');

                //shift to above row
                row--; currentRowScores = aboveRowScores;
              }
            }
          }

          //string reverse 
          std::reverse (cigar.begin(), cigar.end());  

          //shorten the cigar string
          psgl::seqUtils::cigarCompact(cigar);

          //validate if cigar yields best score
          assert ( psgl::seqUtils::cigarScore<ScoreType> (cigar) ==  bestScoreVector[readno].score );

          bestScoreVector[readno].cigar = cigar;

          auto tick2 = __rdtsc();
          time_p2_3 = tick2 - tick1;
        }

        std::cout << "INFO, psgl::alignToDAGLocal_Phase2, aligning read #" << readno + 1 << ", len = " << readLength << ", score " << bestScoreVector[readno].score << ", strand " << bestScoreVector[readno].strand << "\n";
        std::cout << "INFO, psgl::alignToDAGLocal_Phase2, cigar: " << bestScoreVector[readno].cigar << "\n";
        //std::cout << "TIMER, psgl::alignToDAGLocal_Phase2, CPU cycles spent in :  phase 2.1 = " << time_p2_1 << ", phase 2.2 = " << time_p2_2 << ", phase 2.3 = " << time_p2_3 << "\n";
        //std::cout.flush();
      }
    }

  /**
   * @brief                               local alignment routine
   * @tparam[in]  ScoreType               type to store scores in DP matrix
   * @param[in]   readSet
   * @param[in]   graph                   node-labeled directed graph 
   * @param[out]  outputBestScoreVector
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAGLocal( const std::vector<std::string> &readSet,
        const CSR_char_container<VertexIdType, EdgeIdType> &graph,
        std::vector< BestScoreInfo<ScoreType, VertexIdType> > &outputBestScoreVector)
    {
      //create buffer to save best score info for each read and its rev. complement
      std::vector< BestScoreInfo<ScoreType, VertexIdType> > bestScoreVector_P1 (2 * readSet.size() );

      //to save input query sequences for phase 1
      std::vector<std::string> readSet_P1;

      assert (readSet.size() > 0);

      //
      // Phase 1 [get best score values and location]
      //
      {
        auto tick1 = __rdtsc();

        for (size_t readno = 0; readno < readSet.size(); readno++)
        {
          std::string read_reverse (readSet[readno]);
          psgl::seqUtils::reverseComplement( readSet[readno], read_reverse); 

          readSet_P1.push_back (readSet[readno]);
          readSet_P1.push_back (read_reverse);
        }

        assert (bestScoreVector_P1.size() == 2 * readSet.size() );
        assert (readSet_P1.size() == 2 * readSet.size() );


        //align read to ref.
#ifdef __AVX512BW__
        alignToDAGLocal_Phase1_vectorized_wrapper(readSet_P1, graph, bestScoreVector_P1);
#else
        alignToDAGLocal_Phase1_scalar(readSet_P1, graph, bestScoreVector_P1);
#endif

        auto tick2 = __rdtsc();
        std::cout << "TIMER, psgl::alignToDAG, CPU cycles spent in phase 1  = " << tick2 - tick1
                  << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";
      }

      //
      // Phase 2 [comute cigar]
      //
      {
        auto tick1 = __rdtsc();

        std::vector<std::string> readSet_P2;

        for (size_t readno = 0; readno < readSet.size(); readno++)
        {
          if (bestScoreVector_P1[2 * readno].score > bestScoreVector_P1[2 * readno + 1].score)
          {
            outputBestScoreVector.push_back (bestScoreVector_P1[2 * readno]);
            readSet_P2.push_back (readSet_P1[2 * readno]);

            outputBestScoreVector[readno].strand = '+';
          }
          else
          {
            outputBestScoreVector.push_back (bestScoreVector_P1[2 * readno + 1]);
            readSet_P2.push_back (readSet_P1[2 * readno + 1]);

            outputBestScoreVector[readno].strand = '-';
          }
        }

        assert (outputBestScoreVector.size() == readSet.size() );
        assert (readSet_P2.size() == readSet.size() );

        alignToDAGLocal_Phase2(readSet_P2, graph, outputBestScoreVector);

        auto tick2 = __rdtsc();
        std::cout << "TIMER, psgl::alignToDAG, CPU cycles spent in phase 2  = " << tick2 - tick1
                  << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";
      }

    }

  /**
   * @brief                                 alignment routine
   * @tparam[in]  ScoreType                 type to store scores in DP matrix
   * @param[in]   reads                     vector of strings
   * @param[in]   graph
   * @param[out]  outputBestScoreVector
   * @param[in]   mode
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAG(const std::vector<std::string> &reads, 
        const CSR_char_container<VertexIdType, EdgeIdType> &graph,
        std::vector< BestScoreInfo<ScoreType, VertexIdType> > &outputBestScoreVector,
        const MODE mode)  
    {
      static_assert(std::is_signed<ScoreType>::value, 
          "ERROR, psgl::alignToDAG, ScoreType must be a signed type");

      switch(mode)
      {
        //case GLOBAL : alignToDAGGlobal<ScoreType> (reads, graph); break;
        case LOCAL : alignToDAGLocal<ScoreType> (reads, graph, outputBestScoreVector); break;
        //case SEMIGLOBAL: alignToDAGSemiGlobal<ScoreType> (reads, graph); break;
        default: std::cerr << "ERROR, psgl::alignToDAG, Invalid alignment mode"; exit(1);
      }
    }

  /**
   * @brief                                 alignment routine
   * @tparam[in]  ScoreType                 type to store scores in DP matrix
   * @param[in]   qfile                     file name containing reads
   * @param[in]   graph
   * @param[out]  outputBestScoreVector
   * @param[in]   mode
   */
  template <typename ScoreType, typename VertexIdType, typename EdgeIdType>
    void alignToDAG(const std::string &qfile, 
        const CSR_char_container<VertexIdType, EdgeIdType> &graph,
        std::vector< BestScoreInfo<ScoreType, VertexIdType> > &outputBestScoreVector,
        const MODE mode)  
    {
      //Parse all reads into a vector
      std::vector<std::string> reads;

      assert (outputBestScoreVector.empty());

      {
        if( !fileExists(qfile) )
        {
          std::cerr << qfile << " not accessible." << std::endl;
          exit(1);
        }

        //Open the file using kseq
        FILE *file = fopen (qfile.c_str(), "r");
        gzFile fp = gzdopen (fileno(file), "r");
        kseq_t *seq = kseq_init(fp);

        //size of sequence
        int len;

        while ((len = kseq_read(seq)) >= 0) 
        {
          psgl::seqUtils::makeUpperCase(seq->seq.s, len);

          reads.push_back(seq->seq.s);
        }

        //Close the input file
        kseq_destroy(seq);  
        gzclose(fp);  
      }

      std::cout << "INFO, psgl::alignToDAG, total count of reads = " << reads.size() << std::endl;

      alignToDAG<ScoreType> (reads, graph, outputBestScoreVector, mode);
    }
}

#endif
