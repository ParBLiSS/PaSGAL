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
   * @param[in]   readSet           vector of input query sequences to align
   * @param[in]   graph
   * @param[out]  bestScoreVector   vector to keep value and location of best scores,
   *                                vector size is same as count of the reads
   * @note                          reverse complement of the read is not handled here
   */
  void alignToDAGLocal_Phase1_scalar( const std::vector<std::string> &readSet,
                                      const CSR_char_container &graph,
                                      std::vector< BestScoreInfo > &bestScoreVector)
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
      std::vector< std::vector<int32_t> > matrix(2, std::vector<int32_t>(graph.numVertices, 0));

#pragma omp for schedule(dynamic)
      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        //reset buffer
        std::fill(matrix[1].begin(), matrix[1].end(), 0);

        auto readLength = readSet[readno].length();

        int32_t bestScore = 0;
        int32_t bestRow = 0, bestCol = 0;

        //iterate over characters in read
        for (int32_t i = 0; i < readLength; i++)
        {
          //iterate over characters in reference graph
          for (int32_t j = 0; j < graph.numVertices; j++)
          {
            //current reference character
            char curChar = graph.vertex_label[j];

            int32_t currentMax = 0;

            //see if query and ref. character match
            int32_t matchScore = curChar == readSet[readno][i] ? SCORE::match : SCORE::mismatch;

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

            bestScore = psgl_max (bestScore, currentMax);

            //Update best score observed till now
            if (bestScore == currentMax)
            {
              bestScore = currentMax; bestCol = j; bestRow = i;
            }
          } // end of row computation
        } // end of DP

        bestScoreVector[readno].score = bestScore;
        bestScoreVector[readno].refColumnEnd = bestCol;
        bestScoreVector[readno].qryRowEnd = bestRow;

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
   * @brief                         execute first phase of alignment in reverse direction (right to left)
   *                                i.e. compute reverse DP so as to find
   *                                find begin location of the best alignment of each read
   * @param[in]   readSet           vector of input query sequences to align
   * @param[in]   graph
   * @param[out]  bestScoreVector   vector to keep value and location of best scores,
   *                                vector size is same as count of the reads
   * @note                          reverse complement of the read is not handled here
   */
  void alignToDAGLocal_Phase1_rev_scalar( const std::vector<std::string> &readSet,
                                          const CSR_char_container &graph,
                                          std::vector< BestScoreInfo > &bestScoreVector)
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
      std::vector< std::vector<int32_t> > matrix(2, std::vector<int32_t>(graph.numVertices, 0));

#pragma omp for schedule(dynamic)
      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        //reset buffer
        std::fill(matrix[1].begin(), matrix[1].end(), 0);

        auto readLength = readSet[readno].length();

        int32_t bestScore = 0;
        int32_t bestRow = 0, bestCol = 0;

        //iterate over characters in read
        for (int32_t i = 0; i < readLength; i++)
        {
          //iterate over characters in reference graph
          for (int32_t j = graph.numVertices - 1; j >= 0; j--)
          {
            //current reference character
            char curChar = graph.vertex_label[j];

            int32_t currentMax = 0;

            //see if query and ref. character match
            int32_t matchScore = curChar == readSet[readno][i] ? SCORE::match : SCORE::mismatch;

            //match-mismatch edit
            currentMax = psgl_max (currentMax, matchScore);   //local alignment can also start with a match at this char

            for(auto k = graph.offsets_out[j]; k < graph.offsets_out[j+1]; k++)
            {
              //paths with match mismatch edit
              currentMax = psgl_max (currentMax, matrix[(i-1) & 1][ graph.adjcny_out[k] ] + matchScore);
              //'& 1' is same as doing modulo 2

              //paths with deletion edit
              currentMax = psgl_max (currentMax, matrix[i & 1][ graph.adjcny_out[k] ] + SCORE::del);
            }

            //insertion edit
            currentMax = psgl_max( currentMax, matrix[(i-1) & 1][j] + SCORE::ins );

            matrix[i & 1][j] = currentMax;

            bestScore = psgl_max (bestScore, currentMax);

            //Update best score observed till now
            if (bestScore == currentMax)
            {
              bestScore = currentMax; bestCol = j; bestRow = readLength - 1 - i;
            }

            //special handling of the cell where optimal alignment had ended during forward DP
            if (j == bestScoreVector[readno].refColumnEnd && (readLength - 1 - i) == bestScoreVector[readno].qryRowEnd)
            {
              //local alignment needs to end with a match
              assert (currentMax == SCORE::match);

              //add one so that the other end of the optimal alignment can be located without ambuiguity
              matrix[i & 1][j] = SCORE::match + 1;
            }
          } // end of row computation
        } // end of DP

        assert (bestScoreVector[readno].score == bestScore - 1);    //offset by 1
        bestScoreVector[readno].refColumnStart = bestCol;
        bestScoreVector[readno].qryRowStart = bestRow;

      } // all reads done
    } //end of omp parallel

    auto tick2 = __rdtsc();

    std::cout << "TIMER, psgl::alignToDAGLocal_Phase1_rev_scalar, CPU cycles spent in phase 1-R = " << tick2 - tick1 
      << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";

#ifdef VTUNE_SUPPORT
    __itt_pause();
#endif

  }

  /**
   * @brief                         execute second phase of alignment i.e. compute cigar
   * @param[in]   readSet
   * @param[in]   graph
   * @param[in]   bestScoreVector   best score and alignment location for each read
   * @note                          we assume that query sequences are oriented properly
   *                                after executing the alignment phase 1
   */
  void alignToDAGLocal_Phase2(  const std::vector<std::string> &readSet,
                                const CSR_char_container &graph,
                                std::vector< BestScoreInfo > &bestScoreVector)
  {
    assert (bestScoreVector.size() == readSet.size());

    //not using more than 16 threads
    int maxThreadCount = omp_get_max_threads() <= 16 ? omp_get_max_threads() : 16;

    std::vector<double> threadTimings (maxThreadCount, 0);

#pragma omp parallel num_threads(maxThreadCount)  
    {
      threadTimings[omp_get_thread_num()] = omp_get_wtime();

#pragma omp for schedule(dynamic) nowait
      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        //for time profiling within phase 2
        uint64_t time_p2_1, time_p2_2;

        //read length
        auto readLength = readSet[readno].length();

        //
        // PHASE 2.1 : RECOMPUTE DP MATRIX WITH TRACEBACK INFORMATION
        // recomputation is done within selected block of DP matrix
        //

        //width of score matrix that we need in memory
        std::size_t reducedWidth = bestScoreVector[readno].refColumnEnd - bestScoreVector[readno].refColumnStart + 1;

        //for new beginning column
        std::size_t j0 = bestScoreVector[readno].refColumnStart; 

        //height of scoring matrix for re-computation
        std::size_t reducedHeight = bestScoreVector[readno].qryRowEnd - bestScoreVector[readno].qryRowStart + 1; 

        //for new beginning row
        std::size_t i0 = bestScoreVector[readno].qryRowStart; 

        //scores in the last row
        std::vector<int32_t> finalRow(reducedWidth, 0);

        //complete score matrix of size height x width to allow traceback
        //Note: to optimize storge, we only store vertical difference; absolute values of 
        //      which is bounded by gap penalty
#ifdef DEBUG
        std::cout << "INFO, psgl::alignToDAGLocal_Phase2, aligning read #" << readno + 1 << ", memory requested= " << reducedWidth * reducedHeight << " bytes" << std::endl;
#endif

        std::vector< std::vector<int8_t> > completeMatrixLog(reducedHeight, std::vector<int8_t>(reducedWidth, 0));

        {
          auto tick1 = __rdtsc();

          //scoring matrix of size 2 x width, init with zero
          std::vector< std::vector<int32_t> > matrix(2, std::vector<int32_t>(reducedWidth, 0));

          //iterate over characters in read
          for (std::size_t i = 0; i < reducedHeight; i++)
          {
            //iterate over characters in reference graph
            for (std::size_t j = 0; j < reducedWidth; j++)
            {
              //current reference character
              char curChar = graph.vertex_label[j + j0];

              //insertion edit
              int32_t fromInsertion = matrix[(i-1) & 1][j] + SCORE::ins;
              //'& 1' is same as doing modulo 2

              //match-mismatch edit
              int32_t matchScore = curChar == readSet[readno][i + i0] ? SCORE::match : SCORE::mismatch;
              int32_t fromMatch = matchScore;   //also handles the case when in-degree is zero 

              //deletion edit
              int32_t fromDeletion  = -1; 

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

          int32_t bestScoreReComputed = *std::max_element(finalRow.begin(), finalRow.end());

          //the recomputed score and its location should match our original calculation
          assert( bestScoreReComputed == bestScoreVector[readno].score );
          assert( bestScoreReComputed == finalRow[ bestScoreVector[readno].refColumnEnd - j0 ] );

          auto tick2 = __rdtsc();
          time_p2_1 = tick2 - tick1;
        }

        //
        // PHASE 2.2 : COMPUTE CIGAR
        //

        std::string cigar;

        {
          auto tick1 = __rdtsc();

          std::vector<int32_t> currentRowScores = finalRow; 
          std::vector<int32_t> aboveRowScores (reducedWidth);

          int col = reducedWidth - 1;
          int row = reducedHeight - 1;

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
            int32_t fromInsertion = aboveRowScores[col] + SCORE::ins;

            //match-mismatch edit
            int32_t matchScore = curChar == readSet[readno][row + i0] ? SCORE::match : SCORE::mismatch;

            int32_t fromMatch = matchScore;   //also handles the case when in-degree is zero 
            std::size_t fromMatchPos = col;

            //deletion edit
            int32_t fromDeletion = -1; 
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
          assert ( psgl::seqUtils::cigarScore (cigar) ==  bestScoreVector[readno].score );

          bestScoreVector[readno].cigar = cigar;

          auto tick2 = __rdtsc();
          time_p2_2 = tick2 - tick1;
        }

#ifdef DEBUG
        std::cout << "INFO, psgl::alignToDAGLocal_Phase2, aligning read #" << readno + 1 << ", len = " << readLength << ", score " << bestScoreVector[readno].score << ", strand " << bestScoreVector[readno].strand << "\n";
        std::cout << "INFO, psgl::alignToDAGLocal_Phase2, cigar: " << bestScoreVector[readno].cigar << "\n";
        std::cout << "TIMER, psgl::alignToDAGLocal_Phase2, CPU cycles spent in :  phase 2.1 = " << time_p2_1 * 1.0 / ASSUMED_CPU_FREQ << ", phase 2.2 = " << time_p2_2 * 1.0 / ASSUMED_CPU_FREQ << "\n";
        //std::cout.flush();
#endif
      }

      threadTimings[omp_get_thread_num()] = omp_get_wtime() - threadTimings[omp_get_thread_num()];

    }

    std::cout << "TIMER, psgl::alignToDAGLocal_Phase2, individual thread timings (s) : " << printStats(threadTimings) << "\n"; 
  }

  /**
   * @brief                               local alignment routine
   * @param[in]   readSet
   * @param[in]   graph                   node-labeled directed graph 
   * @param[out]  outputBestScoreVector
   */
  void alignToDAGLocal( const std::vector<std::string> &readSet,
      const CSR_char_container &graph,
      std::vector< BestScoreInfo > &outputBestScoreVector)
  {
    //create buffer to save best score info for each read and its rev. complement
    std::vector< BestScoreInfo > bestScoreVector_P1 (2 * readSet.size() );

    //to save input query sequences for phase 1
    std::vector<std::string> readSet_P1;

    assert (readSet.size() > 0);
    assert (outputBestScoreVector.empty());

    //TODO: Dynamically decide the SIMD word type
    typedef SimdInst<int32_t> SIMD;

    //
    // Phase 1 [get best score values and locations]
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
      Phase1_Vectorized<SIMD> obj (readSet_P1, graph); 
      obj.alignToDAGLocal_Phase1_vectorized_wrapper(bestScoreVector_P1);
#else
      alignToDAGLocal_Phase1_scalar (readSet_P1, graph, bestScoreVector_P1);
#endif

      auto tick2 = __rdtsc();
      std::cout << "TIMER, psgl::alignToDAG, CPU cycles spent in phase 1  = " << tick2 - tick1
        << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << std::endl;
    }

#ifdef DEBUG
    {
      for (size_t readno = 0; readno < readSet_P1.size(); readno++)
      {
        std::cout << "INFO, psgl::alignToDAGLocal, read # " << readno + 1 << ", score = " << bestScoreVector_P1[readno].score 
          << ", refColumnEnd = " << bestScoreVector_P1[readno].refColumnEnd
          << ", qryRowEnd = " << bestScoreVector_P1[readno].qryRowEnd << "\n";
      }
    }
#endif

    //
    // Phase 1 - reverse [get begin locations]
    //
    {
      auto tick1 = __rdtsc();

      std::vector<std::string> readSet_P1_R;

      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        if (bestScoreVector_P1[2 * readno].score > bestScoreVector_P1[2 * readno + 1].score)
        {
          outputBestScoreVector.push_back (bestScoreVector_P1[2 * readno]);
          outputBestScoreVector[readno].strand = '+';

          std::string read_reverse (readSet_P1[2 * readno]);
          psgl::seqUtils::reverse( readSet_P1[2 * readno], read_reverse); 
          readSet_P1_R.push_back (read_reverse);
        }
        else
        {
          outputBestScoreVector.push_back (bestScoreVector_P1[2 * readno + 1]);
          outputBestScoreVector[readno].strand = '-';

          std::string read_reverse (readSet_P1[2 * readno + 1]);
          psgl::seqUtils::reverse( readSet_P1[2 * readno + 1], read_reverse); 
          readSet_P1_R.push_back (read_reverse);
        }
      }

      assert (outputBestScoreVector.size() == readSet.size() );
      assert (readSet_P1_R.size() == readSet.size() );

      //align reverse read to ref.
#ifdef __AVX512BW__
      Phase1_Rev_Vectorized<SIMD> obj (readSet_P1_R, graph); 
      obj.alignToDAGLocal_Phase1_rev_vectorized_wrapper(outputBestScoreVector);
#else
      alignToDAGLocal_Phase1_rev_scalar (readSet_P1_R, graph, outputBestScoreVector);
#endif

      auto tick2 = __rdtsc();
      std::cout << "TIMER, psgl::alignToDAG, CPU cycles spent in phase 1-R  = " << tick2 - tick1
        << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << std::endl;
    }

#ifdef DEBUG
    {
      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        std::cout << "INFO, psgl::alignToDAGLocal, read # " << readno + 1 << ", score = " << outputBestScoreVector[readno].score 
          << ", refColumnStart = " << outputBestScoreVector[readno].refColumnStart 
          << ", refColumnEnd = " << outputBestScoreVector[readno].refColumnEnd
          << ", qryRowStart = " << outputBestScoreVector[readno].qryRowStart
          << ", qryRowEnd = " << outputBestScoreVector[readno].qryRowEnd << "\n";
      }
    }
#endif

    //
    // Phase 2 [comute cigar]
    //
    {
      auto tick1 = __rdtsc();

      std::vector<std::string> readSet_P2;

      for (size_t readno = 0; readno < readSet.size(); readno++)
      {
        if (bestScoreVector_P1[2 * readno].score > bestScoreVector_P1[2 * readno + 1].score)
          readSet_P2.push_back (readSet_P1[2 * readno]);
        else
          readSet_P2.push_back (readSet_P1[2 * readno + 1]);
      }

      assert (readSet_P2.size() == readSet.size() );

      alignToDAGLocal_Phase2 (readSet_P2, graph, outputBestScoreVector);

      auto tick2 = __rdtsc();
      std::cout << "TIMER, psgl::alignToDAG, CPU cycles spent in phase 2  = " << tick2 - tick1
        << ", estimated time (s) = " << (tick2 - tick1) * 1.0 / ASSUMED_CPU_FREQ << "\n";
    }
  }

  /**
   * @brief                                 alignment routine
   * @param[in]   reads                     vector of strings
   * @param[in]   graph
   * @param[out]  outputBestScoreVector
   * @param[in]   mode
   */
    void alignToDAG(  const std::vector<std::string> &reads, 
                      const CSR_char_container &graph,
                      std::vector< BestScoreInfo > &outputBestScoreVector,
                      const MODE mode)  
    {
      //TODO: Support other alignment modes: global and semi-global
      switch(mode)
      {
        case LOCAL : alignToDAGLocal (reads, graph, outputBestScoreVector); break;
        default: std::cerr << "ERROR, psgl::alignToDAG, Invalid alignment mode"; exit(1);
      }
    }

  /**
   * @brief                                 alignment routine
   * @param[in]   qfile                     file name containing reads
   * @param[in]   graph
   * @param[out]  outputBestScoreVector
   * @param[in]   mode
   */
    int alignToDAG( const std::string &qfile, 
                    const CSR_char_container &graph,
                    std::vector< BestScoreInfo > &outputBestScoreVector,
                    const MODE mode)  
    {
      //Parse all reads into a vector
      std::vector<std::string> reads;

      assert (outputBestScoreVector.empty());

      {
        //TODO: Read query sequences in batches rather than all at once
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

      alignToDAG (reads, graph, outputBestScoreVector, mode);

      return PSGL_STATUS_OK;
    }

    /**
     * @brief     print alignment results to file
     */
    void printResultsToFile ( const std::string &outputFile,
        const std::vector< BestScoreInfo > &outputBestScoreVector)
    {
      std::ofstream outstrm(outputFile);

      for(auto &e : outputBestScoreVector)
      {
        outstrm << e.score << "\t"
          << e.strand << "\t"
          << e.qryRowStart << "\t" 
          << e.qryRowEnd << "\t"
          << e.refColumnStart << "\t"
          << e.refColumnEnd << "\t"
          << e.cigar << "\n";
      }
    }
}

#endif
