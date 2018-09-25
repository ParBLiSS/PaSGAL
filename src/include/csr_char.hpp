/**
 * @file    csr_char.hpp
 * @brief   routines to store graph in CSR format 
 *          (single character per vertex)
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CSR_CHAR_CONTAINER_HPP
#define CSR_CHAR_CONTAINER_HPP

#include <cassert>
#include <iostream>

//Own includes
#include "csr.hpp"
#include "graph_iter.hpp"

namespace psgl
{

  /**
   * @brief     class to support storage of directed graph in CSR format
   *            each vertex holds a character label
   * @details   CSR is a popular graph storage format-
   *            - vertex numbering starts from 0
   *            - Adjacency list (incoming edges) of vertex i is stored in 
   *              array 'adjcny_in' starting at index offsets_in[i] and 
   *              ending at (but not including) index offsets_in[i+1]
   *            - CSR_char_container should be built using existing 
   *              CSR_container (see class constructor)
   *
   *            - Assumption: count of vertices and edges < 2B because we are
   *              using int32_t type everywhere
   */
  class CSR_char_container
  {
    public:

      //Count of edges and vertices in the graph
      int32_t numVertices;
      int32_t numEdges;

      //contiguous adjacency list of all vertices, size = numEdges
      std::vector<int32_t> adjcny_in;  
      std::vector<int32_t> adjcny_out;  

      //offsets in adjacency list for each vertex, size = numVertices + 1
      std::vector<int32_t> offsets_in;
      std::vector<int32_t> offsets_out;

      //Container to hold character label of all vertices in graph
      std::vector<char> vertex_label;


      /**
       * @brief             build complete CSR_char graph
       * @param[in]   csr   CSR graph container
       */
      void build (CSR_container &csr)
      {
        this->numVertices = csr.totalRefLength();

        vertex_label.reserve (csr.totalRefLength());

        adjcny_in.reserve (csr.totalRefLength() + csr.numEdges - csr.numVertices);
        adjcny_out.reserve (csr.totalRefLength() + csr.numEdges - csr.numVertices);

        offsets_in.reserve (csr.totalRefLength() + 1);
        offsets_out.reserve (csr.totalRefLength() + 1);

        //Save vertex labels
        {
          for (auto &seq : csr.vertex_metadata) 
            for (auto &c : seq)
              vertex_label.push_back(c);
        }

        //Init edges:
        {
          // in edges
          {
            offsets_in.push_back(0);
            std::vector<int32_t> inNeighbors;

            for (graphIterFwd g(csr); !g.end(); g.next())
            {
              //get preceeding dependency offsets from graph
              inNeighbors.clear();
              g.getInNeighborOffsets(inNeighbors);

              for(auto &e : inNeighbors)
                adjcny_in.push_back(e);

              offsets_in.push_back (adjcny_in.size());
            }
          }

          // out edges
          {
            offsets_out.push_back(0);
            std::vector<int32_t> outNeighbors;

            for (graphIterFwd g(csr); !g.end(); g.next())
            {
              //get preceeding dependency offsets from graph
              outNeighbors.clear();
              g.getOutNeighborOffsets(outNeighbors);

              for(auto &e : outNeighbors)
                adjcny_out.push_back(e);

              offsets_out.push_back (adjcny_out.size());
            }
          }
        }

        this->numEdges = adjcny_in.size();

        assert(vertex_label.size() == this->numVertices);

        assert(adjcny_in.size() == csr.totalRefLength() + csr.numEdges - csr.numVertices);
        assert(adjcny_out.size() == csr.totalRefLength() + csr.numEdges - csr.numVertices);

        assert(offsets_in.size() ==  csr.totalRefLength() + 1);
        assert(offsets_out.size() ==  csr.totalRefLength() + 1);

#ifndef NDEBUG
        this->verify();
#endif

        std::cout << "INFO, psgl::CSR_char_container::build, graph converted to CSR format with character labels, n = " << this->numVertices << ", m = " << this->numEdges << std::endl;
      }

      /**
       * @brief             compute and print histogram of degree values
       */
      void printDegreeHistogram () const
      {
        std::cout << "Printing degree distribution ..." << "\n"; 

        int32_t maxDegree = 0;

        //compute maximum degree in the graph
        for(int32_t i = 0; i < this->numVertices; i++)
          for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
            maxDegree = std::max (maxDegree, offsets_in[i+1] - offsets_in[i]);

        std::vector<int32_t> degreeHist (maxDegree + 1, 0);

        //compute histogram
        for(int32_t i = 0; i < this->numVertices; i++)
          for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
            degreeHist [offsets_in[i+1] - offsets_in[i] ]++; 

        for(int32_t i = 0; i <= maxDegree; i++)
          if (degreeHist[i] > 0)
            std::cout << i << " : " << degreeHist[i] << "\n"; 

        std::cout.flush();
      }

      /**
       * @brief             compute and print histogram of hop distances in 
       *                    the sorted order
       */
      void printHopLengthHistogram() const
      {
        std::cout << "Printing hop length distribution  ..." << "\n"; 

        //get maximum hop length
        int32_t maxHopLength = this->directedBandwidth();

        std::vector<int32_t> hopLengthHist (maxHopLength + 1, 0);

        //compute histogram
        for(int32_t i = 0; i < this->numVertices; i++)
          for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
            hopLengthHist [ i - adjcny_in[j] ] ++;

        for(int32_t i = 0; i <= maxHopLength; i++)
          if (hopLengthHist[i] > 0)
            std::cout << i << " : " << hopLengthHist[i] << "\n"; 

        std::cout.flush();
      }

    private:

      /**
       * @brief                   compute maximum distance between connected vertices in the graph (a.k.a. 
       *                          directed bandwidth)
       * @return                  directed graph bandwidth
       */
      std::size_t directedBandwidth() const
      {
        std::size_t bandwidth = 0;   //temporary value 

        //iterate over all vertices in graph to compute bandwidth
        for(int32_t i = 0; i < this->numVertices; i++)
        {
          for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
          {
            auto from_pos = adjcny_in[j];
            auto to_pos = i;

            assert(to_pos > from_pos);

            if (to_pos - from_pos > bandwidth)
              bandwidth = to_pos - from_pos;
          }
        }

        return bandwidth;
      }

      /**
       * @brief     sanity check for correctness of graph storage in CSR format
       */
      void verify() const
      {
        assert(this->numVertices > 0);
        assert(this->numEdges > 0);

        //labels
        {
          assert(vertex_label.size() == this->numVertices);

          for(auto &c : vertex_label)
          {
            assert(std::isupper(c));
            assert(c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'N');
          }
        }

        //adjacency list
        {
          assert(adjcny_in.size() == this->numEdges);
          assert(adjcny_out.size() == this->numEdges);

          for(auto vId : adjcny_in)
            assert(vId >=0 && vId < this->numVertices);

          for(auto vId : adjcny_out)
            assert(vId >=0 && vId < this->numVertices);
        }

        //offset array
        {
          assert(offsets_in.size() == this->numVertices + 1);
          assert(std::is_sorted(offsets_in.begin(), offsets_in.end()));
          assert(offsets_in.front() == 0 && offsets_in.back() == this->numEdges);

          assert(offsets_out.size() == this->numVertices + 1);
          assert(std::is_sorted(offsets_out.begin(), offsets_out.end()));
          assert(offsets_out.front() == 0 && offsets_out.back() == this->numEdges);

          for(auto off : offsets_out)
            assert(off >=0 && off <= this->numEdges); 
        }

        //topologically sorted order
        {
          for(int32_t i = 0; i < this->numVertices; i++)
            for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
              assert( adjcny_in[j] < i);

          for(int32_t i = 0; i < this->numVertices; i++)
            for(auto j = offsets_out[i]; j < offsets_out[i+1]; j++)
              assert( adjcny_out[j] > i);
        }
      }
  };
}

#endif
