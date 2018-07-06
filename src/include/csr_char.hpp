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
   */
  template <typename VertexIdType, typename EdgeIdType>
    class CSR_char_container
    {
      public:

        //Count of edges and vertices in the graph
        VertexIdType numVertices;
        EdgeIdType numEdges;

        //contiguous adjacency list of all vertices, size = numEdges
        std::vector<VertexIdType> adjcny_in;  

        //offsets in adjacency list for each vertex, size = numVertices + 1
        std::vector<EdgeIdType> offsets_in;

        //Container to hold character label of all vertices in graph
        std::vector<char> vertex_label;


        /**
         * @brief             build complete CSR_char graph
         * @param[in]   csr   CSR graph container
         */
        void build (CSR_container<VertexIdType, EdgeIdType> &csr)
        {
          this->numVertices = csr.totalRefLength();
          
          vertex_label.reserve (csr.totalRefLength());
          adjcny_in.reserve (csr.totalRefLength() + csr.numEdges - csr.numVertices);
          offsets_in.reserve (csr.totalRefLength() + 1);

          //Save vertex labels
          {
            for (auto &seq : csr.vertex_metadata) 
              for (auto &c : seq)
                vertex_label.push_back(c);

          }

          //Init edges:
          {
            offsets_in.push_back(0);

            for (graphIterFwd <VertexIdType, EdgeIdType> g(csr); !g.end(); g.next())
            {
              //get preceeding dependency offsets from graph
              std::vector<VertexIdType> inNeighbors;
              g.getNeighborOffsets(inNeighbors);

              for(auto &e : inNeighbors)
                adjcny_in.push_back(e);

              offsets_in.push_back (adjcny_in.size());
            }
          }

          this->numEdges = adjcny_in.size();

          assert(vertex_label.size() == this->numVertices);
          assert(adjcny_in.size() == csr.totalRefLength() + csr.numEdges - csr.numVertices);
          assert(offsets_in.size() ==  csr.totalRefLength() + 1);

#ifndef NDEBUG
          this->verify();
#endif

          std::cout << "INFO, psgl::CSR_char_container::build, graph converted to CSR format with character labels, n = " << this->numVertices << ", m = " << this->numEdges << std::endl;
        }

        /**
         * @brief             compute farthest reachable vertex on left side
         * @param[in]   v     id of starting vertex (where alignment ends)
         * @param[in]   dist  alignment length bound
         * @details           useful to determining traceback range in DP
         */
        VertexIdType computeLeftMostReachableVertex(VertexIdType v, std::size_t dist) const
        {
          assert(v >= 0 && v < this->numVertices);

          //trivial case
          if (v == 0)
            return 0;

          //Initialize a distance vector for vertices {0, 1,... , v-1, v}
          std::vector<std::size_t> distvec(v+1, dist);

          //Initialize distance of starting vertex
          distvec[v] = 0;

          //Traverse vertices in reverse of topological sorted order
          for (VertexIdType i = v; i != (VertexIdType) -1; i--)
          {
            //Update in-neigbors of vertex i
            for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
            {
              assert( adjcny_in[j] < i );

              distvec[ adjcny_in[j] ] = std::min( distvec[ adjcny_in[j] ] , distvec[i] + 1 );
            }
          }

          //Check the leftmost vertex inside dist
          for (VertexIdType i = 0; i <= v; i++)
          {
            if (distvec[i] < dist)
              return i;
          }

          return 0;
        }

      private:


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

            for(auto vId : adjcny_in)
              assert(vId >=0 && vId < this->numVertices);
          }

          //offset array
          {
            assert(offsets_in.size() == this->numVertices + 1);
            assert(std::is_sorted(offsets_in.begin(), offsets_in.end()));
            assert(offsets_in.front() == 0 && offsets_in.back() == this->numEdges);

            for(auto off : offsets_in)
              assert(off >=0 && off <= this->numEdges); 
          }

          //topologically sorted order
          {
            for(VertexIdType i = 0; i < this->numVertices; i++)
              for(auto j = offsets_in[i]; j < offsets_in[i+1]; j++)
                assert( adjcny_in[j] < i);
          }
        }
    };
}

#endif
