/**
 * @file    csr.hpp
 * @brief   routines to store graph in CSR format 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CSR_CONTAINER_HPP
#define CSR_CONTAINER_HPP

#include <cassert>
#include <iostream>

//Own includes
#include "stream.hpp"
#include "vg.pb.h"

namespace psgl
{

  /**
   * @brief     class to support storage of directed graph in CSR format
   * @details   CSR is a popular graph storage format-
   *            - vertex numbering starts from 0
   *            - Adjacency list of vertex i is stored in array 'adjcny'
   *              starting at index offsets[i] and ending at (but not including)
   *              index offsets[i+1]
   */
  template <typename VertexIdType, typename EdgeIdType>
    class CSR_container
    {
      public:

        //Count of edges and vertices in the graph
        VertexIdType numVertices;
        EdgeIdType numEdges;

        //contiguous adjacency list of all vertices, size = numEdges
        std::vector<VertexIdType> adjcny;  

        //offsets in adjacency list for each vertex, size = numVertices + 1
        std::vector<EdgeIdType> offsets;

        //Container to hold metadata (e.g., DNA sequence) of all vertices in graph
        std::vector<std::string> vertex_metadata;

        /**
         * @brief     constructor
         */
        CSR_container()
        {
          numVertices = 0;
          numEdges = 0;
        }

        /**
         * @brief     sanity check for correctness of graph storage in CSR format
         */
        void verify()
        {
          //sequences
          {
            assert(vertex_metadata.size() == this->numVertices);

            for(auto &seq : vertex_metadata)
              assert(seq.length() > 0);
          }

          //adjacency list
          {
            assert(adjcny.size() == this->numEdges);

            for(auto vId : adjcny)
              assert(vId >=0 && vId < this->numVertices);
          }

          //offset array
          {
            assert(offsets.size() == this->numVertices + 1);

            for(auto off : offsets)
              assert(off >=0 && off <= this->numEdges); 

            assert(std::is_sorted(offsets.begin(), offsets.end()));

            assert(offsets.back() == this->numEdges);
          }
        }

        /**
         * @brief           add vertices in graph
         * @param[in]   n   count of vertices to add
         */
        void addVertexCount(VertexIdType n)
        {
          assert(n > 0);  

          this->numVertices += n;
          this->vertex_metadata.resize(this->numVertices);
        }

        /**
         * @brief             add vertex sequence in the graph
         * @param[in]   id    vertex id
         * @param[in]   seq   sequence
         */
        void initVertexSequence(VertexIdType id, const std::string &seq)
        {
          assert(id >= 0 && id < vertex_metadata.size());  
          assert(vertex_metadata[id].length() == 0);

          vertex_metadata[id] = seq;
        }

        /**
         * @brief                   add edges from a vector
         * @param[in] diEdgeVector  vector containing directed edges as pairs <from, to> 
         */
        void initEdges(std::vector<std::pair<VertexIdType, VertexIdType>> &diEdgeVector)
        {
          for(auto &edge: diEdgeVector)
          {
            assert(edge.first >= 0 && edge.first < this->numVertices);
            assert(edge.second >=0 && edge.second < this->numVertices);
          }

          this->numEdges = diEdgeVector.size();
          adjcny.reserve(this->numEdges);
          offsets.resize(this->numVertices + 1);

          //Sort the edge vector before adding it into adjcny
          std::sort(diEdgeVector.begin(), diEdgeVector.end());

          //Compute offsets and adjcny vectors
          offsets[0] = 0;
          auto it_b = diEdgeVector.begin();

          for(VertexIdType i = 0; i < this->numVertices; i++)
          {
            //Range for adjacency list of vertex i
            auto it_e = std::find_if(it_b, diEdgeVector.end(), [i](std::pair<VertexIdType, VertexIdType> &e) { return e.first > i; });

            offsets[i+1] = std::distance(diEdgeVector.begin(), it_e);

            for(auto it = it_b; it != it_e; it++)
              adjcny.push_back(it->second);

            it_b = it_e;   
          }
        }

        /**
         * @brief     print the loaded directed graph to stderr 
         * @details   Format details (assuming n = no. of vertices):
         *              Print n+1 rows in total
         *              First row: value of n and total count of edges
         *              Subsequent rows contain out-neighbors of vertices and their DNA sequences,
         *              one row per vertex
         *            This function is implemented for debugging purpose
         */
        void printGraph()
        {
          std::cerr << this->numVertices << " " << this->numEdges << "\n";

          for (VertexIdType i = 0; i < this->numVertices; i++)
          {
            for (EdgeIdType j = offsets[i]; j < offsets[i+1]; j++)
              std::cerr << adjcny[j] << " ";

            std::cerr << vertex_metadata[i] << "\n";
          }
        }
    };

}

#endif
