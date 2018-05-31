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
   *            - both outgoing and incoming edges are store separately 
   *              - this is redundant but convenient for analysis
   *              - Adjacency (out/in) list of vertex i is stored in array 
   *                'adjcny_' starting at index offsets_[i] and 
   *                ending at (but not including) index offsets_[i+1]
   */
  template <typename VertexIdType, typename EdgeIdType>
    class CSR_container
    {
      public:

        //Count of edges and vertices in the graph
        VertexIdType numVertices;
        EdgeIdType numEdges;

        //contiguous adjacency list of all vertices, size = numEdges
        std::vector<VertexIdType> adjcny_in;  
        std::vector<VertexIdType> adjcny_out;  

        //offsets in adjacency list for each vertex, size = numVertices + 1
        std::vector<EdgeIdType> offsets_in;
        std::vector<EdgeIdType> offsets_out;

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
            assert(offsets_out.size() == this->numVertices + 1);

            for(auto off : offsets_in)
              assert(off >=0 && off <= this->numEdges); 

            for(auto off : offsets_out)
              assert(off >=0 && off <= this->numEdges); 

            assert(std::is_sorted(offsets_in.begin(), offsets_in.end()));
            assert(std::is_sorted(offsets_out.begin(), offsets_out.end()));

            assert(offsets_in.back() == this->numEdges);
            assert(offsets_out.back() == this->numEdges);
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

          //out-edges
          {
            adjcny_out.reserve(this->numEdges);
            offsets_out.resize(this->numVertices + 1);

            //Sort the edge vector before adding it into adjcny
            std::sort(diEdgeVector.begin(), diEdgeVector.end());

            //Compute offsets and adjcny vectors
            offsets_out[0] = 0;
            auto it_b = diEdgeVector.begin();

            for(VertexIdType i = 0; i < this->numVertices; i++)
            {
              //Range for adjacency list of vertex i
              auto it_e = std::find_if(it_b, diEdgeVector.end(), [i](std::pair<VertexIdType, VertexIdType> &e) { return e.first > i; });

              offsets_out[i+1] = std::distance(diEdgeVector.begin(), it_e);

              for(auto it = it_b; it != it_e; it++)
                adjcny_out.push_back(it->second);

              it_b = it_e;   
            }
          }

          //for in-edges  
          {
            adjcny_in.reserve(this->numEdges);
            offsets_in.resize(this->numVertices + 1);

            //reverse the edge vector: <from, to> -> <to, from>
            for(auto &e: diEdgeVector)
              e = std::make_pair(e.second, e.first);

            //Sort the edge vector before adding it into adjcny
            std::sort(diEdgeVector.begin(), diEdgeVector.end());

            //Compute offsets and adjcny vectors
            offsets_in[0] = 0;
            auto it_b = diEdgeVector.begin();

            for(VertexIdType i = 0; i < this->numVertices; i++)
            {
              //Range for adjacency list of vertex i
              auto it_e = std::find_if(it_b, diEdgeVector.end(), [i](std::pair<VertexIdType, VertexIdType> &e) { return e.first > i; });

              offsets_in[i+1] = std::distance(diEdgeVector.begin(), it_e);

              for(auto it = it_b; it != it_e; it++)
                adjcny_in.push_back(it->second);

              it_b = it_e;   
            }
          }

        }

        /**
         * @brief     print the loaded directed graph to stderr 
         * @details   Format details (assuming n = no. of vertices):
         *              Print n+1 rows in total
         *              First row: value of n and total count of edges
         *              Subsequent rows contain vertex id, OUT-neighbors of vertices and their DNA sequences,
         *              one row per vertex
         *            This function is implemented for debugging purpose
         */
        void printGraph()
        {
          std::cerr << this->numVertices << " " << this->numEdges << "\n";

          for (VertexIdType i = 0; i < this->numVertices; i++)
          {
            std::cerr << "[" << i << "] ";
            for (EdgeIdType j = offsets_out[i]; j < offsets_out[i+1]; j++)
              std::cerr << adjcny_out[j] << " ";

            std::cerr << vertex_metadata[i] << "\n";
          }
        }
    };

}

#endif
